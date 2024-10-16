#!/usr/bin/env python3

import os
import sys
import glob
import gzip
import json
import math
import importlib
import numpy as np
from collections import Counter
from collections import defaultdict

ROOT_DIR, MGS_PIPELINE_DIR = sys.argv[1:]
DASHBOARD_DIR = ROOT_DIR + "/dashboard/"

sys.path.insert(0, DASHBOARD_DIR)
import sample_metadata_classifier

# bioproject -> [samples]
bioprojects = defaultdict(set)

# project -> sample -> n_reads
project_sample_reads = defaultdict(Counter)
# project -> sample -> n_nonhuman_reads
project_sample_nonhuman_reads = defaultdict(Counter)

div_samples = defaultdict(set) # sample -> div_samples
div_sample_reads = {} # div_sample -> n_reads
for metadata_fname in glob.glob(
    "%s/bioprojects/*/metadata/metadata.tsv" % ROOT_DIR
):
    project = metadata_fname.split("/")[-3]
    if project in sample_metadata_classifier.skip_projects:
        # didn't finish importing this one, and the dashboard chokes on papers
        # where we don't understand the metadata.
        continue
    with open(metadata_fname) as inf:
        for line in inf:
            div_sample = line.strip().split("\t")[0]
            sample = sample_metadata_classifier.recombine(div_sample, project)
            div_samples[sample].add(div_sample)
            bioprojects[project].add(sample)

            for reads_fname, reads_dict, use_for_div_sample in [
                    ("%s/bioprojects/%s/metadata/%s.n_reads" % (
                        ROOT_DIR, project, div_sample),
                     project_sample_reads,
                     True),
                    ("%s/bioprojects/%s/metadata/%s.nonhuman.n_reads" % (
                        ROOT_DIR, project, div_sample),
                     project_sample_nonhuman_reads,
                     False)]:
                if not os.path.exists(reads_fname):
                    continue
                with open(reads_fname) as readsf:
                    reads_str = readsf.read().strip()
                    if not reads_str:
                        continue
                    n_reads = int(reads_str)

                    reads_dict[project][sample] += n_reads
                    if use_for_div_sample:
                        div_sample_reads[div_sample] = n_reads

projects = list(sorted(project_sample_reads))

papers = {}
for project in projects:
    name_fname = "%s/bioprojects/%s/metadata/name.txt" % (ROOT_DIR, project)
    if os.path.exists(name_fname):
        with open(name_fname) as inf:
            paper_name = inf.read().strip()
    else:
        paper_name = project

    if paper_name not in papers:
        papers[paper_name] = {}
    if "projects" not in papers[paper_name]:
        papers[paper_name]["projects"] = []
    papers[paper_name]["projects"].append(project)

# sample -> {metadata}
sample_metadata = defaultdict(dict)

papers_dir = os.path.join(ROOT_DIR, "papers")
sys.path.insert(0, papers_dir)
for paper_name in papers:
    paper_metadata_name = "%s--metadata" % paper_name.replace(" ", "")
    paper_metadata_fname = paper_metadata_name + ".py"

    if os.path.exists(os.path.join(papers_dir, paper_metadata_fname)):
        paper_metadata_module = importlib.import_module(
            paper_metadata_name)
        papers[paper_name].update(paper_metadata_module.paper_metadata())
    else:
        paper_metadata_module = None
        paper_dir = os.path.join(papers_dir, paper_name.replace(" ", ""))
        for metadata_type in ["link",
                              "na_type",
                              "subset",
                              "mgs-workflow-output"]:
            fname = os.path.join(paper_dir, "%s.txt" % metadata_type)
            if os.path.exists(fname):
                with open(fname) as inf:
                    papers[paper_name][metadata_type] = inf.read().strip()

    if "link" not in papers[paper_name]:
        papers[paper_name]["link"] = "personal communication"

    for project in papers[paper_name]["projects"]:
        with open(
                "%s/bioprojects/%s/metadata/metadata.tsv" % (ROOT_DIR, project)
        ) as inf:
            for line in inf:
                if not line.strip():
                    continue
                bits = line.rstrip("\n").split("\t")

                if paper_metadata_module:
                    div_sample, sample_metadata_dict = \
                        paper_metadata_module.sample_metadata(bits)
                else:
                    if not getattr(sample_metadata_classifier, "interpret"):
                        raise Exception(
                            "sample_metadata_classifier.interpret() not "
                            "found. Did you forget to create %s?" %
                            paper_metadata_fname)
                            
                    div_sample, sample_metadata_dict = \
                        sample_metadata_classifier.interpret(
                            project, papers, bits)

                sample = sample_metadata_classifier.recombine(div_sample, project)

                if sample not in sample_metadata:
                    sample_metadata[sample] = sample_metadata_dict
                elif sample_metadata_dict != sample_metadata[sample]:
                    print(sample, div_sample)
                    import pprint
                    pprint.pprint(sample_metadata_dict)
                    pprint.pprint(sample_metadata[sample])
                    assert False

        for sample in project_sample_reads[project]:
            sample_metadata[sample]["reads"] = \
                project_sample_reads[project][sample]
            if project_sample_nonhuman_reads[project][sample]:
                sample_metadata[sample]["nonhuman_reads"] = \
                    project_sample_nonhuman_reads[project][sample]

# make it json-serializable
for bioproject in bioprojects:
    bioprojects[bioproject] = list(sorted(bioprojects[bioproject]))

for name, val in [
    ("metadata_samples", sample_metadata),
    ("metadata_bioprojects", bioprojects),
    ("metadata_papers", papers),
]:
    with open(DASHBOARD_DIR + name + ".json", "w") as outf:
        json.dump(
            val,
            outf,
            sort_keys=True,
            indent=2,
        )

SKIP_PAPER_KEYS = set(("projects", "summary", "link", "mgs-workflow-output"))
for paper in papers:
    for bioproject in papers[paper]["projects"]:
        bioproject_samples = bioprojects[bioproject]
        with open(DASHBOARD_DIR + "%s.metadata.tsv" % bioproject, "w") as outf:

            keys = set()
            keys.update([key for key in papers[paper]
                         if key not in SKIP_PAPER_KEYS])
            for sample in bioproject_samples:
                keys.update(sample_metadata[sample])

            header = ["sample", *sorted(keys)]
            outf.write("\t".join(header) + "\n")

            for sample in sorted(bioproject_samples):
                outf.write("%s\t%s\n" % (
                    sample, "\t".join(
                        # pull the key from the sample-level metadata if we
                        # have it, otherwise fall back to paper-level.  And use
                        # an empty string if it's unset.
                        str(sample_metadata[sample].get(
                            key, papers[paper].get(key, "")))
                        for key in header[1:])))
