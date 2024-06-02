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
div_samples = defaultdict(set) # sample -> div_samples
div_sample_reads = {} # div_sample -> n_reads
for metadata_fname in glob.glob(
    "%s/bioprojects/*/metadata/metadata.tsv" % ROOT_DIR
):
    project = metadata_fname.split("/")[-3]
    if project in ["PRJEB30546", "PRJNA691135"]:
        # didn't finish importing this one, and the dashboard chokes on papers
        # where we don't understand the metadata.
        continue
    with open(metadata_fname) as inf:
        for line in inf:
            div_sample = line.strip().split("\t")[0]
            sample = sample_metadata_classifier.recombine(div_sample, project)
            div_samples[sample].add(div_sample)
            bioprojects[project].add(sample)

            reads_fname = "%s/bioprojects/%s/metadata/%s.n_reads" % (
                ROOT_DIR,
                project,
                div_sample,
            )
            if not os.path.exists(reads_fname):
                continue
            with open(reads_fname) as readsf:
                reads_str = readsf.read().strip()
                if not reads_str:
                    continue
                n_reads = int(reads_str)

                project_sample_reads[project][sample] += n_reads
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

for paper_name in papers:
    paper_dir = os.path.join(ROOT_DIR, "papers", paper_name.replace(" ", ""))
    link_fname = os.path.join(paper_dir, "link.txt")
    if os.path.exists(link_fname):
        with open(link_fname) as inf:
            papers[paper_name]["link"] = inf.read().strip()
    else:
        papers[paper_name]["link"] = "personal communication"

    na_fname = os.path.join(paper_dir, "na_type.txt")
    if os.path.exists(na_fname):
        with open(na_fname) as inf:
            papers[paper_name]["na_type"] = inf.read().strip()

    subset_fname = os.path.join(paper_dir, "subset.txt")
    if os.path.exists(subset_fname):
        with open(subset_fname) as inf:
            papers[paper_name]["subset"] = inf.read().strip()

# sample -> {metadata}
sample_metadata = defaultdict(dict)

def summarize_readlength_category(rls_category):
    weighted_sum = 0
    weighted_total = 0
    nc = 0
    for rl_category in rls_category:
        for value, count in rl_category.items():
            if value == "NC":
                nc += count
            else:
                weighted_sum += int(value) * count
                weighted_total += count

    if not weighted_total:
        return -1, -1

    if weighted_sum / weighted_total < 5:
        import pprint
        pprint.pprint(rls_category)
        print(weighted_sum)
        print(weighted_total)
        exit(1)

    return [nc / (weighted_total + nc),
            weighted_sum / weighted_total]

def summarize_readlengths(rls):
    all_categories = set()
    for rl in rls:
        for category in rl:
            all_categories.add(category)

    return {category: summarize_readlength_category([
        rl[category] for rl in rls])
            for category in all_categories}

for project in projects:
    with open(
        "%s/bioprojects/%s/metadata/metadata.tsv" % (ROOT_DIR, project)
    ) as inf:
        for line in inf:
            if not line.strip():
                continue
            line = line[:-1]  # drop trailing newline

            (
                div_sample,
                sample_metadata_dict,
            ) = sample_metadata_classifier.interpret(
                project, papers, line.split("\t")
            )
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
        sample_metadata[sample]["reads"] = project_sample_reads[project][
            sample
        ]
        rls = []
        for div_sample in div_samples[sample]:
            rl_fname = "readlengths/%s.rl.json.gz" % div_sample
            try:
                with gzip.open(rl_fname, "rt") as inf:
                    rls.append(json.load(inf))
            except FileNotFoundError:
                continue

            sample_metadata[sample]["readlengths"] = summarize_readlengths(rls)

        ribofracs = []
        weights = []
        for div_sample in div_samples[sample]:
            rf_fname = "ribofrac/%s.ribofrac.txt" % div_sample
            try:
                with open(rf_fname, "r") as file:
                    ribofrac = file.readline()
            except FileNotFoundError:
                continue
            ribofracs.append(float(ribofrac))
            weights.append(div_sample_reads[div_sample])

        if sum(weights):
            ribofrac = np.average(ribofracs, weights=weights)
            if not math.isnan(ribofrac):
                sample_metadata[sample]["ribofrac"] = ribofrac


# make it json-serializable
for bioproject in bioprojects:
    bioprojects[bioproject] = list(sorted(bioprojects[bioproject]))

def round_floats_recursively(val, precision=4):
    if type(val) == type(0.0):
        return round(val, precision)
    if type(val) in [type({}), type(defaultdict())]:
        return {k: round_floats_recursively(v) for k, v in val.items()}
    if type(val) in [type([]), type(())]:
        return [round_floats_recursively(v) for v in val]

    return val

for name, val in [
    ("metadata_samples", sample_metadata),
    ("metadata_bioprojects", bioprojects),
    ("metadata_papers", papers),
]:
    with open(DASHBOARD_DIR + name + ".json", "w") as outf:
        json.dump(
            round_floats_recursively(val),
            outf,
            sort_keys=True,
            indent=2,
        )
