#!/usr/bin/env python3

import os
import sys
import glob
import gzip
import json
import importlib
from collections import Counter
from collections import defaultdict

"""
Generate data for interactive HTML for exploring the taxonomy of human viruses
identified in the samples.

Columns are bioprojects, rows are taxonomic nodes, cells are relative abundances.

Rows can be expanded or collapsed to show child nodes.
"""

ROOT_DIR, MGS_PIPELINE_DIR = sys.argv[1:]
DASHBOARD_DIR = ROOT_DIR + "/dashboard/"

sys.path.insert(0, DASHBOARD_DIR)
import sample_metadata_classifier

human_viruses = set()   # {taxid}
with open("%s/human-viruses.tsv" % MGS_PIPELINE_DIR) as inf:
    for line in inf:
        taxid, name = line.strip().split("\t")
        taxid = int(taxid)
        human_viruses.add(taxid)

parents = {}  # child_taxid -> parent_taxid
with open("%s/nodes.dmp" % DASHBOARD_DIR) as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)
        parents[child_taxid] = parent_taxid

# taxid -> node
human_virus_nodes = {}

mentioned_taxids = set()
for virus_taxid in human_viruses:
    taxid = virus_taxid
    mentioned_taxids.add(taxid)
    while taxid != 1:
        taxid = parents[taxid]
        mentioned_taxids.add(taxid)

for taxid in mentioned_taxids:
    human_virus_nodes[taxid] = [taxid]

for taxid in sorted(mentioned_taxids):
    node = human_virus_nodes[taxid]
    while taxid != 1:
        taxid = parents[taxid]
        if node in human_virus_nodes[taxid]:
            break
        human_virus_nodes[taxid].append(node)
        node = human_virus_nodes[taxid]

# Format: [taxid, children]
# Ex: [1, [10239, [10472,
#       [10473, [46014],
#               [10474, [10475, ...],
#                       [693626, ...],
#                       [1299307, ...]]], ...], ...], ...]
human_virus_tree = human_virus_nodes[1]

# project -> sample -> n_reads
project_sample_reads = defaultdict(dict)
for metadata_fname in glob.glob(
        "%s/bioprojects/*/metadata/metadata.tsv" % ROOT_DIR):
    project = metadata_fname.split("/")[-3]
    if project == "PRJEB30546":
         # didn't finish importing this one, and the dashboard chokes on papers
         # where we don't understand the metadata.
        continue
    with open(metadata_fname) as inf:
        for line in inf:
            sample = line.strip().split("\t")[0]
            reads_fname = "%s/bioprojects/%s/metadata/%s.n_reads" %(
                ROOT_DIR, project, sample)
            if not os.path.exists(reads_fname): continue
            with open(reads_fname) as readsf:
                 project_sample_reads[project][sample] = int(readsf.read())

projects = list(sorted(project_sample_reads))

# paper -> {link, samples, projects, na_type, subset}
papers = {}
for project in projects:
    with open("%s/bioprojects/%s/metadata/name.txt" % (
            ROOT_DIR, project)) as inf:
        paper_name = inf.read().strip()
        if paper_name not in papers:
            papers[paper_name] = {}
        if "projects" not in papers[paper_name]:
            papers[paper_name]["projects"] = []
        papers[paper_name]["projects"].append(project)

for paper_name in papers:
    paper_dir = os.path.join(ROOT_DIR, "papers",
                             paper_name.replace(" ", ""))
    link_fname = os.path.join(paper_dir, "link.txt")
    if os.path.exists(link_fname):
        with open(link_fname) as inf:
            papers[paper_name]["link"] = inf.read().strip()

    na_fname = os.path.join(paper_dir, "na_type.txt")
    if os.path.exists(na_fname):
        with open(na_fname) as inf:
            papers[paper_name]["na_type"] = inf.read().strip()

    subset_fname = os.path.join(paper_dir, "subset.txt")
    if os.path.exists(subset_fname):
        with open(subset_fname) as inf:
            papers[paper_name]["subset"] = inf.read().strip()


# bioproject -> [samples]
bioprojects = defaultdict(set)

project_sample_virus_counts = Counter()
for project in projects:
    for sample in project_sample_reads[project]:
        fname = "humanviruses/%s.humanviruses.tsv" % sample
        if not os.path.exists(fname): continue

        bioprojects[project].add(sample)
        with open(fname) as inf:
            for line in inf:
                line = line.strip()
                if not line: continue

                taxid, count, name = line.split("\t")
                taxid = int(taxid)
                count = int(count)

                project_sample_virus_counts[project, sample, taxid] = count

# comparison taxid -> sample -> clade count
comparison_sample_counts = defaultdict(Counter)
for project in projects:
    for sample in project_sample_reads[project]:
        fname = "top_species_counts/%s.json" % sample
        if not os.path.exists(fname): continue
        with open(fname) as inf:
            comparisons = json.load(inf)
            for taxid, count in comparisons.items():
                comparison_sample_counts[int(taxid)][sample] = count

BACTERIA=2
VIRUS=10239
comparison_taxid_classifications = {
    BACTERIA: [],
    VIRUS: [],
}
for taxid in sorted(comparison_sample_counts):
    p = taxid
    while p not in [1, 0]:
        if p in comparison_taxid_classifications:
            comparison_taxid_classifications[p].append(taxid)
            break
        p = parents[p]

# taxid -> [name]
# first name is scientific name
taxonomic_names = defaultdict(list)
with open("%s/names.dmp" % DASHBOARD_DIR) as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace(
            "\t|\n", "").split("\t|\t")
        taxid = int(taxid)

        if taxid in mentioned_taxids or taxid in comparison_sample_counts:
            if name_class == "scientific name":
                taxonomic_names[taxid].insert(0, name)
            else:
                taxonomic_names[taxid].append(name)

# virus -> sample -> count
virus_sample_counts = defaultdict(Counter)

# sample -> {metadata}
sample_metadata = defaultdict(dict)
for project in projects:
    with open("%s/bioprojects/%s/metadata/metadata.tsv" % (
            ROOT_DIR, project)) as inf:
        for line in inf:
            if not line.strip(): continue
            line = line[:-1]  # drop trailing newline

            sample, sample_metadata_dict = sample_metadata_classifier.interpret(
                project, papers, line.split("\t"))
            sample_metadata[sample] = sample_metadata_dict

    for sample in project_sample_reads[project]:
        sample_metadata[sample]["reads"] = \
            project_sample_reads[project][sample]

for taxid in human_viruses:
    for project in projects:
        for sample in project_sample_reads[project]:
            count = project_sample_virus_counts[project, sample, taxid]
            if count > 0:
                virus_sample_counts[taxid][sample] = count

# make it json-serializable
for bioproject in bioprojects:
    bioprojects[bioproject] = list(sorted(bioprojects[bioproject]))

with open(DASHBOARD_DIR + "data.js", "w") as outf:
    for name, val in [
            ("virus_sample_counts", virus_sample_counts),
            ("comparison_sample_counts", comparison_sample_counts),
            ("comparison_taxid_classifications",
             comparison_taxid_classifications),
            ("sample_metadata", sample_metadata),
            ("bioprojects", bioprojects),
            ("papers", papers),
            ("names", taxonomic_names),
            ("tree", human_virus_tree),
    ]:
        outf.write("%s=%s;\n" % (name, json.dumps(
            val, sort_keys=True, indent=None if val is human_virus_tree else 2)))

for name, val in [
        ("human_virus_sample_counts", virus_sample_counts),
        ("taxonomic_names", taxonomic_names),
        ("human_virus_tree", human_virus_tree),
        ("comparison_sample_counts", comparison_sample_counts),
        ("comparison_taxid_classifications", comparison_taxid_classifications),
        ("metadata_samples", sample_metadata),
        ("metadata_bioprojects", bioprojects),
        ("metadata_papers", papers),
]:
    with open(DASHBOARD_DIR + name + ".json", "w") as outf:
        json.dump(val, outf, sort_keys=True,
                  indent=None if val is human_virus_tree else 2)
