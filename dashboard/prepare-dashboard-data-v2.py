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

"""
Generate data for interactive HTML for exploring the taxonomy of human viruses
identified in the samples.

Columns are bioprojects, rows are taxonomic nodes, cells are relative abundances.

Rows can be expanded or collapsed to show child nodes.

We include a row for a taxid under human viruses if any of:
1. It is a human virus and was observed in any study
2. It is a human virus and it's parent isn't
3. It is an ancestor of something we include
"""

ROOT_DIR, MGS_PIPELINE_DIR = sys.argv[1:]
DASHBOARD_DIR = ROOT_DIR + "/dashboard/"

sys.path.insert(0, DASHBOARD_DIR)
import sample_metadata_classifier

# bioproject -> [samples]
bioprojects = defaultdict(set)

all_human_viruses = set()  # {taxid}
with open("%s/human-viruses.tsv" % MGS_PIPELINE_DIR) as inf:
    for line in inf:
        taxid, name = line.strip().split("\t")
        taxid = int(taxid)
        all_human_viruses.add(taxid)

parents = {}  # child_taxid -> parent_taxid
with open("%s/nodes.dmp" % DASHBOARD_DIR) as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = line.replace("\t|\n", "").split(
            "\t|\t"
        )
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)
        parents[child_taxid] = parent_taxid

for taxid in all_human_viruses:
    if taxid not in parents:
        print("Missing %s" % taxid)

root_human_viruses = set()
for taxid in all_human_viruses:
    if parents[taxid] not in all_human_viruses:
        root_human_viruses.add(taxid)

# project -> sample
project_samples = defaultdict(set)
div_samples = defaultdict(set) # sample -> div_samples
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
            bioprojects[project].add(sample)
            div_samples[sample].add(div_sample)
            project_samples[project].add(sample)

projects = list(sorted(project_samples))

observed_taxids = set()
for fname in glob.glob("hv_hits_putative_collapsed/*.tsv.gz"):
    with gzip.open(fname, "rt") as inf:
        col = None
        for line in inf:
            bits = line.rstrip("\n").split("\t")
            if not col:
                col = bits
                continue

            taxid = int(bits[col.index("taxid")])
            if taxid in all_human_viruses:
                observed_taxids.add(taxid)

# taxid -> node
human_virus_nodes = {}

mentioned_taxids = set()
for virus_taxid in root_human_viruses | observed_taxids:
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

def rc(s):
    return "".join(
        {"T": "A", "G": "C", "A": "T", "C": "G", "N": "N"}[x]
        for x in reversed(s)
    )

project_sample_virus_counts = Counter()
for project in projects:
    fname = "hv_hits_putative_collapsed/%s.tsv.gz" % project
    if not os.path.exists(fname):
        continue

    with gzip.open(fname, "rt") as inf:
        col = None
        for line in inf:
            bits = line.rstrip("\n").split("\t")
            if not col:
                col = bits
                continue

            taxid = int(bits[col.index("taxid")])
            sample = bits[col.index("sample")]
            project_sample_virus_counts[project, sample, taxid] += 1

# comparison taxid -> sample -> clade count
comparison_sample_counts = defaultdict(Counter)
for project in projects:
    fname = "top_species_counts_v2/%s.json" % project
    if not os.path.exists(fname):
        continue
    with open(fname) as inf:
        for sample, taxid_counts in json.load(inf).items():
            for taxid, count in taxid_counts.items():
                comparison_sample_counts[int(taxid)][sample] += count

BACTERIA = 2
VIRUS = 10239
comparison_taxid_classifications = {
    BACTERIA: [],
    VIRUS: [],
}
for taxid in sorted(comparison_sample_counts):
    seen_taxids = []
    p = taxid
    while p not in [1, 0]:
        seen_taxids.append(p)
        if p in comparison_taxid_classifications:
            comparison_taxid_classifications[p].append(taxid)
            break
        if p not in parents:
            raise Exception("Missing Parent: %s" %
                            ",".join(str(x) for x in seen_taxids))
        p = parents[p]

# taxid -> [name]
# first name is scientific name
taxonomic_names = defaultdict(list)
with open("%s/names.dmp" % DASHBOARD_DIR) as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace("\t|\n", "").split(
            "\t|\t"
        )
        taxid = int(taxid)

        if taxid in mentioned_taxids or taxid in comparison_sample_counts:
            if name_class == "scientific name":
                taxonomic_names[taxid].insert(0, name)
            else:
                taxonomic_names[taxid].append(name)

# virus -> sample -> count
virus_sample_counts = defaultdict(Counter)

for taxid in observed_taxids:
    for project in projects:
        for sample in project_samples[project]:
            count = project_sample_virus_counts[project, sample, taxid]
            if count > 0:
                virus_sample_counts[taxid][sample] = count

for name, val in [
    ("human_virus_sample_counts", virus_sample_counts),
    ("taxonomic_names", taxonomic_names),
    ("human_virus_tree", human_virus_tree),
    ("comparison_taxid_classifications", comparison_taxid_classifications),
]:
    if not val:
        continue

    with open(DASHBOARD_DIR + name + "_v2.json", "w") as outf:
        json.dump(
            val,
            outf,
            sort_keys=True,
            indent=None if val is human_virus_tree else 2,
        )

# To make the dashboard load faster, divide counts by bioproject and don't load
# them by default.
for bioproject in bioprojects:
    samples = set(bioprojects[bioproject])

    for full_dict, name in [
        (virus_sample_counts, "human_virus_sample_counts.%s" % bioproject),
        (comparison_sample_counts, "comparison_sample_counts.%s" % bioproject),
    ]:
        bioproject_sample_counts = {}
        for taxid, sample_counts in full_dict.items():
            counts = {
                sample: sample_counts[sample]
                for sample in samples & sample_counts.keys()
            }
            if counts:
                bioproject_sample_counts[taxid] = counts

        if not bioproject_sample_counts:
            continue

        with open(DASHBOARD_DIR + name + "_v2.json", "w") as outf:
            json.dump(bioproject_sample_counts, outf, sort_keys=True, indent=2)
