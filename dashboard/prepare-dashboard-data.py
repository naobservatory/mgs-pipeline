#!/usr/bin/env python3

import os
import glob
import json
from collections import Counter
from collections import defaultdict

"""
Generate data for interactive HTML for exploring the taxonomy of human viruses
identified in the samples.

Columns are bioprojects, rows are taxonomic nodes, cells are relative abundances.

Rows can be expanded or collapsed to show child nodes.
"""

THISDIR=os.path.dirname(__file__)
MGS_PIPELINE_DIR=os.path.abspath(os.path.join(THISDIR, ".."))

human_viruses = set()   # {taxid}
with open("%s/human-viruses.tsv" % MGS_PIPELINE_DIR) as inf:
    for line in inf:
        taxid, name = line.strip().split("\t")
        taxid = int(taxid)
        human_viruses.add(taxid)

parents = {}  # child_taxid -> parent_taxid
with open("nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)
        parents[child_taxid] = parent_taxid

# taxid -> node
nodes = {}

mentioned_taxids = set()
for virus_taxid in human_viruses:
    taxid = virus_taxid
    mentioned_taxids.add(taxid)
    while taxid != 1:
        taxid = parents[taxid]
        mentioned_taxids.add(taxid)

for taxid in mentioned_taxids:
    nodes[taxid] = [taxid]

for virus_taxid in human_viruses:
    taxid = virus_taxid
    node = nodes[taxid]
    while taxid != 1:
        taxid = parents[taxid]
        if node in nodes[taxid]:
            break
        nodes[taxid].append(node)
        node = nodes[taxid]

# Format: [taxid, children]
# Ex: [1, [10239, [10472,
#       [10473, [46014],
#               [10474, [10475, ...],
#                       [693626, ...],
#                       [1299307, ...]]], ...], ...], ...]
tree = nodes[1]

# taxid -> [name]
# first name is scientific name
names = {}
with open("names.dmp") as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace(
            "\t|\n", "").split("\t|\t")
        taxid = int(taxid)
        if taxid in mentioned_taxids:
            if taxid not in names:
                names[taxid] = []
            if name_class == "scientific name":
                names[taxid].insert(0, name)
            else:
                names[taxid].append(name)

# project -> accession -> n_reads
project_sample_reads = defaultdict(dict)
for n_reads_fname in glob.glob(
        "%s/bioprojects/*/metadata/*.n_reads" % MGS_PIPELINE_DIR):
    project = n_reads_fname.split("/")[-3]
    accession = n_reads_fname.split("/")[-1].replace(".n_reads", "")
    with open(n_reads_fname) as inf:
        project_sample_reads[project][accession] = int(inf.read())

projects = list(sorted(project_sample_reads))

bioproject_names = {} # project -> bioproject name
for project in projects:
    with open("%s/bioprojects/%s/metadata/name.txt" % (
            MGS_PIPELINE_DIR, project)) as inf:
        bioproject_names[project] = inf.read().strip()

bioproject_links = {} # project -> paper link
for project, name in bioproject_names.items():
    paper_dir = os.path.join(MGS_PIPELINE_DIR, "papers", name.replace(" ", ""))
    link_fname = os.path.join(paper_dir, "link.txt")
    if os.path.exists(link_fname):
        with open(link_fname) as inf:
            bioproject_links[project] = inf.read().strip()

bioproject_totals = {}

project_accession_virus_counts = {}
for project in projects:
    for accession in project_sample_reads[project]:
        with open("%s.humanviruses.tsv" % accession) as inf:
            for line in inf:
                line = line.strip()
                if not line: continue
            
                taxid, count, name = line.split("\t")
                taxid = int(taxid)
                count = int(count)

                project_accession_virus_counts[project, accession, taxid] = count

# virus -> project -> [count, relab]
virus_project_counts = {}
for virus_taxid in human_viruses:
    virus_project_counts[virus_taxid] = {}
    for project in projects:
        project_count = 0
        project_total = 0
        for accession in project_sample_reads[project]:
            project_total += project_sample_reads[project][accession]
            project_count += project_accession_virus_counts.get((
                project, accession, virus_taxid), 0)

        if project_count > 0:
            virus_project_counts[virus_taxid][project] = (
                project_count, project_count / project_total)
        bioproject_totals[project] = project_total

sorted_projects = [project for key, project in
                   sorted(((bioproject_names[project].split()[-1],
                            bioproject_names[project]),
                           project)
                          for project in projects)]

with open("data.js", "w") as outf:
    for name, val in [
            ("virus_project_counts", virus_project_counts),
            ("projects", sorted_projects),
            ("names", names),
            ("bioproject_names", bioproject_names),
            ("bioproject_links", bioproject_links),
            ("bioproject_totals", bioproject_totals),
            ("tree", tree),
    ]:
        outf.write("%s=%s;\n" % (name, json.dumps(
            val, sort_keys=True, indent=None if name == "tree" else 2)))
