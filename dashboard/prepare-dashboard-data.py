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

# paper -> {link, samples, projects, na_type}
papers = {}
for project in projects:
    with open("%s/bioprojects/%s/metadata/name.txt" % (
            MGS_PIPELINE_DIR, project)) as inf:
        paper_name = inf.read().strip()
        if paper_name not in papers:
            papers[paper_name] = {}
        if "projects" not in papers[paper_name]:
            papers[paper_name]["projects"] = []
        papers[paper_name]["projects"].append(project)

for paper_name in papers:
    paper_dir = os.path.join(MGS_PIPELINE_DIR, "papers",
                             paper_name.replace(" ", ""))
    link_fname = os.path.join(paper_dir, "link.txt")
    if os.path.exists(link_fname):
        with open(link_fname) as inf:
            papers[paper_name]["link"] = inf.read().strip()

    na_fname = os.path.join(paper_dir, "na_type.txt")
    if os.path.exists(na_fname):
        with open(na_fname) as inf:
            papers[paper_name]["na_type"] = inf.read().strip()

project_accession_virus_counts = {}
for project in projects:
    for accession in project_sample_reads[project]:
        with open("humanviruses/%s.humanviruses.tsv" % accession) as inf:
            for line in inf:
                line = line.strip()
                if not line: continue

                taxid, count, name = line.split("\t")
                taxid = int(taxid)
                count = int(count)

                project_accession_virus_counts[project, accession, taxid] = count

# virus -> sample -> count
virus_sample_counts = {}

# sample -> {reads, date, country, location, fine_location}
sample_metadata = {}
for project in projects:
    project_total = 0
    for accession in project_sample_reads[project]:
        if accession not in sample_metadata:
            sample_metadata[accession] = {}
            
        project_total += project_sample_reads[project][accession]
        sample_metadata[accession]["reads"] = \
            project_sample_reads[project][accession]

    with open("%s/bioprojects/%s/metadata/metadata.tsv" % (
            MGS_PIPELINE_DIR, project)) as inf:
        for line in inf:
            if not line.strip(): continue
            line = line[:-1]  # drop trailing newline

            if project in papers["Rothman 2021"]["projects"]:
                accession, date, wtp = line.split("\t")
                sample_metadata[accession]["date"] = date
                sample_metadata[accession]["country"] = "USA"
                sample_metadata[accession]["location"] = "Los Angeles"
                sample_metadata[accession]["fine_location"] = wtp
            elif project in papers["Munk 2022"]["projects"]:
                accession, country, location, date = line.split("\t")
                sample_metadata[accession]["date"] = date
                sample_metadata[accession]["country"] = country
                sample_metadata[accession]["location"] = location
            elif project in papers["Petersen 2015"]["projects"]:
                accession, country, city = line.split("\t")
                sample_metadata[accession]["country"] = country
                sample_metadata[accession]["location"] = city
                # Per Supplementary Table 7 they're all one of
                # 23-08-2013, 27-06-2013, 29-08-2013, 24-08-2013.  But the
                # mapping between samples and dates doesn't seem to be in the
                # paper.
                sample_metadata[accession]["date"] = "Summer 2013"
            elif project in papers["Maritz 2019"]["projects"]:
                sample_metadata[accession]["country"] = "USA"
                sample_metadata[accession]["location"] = "New York City"
                # Paper has "17 raw sewage samples collected from 14 DEP
                # wastewater treatment plants from the five NYC boroughs in
                # November 2014".
                sample_metadata[accession]["date"] = "2014-11"
            else:
                raise Exception("Metadata format for %s unknown" % project)

for virus_taxid in human_viruses:
    virus_sample_counts[virus_taxid] = {}
    for project in projects:
        for accession in project_sample_reads[project]:
            count = project_accession_virus_counts.get((
                project, accession, virus_taxid), 0)
            if count > 0:
                virus_sample_counts[virus_taxid][accession] = count


with open("data.js", "w") as outf:
    for name, val in [
            ("virus_sample_counts", virus_sample_counts),
            ("sample_metadata", sample_metadata),
            ("papers", papers),
            ("names", names),
            ("tree", tree),
    ]:
        outf.write("%s=%s;\n" % (name, json.dumps(
            val, sort_keys=True, indent=None if name == "tree" else 2)))
