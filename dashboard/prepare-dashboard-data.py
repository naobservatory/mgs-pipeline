#!/usr/bin/env python3

import os
import glob
import gzip
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

# taxid -> [name]
# first name is scientific name
human_virus_names = defaultdict(list)
with open("names.dmp") as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace(
            "\t|\n", "").split("\t|\t")
        taxid = int(taxid)

        if taxid in mentioned_taxids:
            if name_class == "scientific name":
                human_virus_names[taxid].insert(0, name)
            else:
                human_virus_names[taxid].append(name)

# project -> sample -> n_reads
project_sample_reads = defaultdict(dict)
for metadata_fname in glob.glob(
        "%s/bioprojects/*/metadata/metadata.tsv" % MGS_PIPELINE_DIR):
    project = metadata_fname.split("/")[-3]
    with open(metadata_fname) as inf:
        for line in inf:
            sample = line.strip().split("\t")[0]
            reads_fname = "%s/bioprojects/%s/metadata/%s.n_reads" %(
                MGS_PIPELINE_DIR, project, sample)
            if not os.path.exists(reads_fname): continue
            with open(reads_fname) as readsf:
                 project_sample_reads[project][sample] = int(readsf.read())

projects = list(sorted(project_sample_reads))

# paper -> {link, samples, projects, na_type, subset}
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

# virus -> sample -> count
virus_sample_counts = defaultdict(Counter)

# sample -> {reads, date, country, location, fine_location, na_type}
sample_metadata = defaultdict(dict)
for project in projects:
    project_total = 0
    for sample in project_sample_reads[project]:
        project_total += project_sample_reads[project][sample]
        sample_metadata[sample]["reads"] = \
            project_sample_reads[project][sample]

    with open("%s/bioprojects/%s/metadata/metadata.tsv" % (
            MGS_PIPELINE_DIR, project)) as inf:
        for line in inf:
            if not line.strip(): continue
            line = line[:-1]  # drop trailing newline

            if project in papers["Rothman 2021"]["projects"]:
                sample, date, wtp = line.split("\t")
                sample_metadata[sample]["date"] = date
                sample_metadata[sample]["country"] = "USA"
                sample_metadata[sample]["location"] = "Los Angeles"
                sample_metadata[sample]["fine_location"] = wtp
            elif project in papers["Crits-Christoph 2021"]["projects"]:
                sample, municipality, date, method, sequencing = line.split("\t")
                sample_metadata[sample]["date"] = date
                sample_metadata[sample]["country"] = "USA"
                sample_metadata[sample]["location"] = "San Francisco"
                sample_metadata[sample]["fine_location"] = municipality
            elif project in papers["Brumfield 2022"]["projects"]:
                sample, na_type, date = line.split("\t")
                sample_metadata[sample]["date"] = date
                sample_metadata[sample]["country"] = "USA"
                sample_metadata[sample]["location"] = "Maryland"
                sample_metadata[sample]["fine_location"] = "Manhole"
                sample_metadata[sample]["na_type"] = na_type
            elif project in papers["Bengtsson-Palme 2016"]["projects"]:
                sample, location, site = line.split("\t")
                sample_metadata[sample]["date"] = "2012-09"
                sample_metadata[sample]["country"] = "Sweden"
                sample_metadata[sample]["location"] = location
                sample_metadata[sample]["fine_location"] = site
            elif project in papers["Brinch 2020"]["projects"]:
                sample, loc, date = line.split("\t")
                sample_metadata[sample]["date"] = date
                sample_metadata[sample]["country"] = "Denmark"
                sample_metadata[sample]["location"] = "Copenhagen"
                sample_metadata[sample]["fine_location"] = loc
            elif project in papers["Munk 2022"]["projects"]:
                sample, country, location, date = line.split("\t")
                sample_metadata[sample]["date"] = date
                sample_metadata[sample]["country"] = country
                sample_metadata[sample]["location"] = location
            elif project in papers["Petersen 2015"]["projects"]:
                sample, country, city = line.split("\t")
                sample_metadata[sample]["country"] = country
                sample_metadata[sample]["location"] = city
                # Per Supplementary Table 7 they're all one of
                # 23-08-2013, 27-06-2013, 29-08-2013, 24-08-2013.  But the
                # mapping between samples and dates doesn't seem to be in the
                # paper.
                sample_metadata[sample]["date"] = "Summer 2013"
            elif project in papers["Maritz 2019"]["projects"]:
                sample = line.split("\t")[0]
                sample_metadata[sample]["country"] = "USA"
                sample_metadata[sample]["location"] = "New York City"
                # Paper has "17 raw sewage samples collected from 14 DEP
                # wastewater treatment plants from the five NYC boroughs in
                # November 2014".
                sample_metadata[sample]["date"] = "2014-11"
            elif project in papers["Fierer 2022"]["projects"]:
                sample = line.strip().split("\t")[0]
                sample_metadata[sample]["country"] = "USA"
                sample_metadata[sample]["location"] = "Boulder, CO"
                sample_metadata[sample]["date"] = "2020-09"
                # I can't find metadata on which samples are from which days or
                # which spots on campus.
            elif project in papers["Ng 2019"]["projects"]:
                sample, stage, date = line.strip().split("\t")
                sample_metadata[sample]["country"] = "Singapore"
                sample_metadata[sample]["location"] = "Singapore"
                sample_metadata[sample]["date"] = date
                sample_metadata[sample]["fine_location"] = {
                    "Effluent from Membrane Biorector (MBR)": "MBR",
                    "Effluent from Primary Settling Tank (PST)": "PST",
                    "Effluent from Secondary Settling Tank (SST)": "SST",
                    "Effluent from Wet Well (WW)": "WW",
                    "Influent": "Influent",
                    "Sludge (SLUDGE)": "Sludge",
                }[stage]
            elif project in papers["Hendriksen 2019"]["projects"]:
                sample, date, cluster = line.strip().split("\t")
                sample_metadata[sample]["country"] = "Kenya"
                sample_metadata[sample]["location"] = "Kibera"
                sample_metadata[sample]["fine_location"] = cluster
                sample_metadata[sample]["date"] = date
            elif project in papers["Yang 2020"]["projects"]:
                sample, city = line.strip().split("\t")
                sample_metadata[sample]["country"] = "China"
                sample_metadata[sample]["location"] = city
                sample_metadata[sample]["date"] = "2018"
            elif project in papers["Wang 2022"]["projects"]:
                sample, date, hospital = line.strip().split("\t")
                sample_metadata[sample]["country"] = "Saudi Arabia"
                sample_metadata[sample]["location"] = "Jeddah"
                sample_metadata[sample]["date"] = date
                sample_metadata[sample]["fine_location"] = hospital
            elif project in papers["Johnson 2023"]["projects"]:
                sample, _, _ = line.strip().split("\t")
                sample_metadata[sample]["country"] = "United States"
                sample_metadata[sample]["location"] = "Municipal"
                sample_metadata[sample]["date"] = "2022"
            elif project in papers["Cui 2023"]["projects"]:
                sample = line.strip()
                sample_metadata[sample]["country"] = "China"
                # Also possible this was Changchun
                sample_metadata[sample]["location"] = "Harbin"
                # Also possible this was 2020-10-15
                sample_metadata[sample]["date"] = "2022-10-19"
            else:
                raise Exception("Metadata format for %s unknown" % project)

for taxid in human_viruses:
    for project in projects:
        for sample in project_sample_reads[project]:
            count = project_sample_virus_counts[project, sample, taxid]
            if count > 0:
                virus_sample_counts[taxid][sample] = count

# make it json-serializable
for bioproject in bioprojects:
    bioprojects[bioproject] = list(sorted(bioprojects[bioproject]))

with open("data.js", "w") as outf:
    for name, val in [
            ("virus_sample_counts", virus_sample_counts),
            ("sample_metadata", sample_metadata),
            ("bioprojects", bioprojects),
            ("papers", papers),
            ("names", human_virus_names),
            ("tree", human_virus_tree),
    ]:
        outf.write("%s=%s;\n" % (name, json.dumps(
            val, sort_keys=True, indent=None if val is human_virus_tree else 2)))

for name, val in [
        ("human_virus_sample_counts", virus_sample_counts),
        ("human_virus_names", human_virus_names),
        ("human_virus_tree", human_virus_tree),
        ("metadata_samples", sample_metadata),
        ("metadata_bioprojects", bioprojects),
        ("metadata_papers", papers),
]:
    with open(name + ".json", "w") as outf:
        json.dump(val, outf, sort_keys=True,
                  indent=None if val is human_virus_tree else 2)
