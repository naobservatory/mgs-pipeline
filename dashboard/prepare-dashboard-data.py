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
with open("names.dmp") as inf:
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
            MGS_PIPELINE_DIR, project)) as inf:
        for line in inf:
            if not line.strip(): continue
            line = line[:-1]  # drop trailing newline

            if project in papers["Rothman 2021"]["projects"]:
                sample, date, wtp, is_enriched = line.split("\t")
                if wtp == "JW":
                    # Rothman confirmed over email that JW = JWPCP.
                    wtp = "JWPCP"

                sample_metadata[sample] = dict(
                    date=date,
                    country="USA",
                    location="Los Angeles",
                    county={
                        # Hyperion
                        "HTP": "Los Angeles County",
                        # San Jose Creek
                        "SJ": "Los Angeles County",
                        # Joint Water Pollution Control Plant
                        "JWPCP": "Los Angeles County",
                        # Orange County
                        "OC": "Orange County",
                        # Point Loma
                        "PL": "San Diego County",
                        # South Bay
                        "SB": "San Diego County",
                        # North City
                        "NC": "San Diego County",
                        # Escondido Hale Avenue Resource Recovery Facility
                        "ESC": "San Diego County",
                    }[wtp],
                    fine_location=wtp,
                    enrichment="panel" if is_enriched == "1" else "viral")
            elif project in papers["Crits-Christoph 2021"]["projects"]:
                sample, municipality, date, method, sequencing = line.split("\t")
                sample_metadata[sample] = dict(
                    date=date,
                    country="USA",
                    location="San Francisco",
                    county={
                        "Berkeley": "Alameda County",
                        "Marin": "Marin County",
                        "Oakland": "Alameda County",
                        "SF": "San Francisco County",
                    }[municipality],
                    fine_location=municipality,
                    method=method,
                    enrichment="panel" if sequencing == "enriched" else "viral")
            elif project in papers["Brumfield 2022"]["projects"]:
                sample, na_type, date = line.split("\t")
                sample_metadata[sample] = dict(
                    date=date,
                    country="USA",
                    location="Maryland",
                    fine_location="Manhole",
                    na_type=na_type)
            elif project in papers["Bengtsson-Palme 2016"]["projects"]:
                sample, location, site = line.split("\t")
                sample_metadata[sample] = dict(
                    date="2012-09",
                    country="Sweden",
                    location=location,
                    fine_location=site)
            elif project in papers["Brinch 2020"]["projects"]:
                sample, loc, date = line.split("\t")
                sample_metadata[sample] = dict(
                    date=date,
                    country="Denmark",
                    location="Copenhagen",
                    fine_location=loc)
            elif project in papers["Spurbeck 2023"]["projects"]:
                sample, loc, date = line.split("\t")
                sample_metadata[sample] = dict(
                    date=date,
                    country="USA",
                    location="Ohio",
                    # https://github.com/naobservatory/mgs-pipeline/issues/9
                    county={
                        "A": "Summit County",
                        "B": "Trumbull County",
                        "C": "Lucas County",
                        "D": "Lawrence County",
                        "E": "Sandusky County",
                        "F": "Franklin County",
                        "G": "Licking County",
                        "H": "Franklin County",
                        "I": "Greene County",
                        "J": "Montogmery County",
                    }[loc],
                    fine_location=loc,
                    enrichment="viral",
                    method={
                        "A": "AB",
                        "B": "AB",
                        "C": "C",
                        "D": "D",
                        "E": "EFGH",
                        "F": "EFGH",
                        "G": "EFGH",
                        "H": "EFGH",
                        "I": "IJ",
                        "J": "IJ",
                    }[loc])

            elif project in papers["Munk 2022"]["projects"]:
                sample, country, location, date = line.split("\t")
                sample_metadata[sample] = dict(
                    date=date,
                    country=country,
                    location=location)
            elif project in papers["Petersen 2015"]["projects"]:
                sample, country, city = line.split("\t")
                sample_metadata[sample] = dict(
                    country=country,
                    location=city,
                    # Per Supplementary Table 7 they're all one of 23-08-2013,
                    # 27-06-2013, 29-08-2013, 24-08-2013.  But the mapping
                    # between samples and dates doesn't seem to be in the
                    # paper.
                    date="Summer 2013")
            elif project in papers["Maritz 2019"]["projects"]:
                sample = line.split("\t")[0]
                sample_metadata[sample] = dict(
                    country="USA",
                    location="New York City",
                    # Paper has "17 raw sewage samples collected from 14 DEP
                    # wastewater treatment plants from the five NYC boroughs in
                    # November 2014".
                    date="2014-11")
            elif project in papers["Fierer 2022"]["projects"]:
                sample = line.strip().split("\t")[0]
                sample_metadata[sample] = dict(
                    country="USA",
                    location="Boulder, CO",
                    # I can't find metadata on which samples are from which
                    # days or which spots on campus.
                    date="2020-09")
            elif project in papers["Ng 2019"]["projects"]:
                sample, stage, date = line.strip().split("\t")
                sample_metadata[sample] = dict(
                    country="Singapore",
                    location="Singapore",
                    date=date,
                    fine_location={
                        "Effluent from Membrane Biorector (MBR)": "MBR",
                        "Effluent from Primary Settling Tank (PST)": "PST",
                        "Effluent from Secondary Settling Tank (SST)": "SST",
                        "Effluent from Wet Well (WW)": "WW",
                        "Influent": "Influent",
                        "Sludge (SLUDGE)": "Sludge",
                    }[stage])
            elif project in papers["Hendriksen 2019"]["projects"]:
                sample, date, cluster = line.strip().split("\t")
                sample_metadata[sample] = dict(
                    country="Kenya",
                    location="Kibera",
                    fine_location=cluster,
                    date=date)
            elif project in papers["Yang 2020"]["projects"]:
                sample, city = line.strip().split("\t")
                sample_metadata[sample] = dict(
                    country="China",
                    location=city,
                    date="2018",
                    enrichment="viral")
            elif project in papers["Wang 2022"]["projects"]:
                sample, date, hospital = line.strip().split("\t")
                sample_metadata[sample] = dict(
                    country="Saudi Arabia",
                    location="Jeddah",
                    date=date,
                    fine_location=hospital)
            elif project in papers["Johnson 2023"]["projects"]:
                sample, _, _ = line.strip().split("\t")
                sample_metadata[sample] = dict(
                    country="United States",
                    location="Municipal",
                    date="2022",
                    enrichment="viral")
            elif project in papers["Cui 2023"]["projects"]:
                sample = line.strip()
                sample_metadata[sample] = dict(
                    country="China",
                    # Also possible this was Changchun
                    location="Harbin",
                    # Also possible this was 2020-10-15
                    date="2022-10-19")
            else:
                raise Exception("Metadata format for %s unknown" % project)
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

with open("data.js", "w") as outf:
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
    with open(name + ".json", "w") as outf:
        json.dump(val, outf, sort_keys=True,
                  indent=None if val is human_virus_tree else 2)
 
