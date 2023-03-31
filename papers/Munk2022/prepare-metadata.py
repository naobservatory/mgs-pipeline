#!/usr/bin/env python3

import re
import os
import glob
from collections import defaultdict

THISDIR = os.path.dirname(__file__)
root = os.path.join(THISDIR, "..", "..")

records = {}
with open("monk-2022-supplementary-data-1.tsv") as inf:
    for n, line in enumerate(inf):
        cols = line.strip().split("\t")
        if n == 0:
            colnames = cols
        else:
            records[cols[colnames.index("ena_run_acc")]] = (
                cols[colnames.index("country")],
                cols[colnames.index("city")],
                cols[colnames.index("collection_date")])

sample_alias_to_run_accessions = defaultdict(set)
with open("Hendriksen2019-metadata2.tsv") as inf:
    for line in inf:
        if line.startswith("sample_accession"): continue

        _, run_accession, sample_alias = line.strip().split()
        sample_alias = re.sub("_[0-9]+$", "", sample_alias).replace(".", "-")

        sample_alias_to_run_accessions[sample_alias].add(run_accession)

with open("Hendriksen2019-metadata3.tsv") as inf:
    for line in inf:
        if line.startswith("sample_accession"): continue

        _, run_accession, sample_alias = line.strip().split()
        sample_alias = re.sub("_[0-9]+$", "", sample_alias).replace(".", "-")

        sample_alias_to_run_accessions[sample_alias].add(run_accession)

with open("Hendriksen2019-metadata1.tsv") as inf:
    for n, line in enumerate(inf):
         cols = line.strip().split("\t")
         if n == 0:
             colnames = cols
         else:
            sample_alias = cols[colnames.index("sample_ID")].replace(
                ".", "-")
            run_accessions = sample_alias_to_run_accessions[sample_alias]
            if not run_accessions:
                raise Exception("Missing data for %s" % sample_alias)

            for run_accession in run_accessions:
                record = (cols[colnames.index("country")],
                          cols[colnames.index("city")],
                          cols[colnames.index("sample.date")])
                
                if run_accession in records:
                    continue
                records[run_accession] = record

bioproject_dir = os.path.join(root, "bioprojects")
for bioproject in os.listdir(bioproject_dir):
    with open(os.path.join(
        bioproject_dir, bioproject, "metadata", "name.txt")) as inf:
        if inf.read().strip() != "Munk 2022": continue

    metadata_fname = os.path.join(
        bioproject_dir, bioproject, "metadata", "metadata.tsv")
    out = []
    with open(metadata_fname) as inf:
        for line in inf:
            ena_run_acc = line.strip()
            if ena_run_acc not in records:
                print ("dropping", bioproject, ena_run_acc)
                continue

            out.append([ena_run_acc, *records[ena_run_acc]])

    with open(metadata_fname, "w") as outf:
        for row in out:
            outf.write("\t".join(row) + "\n")
