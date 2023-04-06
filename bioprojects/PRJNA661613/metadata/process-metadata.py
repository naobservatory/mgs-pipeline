#!/usr/bin/env python3

accession_to_sample_id = {}
with open("raw-metadata.tsv") as inf:
    for line in inf:
        if not line.strip(): continue
        line = line[:-1]  # drop final newline
        accession, sample_id = line.split("\t")
        accession_to_sample_id[accession] = sample_id

# sample -> info
sample_info = {}
with open("sample_metadata.tsv") as inf:
    for line in inf:
        if not line.strip(): continue
        if line.startswith("sample"): continue  # header row
        line = line[:-1]  # drop final newline

        sample, core, location, method, date, mean_ct, ct_from_replicate, \
            sequencing = line.split("\t")

        if sequencing != "unenriched":
            continue

        sample = sample.strip()
        location = location.strip()
        method = method.strip()
        date = date.strip()

        mm, dd, yy = date.split("/")
        date = "20%s-%s-%s" % (yy, mm, dd)

        sample_info[sample] = {
            "location": location,
            "method": method,
            "date": date,
        }

with open("metadata.tsv", "w") as outf:
    for accession, sample_id in sorted(accession_to_sample_id.items()):
        outf.write("%s\t%s\t%s\t%s\n" % (
            accession,
            sample_info[sample_id]["location"],
            sample_info[sample_id]["date"],
            sample_info[sample_id]["method"]))
