#!/usr/bin/env python3

accession_to_sample_id = {}
with open("raw-metadata.tsv") as inf:
    for line in inf:
        if not line.strip(): continue
        if line.startswith("run_accession"): continue
        line = line[:-1]  # drop final newline
        accession, sample_id, upload_date = line.split("\t")
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

        sample = sample.strip()
        if sample == "9_09_S1":
            # Based on the date column having 06/09/20 and this sample not
            # existing in SRA, I think this is very likely a typo for 6_09_S1
            sample = "6_09_S1"
        location = location.strip()
        method = method.strip()
        date = date.strip()

        mm, dd, yy = date.split("/")
        date = "20%s-%s-%s" % (yy, mm, dd)

        sample_info[sample] = {
            "location": location,
            "method": method,
            "sequencing": sequencing,
            "date": date,
        }

with open("metadata.tsv", "w") as outf:
    for accession, sample_id in sorted(
            accession_to_sample_id.items()):
        outf.write("%s\t%s\t%s\t%s\t%s\n" % (
            accession,
            sample_info[sample_id]["location"],
            sample_info[sample_id]["date"],
            sample_info[sample_id]["method"],
            sample_info[sample_id]["sequencing"]))
