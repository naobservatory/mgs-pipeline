#!/usr/bin/env python3

sample_info = {}
with open("wastewater_sample_metadata_table1.tsv") as inf:
    for line in inf:
        line = line.rstrip("\n")
        if line.startswith("sample_ID"):
            continue
        sample_id, site, city, date, flow_rate = line.split("\t")
        sample_info[sample_id] = [site, city, date, flow_rate]

with open("metadata-raw.tsv") as inf:
    with open("metadata.tsv", "w") as outf:
        for line in inf:
            line = line.rstrip("\n")
            if line.startswith("run_accession"):
                continue
            run_accession, _, sample_id = line.split("\t")

            record = [run_accession, sample_id]
            record.extend(sample_info[sample_id])
            outf.write("\t".join(record) + "\n")
