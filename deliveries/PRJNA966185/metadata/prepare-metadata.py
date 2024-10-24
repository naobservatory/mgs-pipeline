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
            run_accession, run_alias, _, sample_id = line.split("\t")

            record = [
                run_accession,
                sample_id,
                "0" if "unenriched" in run_alias else "1",
            ]
            record.extend(sample_info[sample_id])
            outf.write("\t".join(record) + "\n")
