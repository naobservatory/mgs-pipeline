#!/usr/bin/env python3

rows = []
with open("run_aliases.txt") as inf:
    for line in inf:
        run_accession, run_alias = line.strip().split("\t")
        if run_accession == "run_accession": continue

        location = run_alias.split("_")[-3]

        # Cities from Supplementary Table 7.
        country, city = {
            "China": ("China", "Beijing"),
            # city from looking up where GL 782 flies from.
            "Greenland": ("Greenland", "Kangerlussuaq"),
            "Japan": ("Japan", "Tokyo"),
            "Newark": ("USA", "Newark, NJ"),  # Could call this NYC
            "Pakistan": ("Pakistan", "Islamabad"),
            "Singapore": ("Singapore", "Singapore"),
            "Thailand": ("Thailand", "Bankok"),
            "Toronto": ("Canada", "Toronto"),
            "Washington": ("USA", "Washington DC"),
        }[location]
        
        # I don't see how to get sample dates from the paper.  Supplementary
        # table 7 seems like it should have this information, but I don't see
        # any mapping from sample IDs ("Thailand 2") to rows in the table.

        rows.append((run_accession, country, city))

rows.sort()

with open("../../bioprojects/PRJEB12466/metadata/metadata.tsv", "w") as outf:
    for row in rows:
        outf.write("\t".join(row))
        outf.write("\n")
