#!/usr/bin/env python3

# library id -> plant, type
sample_info = {}

with open("table-s4.tsv") as inf:
    for line in inf:
        bits = line.strip().split("\t")
        if not bits[0].isdigit():
            continue

        sample_info[bits[0]] = bits[1], bits[2]

with open("raw_metadata.tsv") as inf:
    with open("metadata.tsv", "w") as outf:
        for line in inf:
            if line.startswith("run_accession"):
                continue

            run_accession, library_name = line.strip().split("\t")

            outf.write(
                "%s\t%s\t%s\n"
                % (
                    run_accession,
                    sample_info[library_name][0],
                    sample_info[library_name][1],
                )
            )
