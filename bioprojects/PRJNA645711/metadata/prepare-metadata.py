#!/usr/bin/env python3

with open("raw_metadata.tsv") as inf:
    with open("metadata.tsv", "w") as outf:
        for line in inf:
            (
                run_accession,
                sample_accession,
                run_alias,
                sample_alias,
                sample_title,
            ) = line.strip().split("\t")

            if run_accession == "run_accession":
                continue

            outf.write("%s\t%s\n" % (run_accession, sample_title.split()[-1]))
