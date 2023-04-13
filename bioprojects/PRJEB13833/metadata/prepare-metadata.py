#!/usr/bin/env python3
with open("raw_metadata.tsv") as inf:
    with open("metadata.tsv", "w") as outf:
        for line in inf:
            line = line.strip()
            if line.startswith("run"): continue

            run_accession, sample_accession, sample_alias, sample_title = \
                line.split("\t")

            _, yyyy, mm, dd, site = sample_title.split("_")

            outf.write("%s\t%s-%s-%s\tCluster %s\n" % (
                run_accession, yyyy, mm, dd, site))
