#!/usr/bin/env python3

with open("raw_metadata.tsv") as inf:
    with open("metadata.tsv", "w") as outf:
        for line in inf:
            (
                run_accession,
                sample_accession,
                read_counts,
                sample_alias,
            ) = line.strip().split("\t")
            if run_accession == "run_accession":
                continue

            location_id, date = sample_alias.rsplit("_", 1)

            # Their docs don't say which is hospital A and which is B, but
            # looking at Figures S1 and S2 the number of samples and dates
            # match.
            location = {
                "EJH_MT": "A",
                "ALA_IN": "B",
            }[location_id]

            mm = date[0:2]
            dd = date[2:4]
            yy = date[4:6]

            outf.write(
                "%s\t20%s-%s-%s\t%s\n" % (run_accession, yy, mm, dd, location)
            )
