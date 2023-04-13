#!/usr/bin/env python3

sample_metadata = {}

with open("table-s1.tsv") as inf:
    for line in inf:
        line = line.strip()
        if line.startswith("Sample type"): continue

        bits = line.split("\t")
        sample_type, date_sampled, _, _, _, _, sra_experiment_accession = bits

        month, year = date_sampled.split()

        month = {
            "August": "08",
            "October": "10",
            "November": "11",
            "January": "01",
            "March": "03",
            "May": "05",
        }[month]

        date_sampled = "%s-%s" % (year, month)
        
        sample_metadata[sra_experiment_accession] = sample_type, date_sampled

with open("raw_metadata.tsv") as inf:
    with open("metadata.tsv", "w") as outf:
        for line in inf:
            run_accession, sra_experiment_accession = line.strip().split("\t")
            outf.write("%s\t%s\t%s\n" % (
                run_accession, *sample_metadata[sra_experiment_accession]))
        
