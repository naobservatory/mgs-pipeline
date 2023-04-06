#!/usr/bin/env python3

sample_accession_to_collection_date = {}

with open("dataset-s4.tsv") as inf:
    for line in inf:
        if line.startswith("SI no."): continue # header

        bits = line.strip().split("\t")
        sample_accession = bits[2]
        collection_date = bits[20]

        sample_accession_to_collection_date[sample_accession] = collection_date

with open("raw_metadata.tsv") as inf:
    with open("metadata.tsv", "w") as outf:
        for line in inf:
            run_accession, sample_accession = line.strip().split("\t")
            collection_date = sample_accession_to_collection_date[
                sample_accession]
            outf.write("%s\t%s\n" % (
                run_accession, collection_date))


    
