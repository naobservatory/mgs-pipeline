#!/usr/bin/env python3
with open("metadata.tsv", "w") as outf:
    with open("raw_metadata.tsv") as inf:
        for line in inf:
            (run_accession,
             sample_accession,
             instrument_model,
             fastq_ftp,
             sample_alias,
             first_created) = line.strip().split("\t")

            if run_accession == "run_accession": continue
            
            if instrument_model == "MinION":
                continue  # Pipeline can't handle long-read yet

            location, raw_date = sample_alias.split("_")
            
            if len(raw_date) == 5:
                mm = "0" + raw_date[0]
                dd = raw_date[1:3]
                yy = raw_date[3:5]
            else:
                mm = raw_date[0:2]
                dd = raw_date[2:4]
                yy = raw_date[4:6]
                
            date = "20%s-%s-%s" % (yy, mm, dd)
            
            outf.write("%s\t%s\t%s\n" % (run_accession, location, date))
