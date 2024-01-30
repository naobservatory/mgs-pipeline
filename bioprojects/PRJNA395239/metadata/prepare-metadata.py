#!/usr/bin/env python3
import csv

sra_run_path = "sra_run_table.csv"
# Download from here: 
# https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=2&WebEnv=MCID_65b917aa62f9732e9951754e&f=run_file_create_date_dt%3An&o=acc_s%3Aa#:~:text=Bases-,Download,-Cloud%20Data%20Delivery
# Click on Download > "Metadata"


with open("metadata.tsv", "w") as outf:
    with open(sra_run_path, mode='r') as inf:
        reader = csv.reader(inf)
        writer = csv.writer(outf, delimiter='\t')
        next(reader)

        for line in reader:
            run = line[0]
            sample_origin, na_type = line[-4].strip().split("_") 
            print(run, sample_origin, na_type)
            writer.writerow([run, sample_origin, na_type])
