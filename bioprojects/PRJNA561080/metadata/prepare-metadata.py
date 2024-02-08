#!/usr/bin/env/ python3
import csv

sra_run_path = "sra_run_table.csv"
# Download from here:
# https://www.ncbi.nlm.nih.gov/Traces/study/?page=6&query_key=2&WebEnv=MCID_65c52c5be69ed911ee3b3821&o=acc_s%3Aa#:~:text=Bases-,Download,-Cloud%20Data%20Delivery
# Click on Download > "Metadata" in the row "Total"




with open("metadata.tsv", "w") as outf:
    with open(sra_run_path , mode='r') as inf:
        reader = csv.DictReader(inf)
        writer = csv.writer(outf, delimiter='\t')
        next(reader) # skipping header

        for line in reader:
            run = line['Run'] 
            country = line['geo_loc_name_country']
            continent = line['geo_loc_name_country_continent']
            sampling_device = line['samp_collect_device']
            raw_location =  row['geo_loc_name']
            
            fine_location = 
            
