#!/usr/bin/env/ python3
import csv
import re
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
            run = line['Run'].strip() 
            country = line['geo_loc_name_country'].strip()
            continent = line['geo_loc_name_country_continent'].strip()
            sampling_device = line['samp_collect_device'].strip()
            raw_location =  line['geo_loc_name'].strip()
            date = line["Collection_Date"].strip()
            location_bits = re.split(r":|\s*,\s*",raw_location)
            if len(location_bits) == 3:
                location, sub_city_location = location_bits[1].replace("\\", ''), location_bits[2]
            elif len(location_bits) == 4:
                location, sub_city_location = location_bits[1].replace("\\", ''), location_bits[3]
            else:
                location, sub_city_location = None, None
            if sub_city_location:
                fine_location = f"{location}, {sub_city_location}"
            writer.writerow([run, country, continent, location, fine_location, sampling_device, date])
