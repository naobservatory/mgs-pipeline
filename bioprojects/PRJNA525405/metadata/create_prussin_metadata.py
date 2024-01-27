import re
import csv
from collections import defaultdict
sra_run_path = "sraruninfo.csv"
disease_tsv_path = "biosample.txt" 

metadata_file_path = "metadata.tsv"

with open(disease_tsv_path, mode='r') as infile:
    lines = infile.readlines()
    biosample_metadata = defaultdict(list)
    current_bio_sample = None
    next_line_is_description = None
    for line in lines:
        if 'BioSample:' in line:
            current_biosample = line.split('BioSample:')[1].split(';')[0].strip()
        if 'season' in line:
            
            season = re.search(r'/season="(.*?)"',line).group(1).strip()
        if 'type_of_nucleic_acid' in line:
            na_type = re.search(r'/type_of_nucleic_acid="AP-(DNA|RNA)', line).group(1).strip()
        if next_line_is_description:
            description = line.strip()
            next_line_is_description = False
        if 'Description' in line:
            next_line_is_description = True
        if 'Accession' in line:
            biosample_metadata[current_biosample].extend([
                season,
                na_type,
                description])


with open(sra_run_path, mode='r', newline='') as infile, open(metadata_file_path, mode='w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile, delimiter='\t')

    next(reader)

    for line in reader:
        run = line[0]
        biosample = line[6]
        date = line[10]
        season, na_type, description = biosample_metadata[biosample]

        writer.writerow([run, na_type, date, description, season])
