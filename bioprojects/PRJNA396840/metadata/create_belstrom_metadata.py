import csv
from collections import defaultdict
input_file_path = "sraruninfo.tsv"
disease_tsv_path = "biosample.txt" 

metadata_file_path = "metadata.tsv"

with open(disease_tsv_path, mode='r') as infile:
    lines = infile.readlines()
    metadata_biosample = defaultdict(dict) 
    for line in lines:
        if 'BioSample:' in line:
            current_biosample = line.split('BioSample:')[1].split(';')[0].strip()
        if line.strip() in ["Caries", "Periodontitis", "Healthy control"]:
            disease_state = line.strip()
        if 'Accession:' in line:
            metadata_biosample[current_biosample] = disease_state 


with open(input_file_path, mode='r', newline='') as infile, open(metadata_file_path, mode='w', newline='') as outfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')

    next(reader)

    for line in reader:
        run = line[0]
        library_source = line[14]
        if library_source == "METATRANSCRIPTOMIC":
            na_type="RNA"
        else:
            na_type="DNA"
        bio_sample = line[24]
        disease_status = metadata_biosample[bio_sample] 

        writer.writerow([run, na_type, disease_status])
