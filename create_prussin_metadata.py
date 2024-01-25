import csv
from collections import defaultdict
input_file_path = "sraruninfo.tsv"
disease_tsv_path = "biosample.txt" 

metadata_file_path = "bioprojects/PRJNA525405/metadata/metadata.tsv"
sras = defaultdict(list)


with open(metadata_file_path, mode='r', newline='') as infile:
    reader = csv.reader(infile, delimiter='\t') 
    for line in reader:
        sra_id = line[0]
        sras[sra_id] = []



with open(disease_tsv_path, mode='r') as infile:
    lines = infile.readlines()

    current_bio_sample = None
    for line in lines:
        if 'BioSample:' in line:
            current_biosample = line.split('BioSample:')[1].split(';')[0].strip()
        if line in ["Caries", "Periodontitis", "Healthy Control"]:
            sras[current_biosample].append(line.strip())
with open(input_file_path, mode='r', newline='') as infile, open(metadata_file_path, mode='a', newline='') as outfile:
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

        bio_sample = line [25]
        if run in sras:
            sras[run].append(na_type)
            sras[run].append(bio_sample)




writer.writerow([run, na_type, bio_sample])
