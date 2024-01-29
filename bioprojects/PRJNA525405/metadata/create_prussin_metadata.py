import re
import csv
from collections import defaultdict
sra_run_path = "sraruninfo.csv"
# Download from here: https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=4&WebEnv=MCID_65b7f332bd5d126b3ce13457&o=acc_s%3Aa#:~:text=Bases-,Download,-Cloud%20Data%20Delivery 
# Click on Download > "Metadata"
disease_tsv_path = "biosample.txt" 
# Download from here: https://www.ncbi.nlm.nih.gov/biosample?Db=biosample&DbFrom=bioproject&Cmd=Link&LinkName=bioproject_biosample&LinkReadableName=BioSample&ordinalpos=1&IdsFromResult=525405
# First, show all samples on one page. Then, next to the "{n} per page" button, 
# click on "Summary" and check "Full(txt)". Copy paste that data into new file 
# with name "biosample.txt".
metadata_file_path = "metadata.tsv"

with open(disease_tsv_path, mode='r') as infile:
    lines = infile.readlines()
    biosample_metadata = defaultdict(list)
    current_bio_sample = None
    next_line_is_sampling_range = None
    for line in lines:
        if 'BioSample:' in line:
            current_biosample = line.split('BioSample:')[1].split(';')[0].strip()
        if 'season' in line:
            season = re.search(r'/season="(.*?)"',line).group(1).strip()
        if 'type_of_nucleic_acid' in line:
            na_type = re.search(r'/type_of_nucleic_acid="AP-(DNA|RNA)', line).group(1).strip()
        if next_line_is_sampling_range:
            sampling_range = line.strip()
            next_line_is_sampling_range = False
        if 'Description' in line:
            next_line_is_sampling_range = True
        if 'Accession' in line:
            biosample_metadata[current_biosample].extend([
                season,
                na_type,
                sampling_range])


with open(sra_run_path, mode='r', newline='') as infile, open(metadata_file_path, mode='w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile, delimiter='\t')

    next(reader)

    for line in reader:
        run = line[0]
        biosample = line[6]
        date = line[10]
        season, na_type, sampling_range = biosample_metadata[biosample]

        writer.writerow([run, na_type, date, sampling_range, season])
