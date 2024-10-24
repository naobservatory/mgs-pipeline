curl -sS 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA438174&result=read_run&fields=experiment_accession&format=tsv&download=true&limit=0' | grep -v "run" > raw_metadata.tsv

