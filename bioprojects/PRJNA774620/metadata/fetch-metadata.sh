curl -sS 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA774620&result=read_run&fields=library_strategy,read_count,run_alias,fastq_ftp,sample_alias&format=tsv&download=true&limit=0' > raw_metadata.tsv
