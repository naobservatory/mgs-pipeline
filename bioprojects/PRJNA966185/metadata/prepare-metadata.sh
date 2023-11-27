#!/usr/bin/env bash

if [ ! -e wastewater_sample_metadata_table1.tsv ]; then
    wget https://zenodo.org/record/7884454/files/data.zip
    unzip data.zip
    rm data.zip
    rm -r __MACOSX
    mv data/wastewater_sample_metadata_table1.tsv .
    mv data/no-probe_with-TWIST-probe_control_read_alignment_table1.tsv .
    rm -r data/    
fi

if  [ ! -e metadata-raw.tsv ]; then
  curl -sS \
  'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA966185&result=read_run&fields=run_accession,run_alias,sample_alias&format=tsv&download=true&limit=0' \
  | sort > metadata-raw.tsv
fi

./prepare-metadata.py
