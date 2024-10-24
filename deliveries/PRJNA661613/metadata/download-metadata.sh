#!/usr/bin/env bash
curl -sS 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA661613&result=read_run&fields=run_accession,read_count,sample_alias,first_created&format=tsv&download=true&limit=0' | awk -F'\t' '{print $2"\t"$4"\t"$NF}' > raw-metadata.tsv
wget 'https://raw.githubusercontent.com/alexcritschristoph/wastewater_sarscov2/master/data/sample_metadata.tsv'
./process-metadata.py
