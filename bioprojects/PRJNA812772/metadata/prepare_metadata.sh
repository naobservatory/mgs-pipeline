#!/usr/bin/env bash

curl -sS 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA812772&result=read_run&fields=run_accession,read_count,sample_alias,first_created&format=tsv&download=true&limit=0' | grep RNA | awk -F'\t' '{print $2"\t"$1}' > raw_metadata.tsv

./prepare_metadata.py
