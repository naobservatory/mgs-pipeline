#!/usr/bin/env bash

curl -sS 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA812772&result=read_run&fields=run_accession,sample_alias&format=tsv&download=true&limit=0' | grep -v run_accession > raw_metadata.tsv

./prepare_metadata.py
