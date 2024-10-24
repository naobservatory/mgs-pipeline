#!/usr/bin/env bash

curl -sS 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB49260&result=read_run&fields=read_count,sample_alias&format=tsv&download=true&limit=0' > raw_metadata.tsv
