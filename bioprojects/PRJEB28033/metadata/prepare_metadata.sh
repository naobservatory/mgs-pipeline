#!/usr/bin/env bash

ENDPOINT="https://www.ebi.ac.uk/ena/portal/api/filereport"
ACCESSION="PRJEB28033"
FIELDS="fastq_ftp,run_alias"
PARAMS="?accession=${ACCESSION}&fields=${FIELDS}&format=tsv"
PARAMS="${PARAMS}&result=read_run&download=true&limit=0"

RAW_METADATA=raw_metadata.tsv
COMMON_METADATA=metadata.tsv

if [[ ! -e $RAW_METADATA ]]; then
    wget "${ENDPOINT}${PARAMS}" -O raw_metadata.tsv
fi

if [[ ! -e "$COMMON_METADATA" ]]; then
    ./parse_metadata.py $RAW_METADATA $COMMON_METADATA
fi

echo "Parsed output in $COMMON_METADATA"
