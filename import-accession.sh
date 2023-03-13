#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    echo "Usage: ./import-accession.sh <accessions>"
fi

for ACCESSION in "$@"; do
    ENDPOINT="https://www.ebi.ac.uk/ena/portal/api/filereport"
    FIELDS="fastq_ftp"
    PARAMS="?accession=${ACCESSION}&fields=${FIELDS}&format=tsv"
    PARAMS="${PARAMS}&result=read_run&download=true&limit=0"

    curl -sS "${ENDPOINT}${PARAMS}" | \
        grep -v fastq_ftp | \
        awk '{print $NF}' | \
        tr ';' '\n' | \
        xargs -P 4 -I {} bash -c \
          'if ! aws s3 ls s3://nao-mgs/'$ACCESSION'/raw/$(basename {}) \
               &>/dev/null; then \
             curl -sS {} | \
               aws s3 cp - s3://nao-mgs/'$ACCESSION'/raw/$(basename {}); \
           fi'
done
