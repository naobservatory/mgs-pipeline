#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    echo "Usage: ./import-accession.sh <accessions>"
fi

for ACCESSION in "$@"; do
    ENDPOINT="https://www.ebi.ac.uk/ena/portal/api/filereport"
    PARAMS="?accession=${ACCESSION}&fields=fastq_ftp&format=tsv"
    PARAMS="${PARAMS}&result=read_run&download=true&limit=0"

    curl -sS "${ENDPOINT}${PARAMS}" | \
        grep -v fastq_ftp | \
        awk '{print $NF}' | \
        tr ';' '\n' | \
        xargs -P 4 -I {} ./import-accession-single.sh $ACCESSION {}
done
