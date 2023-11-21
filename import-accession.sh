#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    echo "Usage: ./import-accession.sh <accessions>"
fi

metadata_file="metadata-ncbi.tsv"

for ACCESSION in "$@"; do
    # Define the directory based on BioProject ID
    metadata_dir="bioprojects/${ACCESSION}/metadata"
    mkdir -p "$metadata_dir"
    metadata_file="${metadata_dir}/metadata-ncbi.tsv"

    ENDPOINT="https://www.ebi.ac.uk/ena/portal/api/filereport"
    PARAMS="?accession=${ACCESSION}&fields=fastq_ftp&format=tsv"
    PARAMS="${PARAMS}&result=read_run&download=true&limit=0"

    curl -sS "${ENDPOINT}${PARAMS}" | \
        grep -v fastq_ftp | \
        awk '{print $NF}' | \
        tr ';' '\n' | \
        xargs -P 4 -I {} ./import-accession-single.sh $ACCESSION {}

    # Fetch and append metadata
    echo "Fetching metadata for $ACCESSION"
    (esearch -db sra -query "$ACCESSION" | efetch -format runinfo || echo "$ACCESSION") >> "$metadata_file"
done
