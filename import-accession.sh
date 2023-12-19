#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    echo "Usage: ./import-accession.sh <accessions>"
    exit 1
fi

# Function to check if E-utilities are installed
check_eutilities() {
    if ! command -v esearch &> /dev/null; then
        echo "E-utilities are not installed. Please refer to the dependencies section of the README."
        return 1
    else
        return 0
    fi
}

for ACCESSION in "$@"; do
    ENDPOINT="https://www.ebi.ac.uk/ena/portal/api/filereport"
    PARAMS="?accession=${ACCESSION}&fields=fastq_ftp&format=tsv"
    PARAMS="${PARAMS}&result=read_run&download=true&limit=0"

    # Fetch SRA data and check if any files are returned
    SRA_FILES=$(curl -sS "${ENDPOINT}${PARAMS}" | grep -v fastq_ftp)

    if [ -n "$SRA_FILES" ]; then
        # Define the directory based on BioProject ID
        metadata_dir="bioprojects/${ACCESSION}/metadata"
        mkdir -p "$metadata_dir"

        # Download SRA data
        echo "$SRA_FILES" | awk '{print $NF}' | tr ';' '\n' | xargs -P 4 -I {} ./import-accession-single.sh $ACCESSION {}

        # Create metadata.tsv with SRA run names
        echo "$SRA_FILES" | awk -F'[/.]' '{print $(NF-2)}' > "${metadata_dir}/metadata.tsv"

        if check_eutilities; then
            # Fetch metadata from NCBI
            echo "Fetching metadata for $ACCESSION from NCBI"
            METADATA=$(esearch -db sra -query "$ACCESSION" | efetch -format runinfo)

            if [ -n "$METADATA" ] && [ "$(echo "$METADATA" | grep -v '^Run,' | wc -l)" -gt 0 ]; then
                echo "$METADATA" > "${metadata_dir}/metadata-ncbi.tsv"
            else
                echo "No metadata was found for $ACCESSION on NCBI."
            fi
        fi
    else
        echo "No files found for accession $ACCESSION"
    fi
done

