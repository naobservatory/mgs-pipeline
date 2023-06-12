#!/usr/bin/env bash

set -e  # exit on error

ROOT_DIR="$PWD"
if [ "$(basename $PWD)" == "mgs-pipeline" ]; then
    S3_DIR="s3://nao-mgs/"
    MGS_PIPELINE_DIR="$ROOT_DIR"
elif [ "$(basename $PWD)" == "mgs-restricted" ]; then
    S3_DIR="s3://nao-restricted/"
    MGS_PIPELINE_DIR=$ROOT_DIR/../mgs-pipeline
else
    echo "Run this from either mgs-pipeline or mgs-restricted"
    exit 1
fi

DASHBOARD_CODE_DIR=$MGS_PIPELINE_DIR/dashboard

cd dashboard

mkdir -p top_species_scratch/
for bioproject in $(ls ../bioprojects/); do
    CLADECOUNTS_DIR="$S3_DIR$bioproject/cladecounts"
    for cladecounts in $(aws s3 ls "$CLADECOUNTS_DIR/" | awk '{print $NF}'); do
        if [ ! -s top_species_scratch/${cladecounts/.gz} ]; then
            echo "$CLADECOUNTS_DIR/$cladecounts"
        fi
    done
done | xargs -P 32 -I {} $DASHBOARD_CODE_DIR/determine_comparison_species_single.sh {}

cat top_species_scratch/*.tsv | \
    awk '{print $2}' | \
    sort -n | \
    uniq \
        > comparison_taxids.txt

$DASHBOARD_CODE_DIR/determine_key_clades.py > key_clade_taxids.txt
                                                                  
mkdir -p top_species_counts/
for bioproject in $(ls ../bioprojects/); do
    CLADECOUNTS_DIR="$S3_DIR$bioproject/cladecounts"
    for cladecounts in $(aws s3 ls "$CLADECOUNTS_DIR/" | awk '{print $NF}'); do
        if [ ! -s top_species_counts/${cladecounts/.tsv.gz/.json} ]; then
            echo "$CLADECOUNTS_DIR/$cladecounts"
        fi
    done
done | xargs -P 32 -I {} $DASHBOARD_CODE_DIR/count_top_species_and_key_clades.sh {}
