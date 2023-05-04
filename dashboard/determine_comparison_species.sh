#!/usr/bin/env bash

set -e  # exit on error

# always run in dashboard dir
cd -- "$( dirname -- "${BASH_SOURCE[0]}")"

mkdir -p top_species_scratch/
for bioproject in $(ls ../bioprojects/); do
    CLADECOUNTS_DIR="s3://nao-mgs/$bioproject/cladecounts"
    for cladecounts in $(aws s3 ls "$CLADECOUNTS_DIR/" | awk '{print $NF}'); do
        if [ ! -s top_species_scratch/${cladecounts/.gz} ]; then
            echo "$CLADECOUNTS_DIR/$cladecounts"
        fi
    done
done | xargs -P 32 -I {} ./determine_comparison_species_single.sh {}

cat top_species_scratch/*.tsv | \
    awk '{print $2}' | \
    sort -n | \
    uniq \
        > comparison_taxids.txt

./determine_key_clades.py > key_clade_taxids.txt
                                                                  
mkdir -p top_species_counts/
for bioproject in $(ls ../bioprojects/); do
    CLADECOUNTS_DIR="s3://nao-mgs/$bioproject/cladecounts"
    for cladecounts in $(aws s3 ls "$CLADECOUNTS_DIR/" | awk '{print $NF}'); do
        if [ ! -s top_species_counts/${cladecounts/.tsv.gz/.json} ]; then
            echo "$CLADECOUNTS_DIR/$cladecounts"
        fi
    done
done | xargs -P 32 -I {} ./count_top_species_and_key_clades.sh {}

echo "Now run ./prepare-dashboard-data.sh"
