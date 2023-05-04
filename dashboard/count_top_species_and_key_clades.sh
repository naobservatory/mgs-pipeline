#!/usr/bin/env bash
echo "count_top_species_and_key_clades.sh $1"
run_accession=$(basename $1 | sed s/.tsv.gz//)
aws s3 cp "$1" - \
    | gunzip \
    | ./count_top_species_and_key_clades.py \
        > top_species_counts/$run_accession.json
