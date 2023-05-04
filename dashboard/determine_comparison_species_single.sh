#!/usr/bin/env bash
echo "determine_comparison_species_single.sh $1"
run_accession=$(basename $1 | sed s/.tsv.gz//)
aws s3 cp "$1" - \
    | gunzip \
    | ./determine_comparison_species.py \
        > top_species_scratch/$run_accession.tsv
