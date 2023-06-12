#!/usr/bin/env bash
echo "determine_comparison_species_single.sh $1"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

run_accession=$(basename $1 | sed s/.tsv.gz//)
aws s3 cp "$1" - \
    | gunzip \
    | $SCRIPT_DIR/determine_comparison_species.py \
        > top_species_scratch/$run_accession.tsv
