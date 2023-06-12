#!/usr/bin/env bash
echo "count_top_species_and_key_clades.sh $1"
run_accession=$(basename $1 | sed s/.tsv.gz//)

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

aws s3 cp "$1" - \
    | gunzip \
    | $SCRIPT_DIR/count_top_species_and_key_clades.py \
        > top_species_counts/$run_accession.json
