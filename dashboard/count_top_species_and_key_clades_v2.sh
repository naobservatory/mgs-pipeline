#!/usr/bin/env bash
echo "count_top_species_and_key_clades_v2.sh $1"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

DELIVERY=$(echo "$1" | awk '{print $1}')
KRAKEN_REPORTS_MERGED=$(echo "$1" | awk '{print $2}')

if [[ -z "$DELIVERY" ]] || [[ ! -e "../bioprojects/$DELIVERY/" ]]; then
    echo "Couldn't calculate delivery from $1; got $DELIVERY"
    exit 1
fi

aws s3 cp "$KRAKEN_REPORTS_MERGED" - \
    | gunzip \
    | $SCRIPT_DIR/count_top_species_and_key_clades_v2.py \
        > top_species_counts_v2/$DELIVERY.json
