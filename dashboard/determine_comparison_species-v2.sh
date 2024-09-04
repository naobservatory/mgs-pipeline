#!/usr/bin/env bash

set -e  # exit on error

REFSUFFIX="$1"

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

function list_bioprojects_and_merged_kraken_reports() {
    for bioproject in $(ls ../bioprojects/); do
        PAPER=$bioproject
        NAME_FNAME=../bioprojects/$bioproject/metadata/name.txt
        if [[ -e $NAME_FNAME ]]; then
            PAPER=$(cat $NAME_FNAME)
        fi
        
        MGS_WORKFLOW_OUTPUT_FNAME=../papers/$PAPER/mgs-workflow-output.txt
        if [[ ! -e $MGS_WORKFLOW_OUTPUT_FNAME ]]; then
            continue
        fi
        
        echo $bioproject $(cat $MGS_WORKFLOW_OUTPUT_FNAME)/output/results/taxonomy/kraken_reports_merged.tsv.gz
    done
}

mkdir -p top_species_scratch_v2/
list_bioprojects_and_merged_kraken_reports | \
    xargs -P 32 -I {} $DASHBOARD_CODE_DIR/determine_comparison_species_single_v2.sh {}

cat top_species_scratch_v2/*.tsv | \
    awk '{print $2}' | \
    sort -n | \
    uniq \
        > comparison_taxids_v2.txt

$DASHBOARD_CODE_DIR/determine_key_clades.py > key_clade_taxids.txt
                                                                  
mkdir -p top_species_counts_v2/
list_bioprojects_and_merged_kraken_reports | \
    xargs -P 32 -I {} $DASHBOARD_CODE_DIR/count_top_species_and_key_clades_v2.sh {}
