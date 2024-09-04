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

METADATA_ONLY=false
if [ $# -gt 0 ]; then
    if [[ "$1" = "--metadata-only" ]]; then
        METADATA_ONLY=true
    else
        echo "Argument $1 not understood"
    fi
fi

$MGS_PIPELINE_DIR/collect-n-reads.sh

REFSUFFIX=$(cat $MGS_PIPELINE_DIR/reference-suffix.txt)

cd $ROOT_DIR/dashboard

$MGS_PIPELINE_DIR/dashboard/prepare-dashboard-metadata.py $ROOT_DIR $MGS_PIPELINE_DIR

if $METADATA_ONLY; then
    echo "Metadata complete; exiting early as requested."
    exit 0
fi

cd $ROOT_DIR
$MGS_PIPELINE_DIR/dashboard/determine_comparison_species-v2.sh $REFSUFFIX

cd $ROOT_DIR/dashboard

mkdir -p hv_hits_putative_collapsed
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

    aws s3 cp \
        $(cat $MGS_WORKFLOW_OUTPUT_FNAME)/output/results/hv/hv_hits_putative_collapsed.tsv.gz \
        hv_hits_putative_collapsed/$bioproject.tsv.gz
done

$MGS_PIPELINE_DIR/dashboard/prepare-dashboard-data-v2.py \
    $ROOT_DIR $MGS_PIPELINE_DIR

echo
echo "Now check in your changes and send for review."
echo
echo "After review is complete and the changes are merged, you'll need to"
echo "ssh into data.securebio.org and run:"
echo "    cd $(basename $ROOT_DIR)/"
echo "    git pull"
echo "To update the dashboard."
