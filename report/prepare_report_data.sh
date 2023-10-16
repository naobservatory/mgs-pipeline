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

cd $ROOT_DIR/report

mkdir -p cladecounts

for study in $(aws s3 ls $S3_DIR | awk '{print $NF}'); do
    for cc in $(aws s3 ls $S3_DIR${study}cladecounts/ | \
                    awk '{print $NF}'); do
        if [ ! -s cladecounts/$cc ]; then
            echo $S3_DIR${study}cladecounts/$cc
        fi
    done
done | xargs -I {} -P 32 aws s3 cp {} cladecounts/

$MGS_PIPELINE_DIR/report/prepare-report-data.py $ROOT_DIR $MGS_PIPELINE_DIR

# TODO 
# 1. make prepare-report_data take a bioproject and sample arg. No sample arg is given, run report on all samples in biorpoject
# 2. make a temporary directory where tmp data like cladecounts can be stored
# 3. copy neccessary data over from aws, generate stat, and send it to a sample or bioproject specific metadata or report statistics directory
#
