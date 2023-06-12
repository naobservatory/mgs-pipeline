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

cd $ROOT_DIR/dashboard
mkdir -p humanviruses/
mkdir -p cladecounts/

if [ ! -e names.dmp ] ; then
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2022-12-01.zip
    unzip taxdmp_2022-12-01.zip
fi

for study in $(aws s3 ls $S3_DIR | awk '{print $NF}'); do
    for hv in $(aws s3 ls $S3_DIR${study}humanviruses/ | \
                    awk '{print $NF}'); do
        if [ ! -s humanviruses/$hv ]; then
            echo $S3_DIR${study}humanviruses/$hv
        fi
    done
done | xargs -I {} -P 32 aws s3 cp {} humanviruses/

$MGS_PIPELINE_DIR/dashboard/prepare-dashboard-data.py $ROOT_DIR $MGS_PIPELINE_DIR

echo "Now check in data.js and the json files and check out on prod"
echo "Then run copy-down-hvreads.sh and send the result to prod out of band"
echo
echo "(Paper not showing up?  Did you remember to run ./collect-n-reads.sh ?)"
