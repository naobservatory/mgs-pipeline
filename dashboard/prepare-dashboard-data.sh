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

$MGS_PIPELINE_DIR/collect-n-reads.sh

cd $ROOT_DIR/dashboard
mkdir -p allmatches/
mkdir -p hvreads/
mkdir -p hvrfull/
mkdir -p ribocounts/

if [ ! -e names.dmp ] ; then
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2022-12-01.zip
    unzip taxdmp_2022-12-01.zip
fi

cd $ROOT_DIR
$MGS_PIPELINE_DIR/dashboard/determine_comparison_species.sh

cd $ROOT_DIR/dashboard
for study in $(aws s3 ls $S3_DIR | awk '{print $NF}'); do
    for am in $(aws s3 ls $S3_DIR${study}allmatches/ | \
                    awk '{print $NF}'); do
        if [ ! -s allmatches/$am ]; then
            echo $S3_DIR${study}allmatches/$am
        fi
    done
done | xargs -I {} -P 32 aws s3 cp {} allmatches/

for study in $(aws s3 ls $S3_DIR | awk '{print $NF}'); do
    for hvr in $(aws s3 ls $S3_DIR${study}hvreads/ | \
                    awk '{print $NF}'); do
        if [ ! -s hvreads/$hvr ]; then
            echo $S3_DIR${study}hvreads/$hvr
        fi
    done
done | xargs -I {} -P 32 aws s3 cp {} hvreads/

for study in $(aws s3 ls $S3_DIR | awk '{print $NF}'); do
    for rc in $(aws s3 ls $S3_DIR${study}ribocounts/ | \
                    awk '{print $NF}'); do
    	if [ ! -s ribocounts/$rc ]; then
	    echo $S3_DIR${study}ribocounts/$rc
	fi
     done
done | xargs -I {} -P 32 aws s3 cp {} ribocounts/

$MGS_PIPELINE_DIR/dashboard/prepare-dashboard-data.py $ROOT_DIR $MGS_PIPELINE_DIR

echo "Now check in data.js and the json files and check out on prod"
echo "Then run copy-down-hvreads.sh and send the result to prod out of band"
