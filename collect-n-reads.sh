#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [ "$(basename $PWD)" == "mgs-pipeline" ]; then
    S3_DIR="s3://nao-mgs/"
elif [ "$(basename $PWD)" == "mgs-restricted" ]; then
    S3_DIR="s3://nao-restricted/"
else
    echo "Run this from either mgs-pipeline or mgs-restricted"
    exit 1
fi

for bioproject in $(aws s3 ls $S3_DIR | awk '{print $NF}'); do
    for raw in $(aws s3 ls $S3_DIR${bioproject}raw/ | \
                     awk '{print $NF}' | \
                     grep "_1.fastq"); do
        echo $S3_DIR${bioproject}raw/$raw
    done | xargs -I {} -P 32 $SCRIPT_DIR/collect-n-reads-single.sh $bioproject {}
done 
