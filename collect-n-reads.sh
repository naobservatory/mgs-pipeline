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
    if [ ! -d bioprojects/${bioproject}metadata/ ] ; then
        continue
    fi

    for indir in $S3_DIR${bioproject}raw/ $S3_DIR${bioproject}nonhuman/; do
        n_1=$(aws s3 ls $indir | grep "_1.fastq" | wc -l)
        if [ $n_1 == 0 ]; then
            for fname in $(aws s3 ls $indir | awk '{print $NF}'); do
                echo $indir$fname
            done
        else
            for fname in $(aws s3 ls $indir | awk '{print $NF}' | \
                               grep "_1.fastq"); do
                echo $indir$fname
            done
        fi
    done | xargs -I {} -P 32 $SCRIPT_DIR/collect-n-reads-single.sh $bioproject {}
done 
