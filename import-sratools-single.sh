#!/usr/bin/env bash

BIOPROJECT=$1
SAMPLE=$2

if [[ "$SAMPLE" != SRR* ]]; then
    echo "Bad sample '$SAMPLE' for $BIOPROJECT"
    exit 1
fi

~/sratoolkit.3.0.5-centos_linux64/bin/prefetch $SAMPLE
~/sratoolkit.3.0.5-centos_linux64/bin/fasterq-dump $SAMPLE

for fastq in $SAMPLE*.fastq; do
    gzip $fastq
    aws s3 cp $fastq.gz s3://nao-mgs/$BIOPROJECT/raw/
    rm $fastq.gz
done

rm -r $SAMPLE
