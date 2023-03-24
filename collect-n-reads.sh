#!/usr/bin/env bash

for study in $(aws s3 ls s3://nao-mgs/ | awk '{print $NF}'); do
    for raw in $(aws s3 ls s3://nao-mgs/${study}raw/ | \
                     awk '{print $NF}' | \
                     grep "_1.fastq"); do
        echo s3://nao-mgs/${study}raw/$raw
    done | xargs -I {} -P 32 ./collect-n-reads-single.sh $study {}
done 
