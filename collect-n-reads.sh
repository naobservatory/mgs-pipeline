#!/usr/bin/env bash

for bioproject in $(aws s3 ls s3://nao-mgs/ | awk '{print $NF}'); do
    for raw in $(aws s3 ls s3://nao-mgs/${bioproject}raw/ | \
                     awk '{print $NF}' | \
                     grep "_1.fastq"); do
        echo s3://nao-mgs/${bioproject}raw/$raw
    done | xargs -I {} -P 32 ./collect-n-reads-single.sh $bioproject {}
done 
