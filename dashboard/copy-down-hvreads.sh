#!/usr/bin/env bash

# always run in dashboard dir
cd -- "$( dirname -- "${BASH_SOURCE[0]}")"

mkdir -p hvreads/
for bioproject in $(ls ../bioprojects/); do
    for hvread in $(aws s3 ls s3://nao-mgs/$bioproject/hvreads/ | \
                        awk '{print $NF}'); do
        if [ ! -e hvreads/$hvread ]; then
            echo s3://nao-mgs/$bioproject/hvreads/$hvread
        fi
    done
done | xargs -I {} -P 16 aws s3 cp {} hvreads/
    
                       
