#!/usr/bin/env bash

if [ ! -e names.dmp ] ; then
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
    unzip taxdmp.zip
fi

for study in $(aws s3 ls s3://nao-mgs/ | awk '{print $NF}'); do
    for hv in $(aws s3 ls s3://nao-mgs/${study}humanviruses/ | \
                    awk '{print $NF}'); do
        if [ ! -e $hv ]; then
            echo s3://nao-mgs/${study}humanviruses/$hv
        fi
    done
done | xargs -I {} -P 32 aws s3 cp {} .

./prepare-dashboard-data.py
