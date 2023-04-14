#!/usr/bin/env bash

set -e  # exit on error

# always run in dashboard dir
cd -- "$( dirname -- "${BASH_SOURCE[0]}")"

if [ ! -e names.dmp ] ; then
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2022-12-01.zip
    unzip taxdmp_2022-12-01.zip
fi

for study in $(aws s3 ls s3://nao-mgs/ | awk '{print $NF}'); do
    for hv in $(aws s3 ls s3://nao-mgs/${study}humanviruses/ | \
                    awk '{print $NF}'); do
        if [ ! -e humanviruses/$hv ]; then
            echo s3://nao-mgs/${study}humanviruses/$hv
        fi
    done
done | xargs -I {} -P 32 aws s3 cp {} humanviruses/

./prepare-dashboard-data.py

echo "Now check in data.js and the json files and check out on prod"
echo "Then run copy-down-hvreads.sh and send the result to prod out of band"
echo
echo "(Paper not showing up?  Did you remember to run ./collect-n-reads.sh ?)"
