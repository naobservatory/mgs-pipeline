#!/usr/bin/env bash

bioproject="$1"
fname="$2"

accession=$(echo "$fname" | \
                awk -F/ '{print $NF}' | \
                sed s/_1[.]fastq[.]gz// | \
                sed s/[.]fastq[.]gz//)

variant=""
if [[ "$fname" == *"/nonhuman/"* ]]; then
    variant=".nonhuman"
fi

n_reads="bioprojects/${bioproject}metadata/$accession$variant.n_reads"

if [ -e $n_reads ]; then
    exit 0
fi

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
N_READS=$(aws s3 cp "$fname" - | gunzip | $SCRIPT_DIR/collect-n-reads-single.py)
if [ -z "${N_READS}" ]; then
   echo "Failed to count n_reads for $fname in $bioproject"
   exit 1
fi

echo "$N_READS" > "$n_reads"
