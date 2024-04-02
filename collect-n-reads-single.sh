#!/usr/bin/env bash

bioproject="$1"
raw="$2"

accession=$(echo "$raw" | \
                awk -F/ '{print $NF}' | \
                sed s/_1[.]fastq[.]gz// | \
                sed s/[.]fastq[.]gz//)
n_reads="bioprojects/${bioproject}metadata/$accession.n_reads"

if [ -e $n_reads ]; then
    exit 0
fi

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
N_READS=$(aws s3 cp "$raw" - | gunzip | $SCRIPT_DIR/collect-n-reads-single.py)
if [ -z "${N_READS}" ]; then
   echo "Failed to count n_reads for $raw in $bioproject"
   exit 1
fi

echo "$N_READS" > "$n_reads"
