#!/usr/bin/env bash

bioproject="$1"
raw="$2"

accession=$(echo "$raw" | \
                awk -F/ '{print $NF}' | \
                awk -F_ '{print $1}' | \
                sed s/.fastq.gz//)
n_reads="bioprojects/${bioproject}metadata/$accession.n_reads"

if [ -e $n_reads ]; then
    exit 0
fi

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
aws s3 cp "$raw" - | gunzip | $SCRIPT_DIR/collect-n-reads-single.py > "$n_reads"
