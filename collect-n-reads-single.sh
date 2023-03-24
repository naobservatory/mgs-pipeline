#!/usr/bin/env bash

study="$1"
raw="$2"

accession=$(echo "$raw" | awk -F/ '{print $NF}' | awk -F_ '{print $1}')
n_reads="studies/${study}metadata/$accession.n_reads"

if [ -e $n_reads ]; then
    exit 0
fi

aws s3 cp "$raw" - | gunzip | grep -c ^@ > "$n_reads"
