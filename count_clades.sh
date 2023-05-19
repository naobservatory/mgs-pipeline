#!/usr/bin/env bash

set -euo pipefail

bioproject="$1"
sample="$2"

echo "Counting by clade for $bioproject $sample"
for kraken_file in $(
   aws s3 ls "s3://nao-mgs/$bioproject/processed/$sample" |
   awk '{print $NF}' | grep -v discarded); do
      aws s3 cp "s3://nao-mgs/$bioproject/processed/$kraken_file" - | gunzip
done | ./count_clades.py | \
    gzip | \
    aws s3 cp - "s3://nao-mgs/$bioproject/cladecounts/$sample.tsv.gz"
