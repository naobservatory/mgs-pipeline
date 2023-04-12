#!/usr/bin/env bash

set -euo pipefail

bioproject="$1"
sample="$2"

echo "Counting taxids for $bioproject $sample"

for kraken_file in $(
   aws s3 ls "s3://nao-mgs/$bioproject/processed/$sample" |
   awk '{print $NF}'); do
      aws s3 cp "s3://nao-mgs/$bioproject/processed/$kraken_file" - | gunzip
done | ./count_taxonomic_ids.py | \
    gzip | \
    aws s3 cp - "s3://nao-mgs/$bioproject/allcounts/$sample.tsv.gz"
