#!/usr/bin/env bash

set -euo pipefail

s3_bucket="$1"
delivery="$2"
sample="$3"
reference="$4"

./download-taxonomy.sh

echo "Counting by clade for $delivery $sample"
for kraken_file in $(
   aws s3 ls "$s3_bucket/$delivery/processed-$reference/$sample" |
   awk '{print $NF}' | grep -v discarded); do
    aws s3 cp "$s3_bucket/$delivery/processed-$reference/$kraken_file" - \
        | gunzip
done | ./count_clades.py | \
    gzip | \
    aws s3 cp - "$s3_bucket/$delivery/cladecounts-$reference/$sample.tsv.gz"
