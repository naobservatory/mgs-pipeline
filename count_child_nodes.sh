#!/usr/bin/env bash

set -euo pipefail

bioproject="$1"
sample="$2"

echo "Counting child taxids for $bioproject $sample"

aws s3 cp "s3://nao-mgs/$bioproject/allcounts/$sample.tsv.gz" - | gunzip | \
    ./count_child_nodes.py | gzip | \
    aws s3 cp - "s3://nao-mgs/$bioproject/childcounts/$sample.tsv.gz"
