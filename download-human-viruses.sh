#!/usr/bin/env bash

# TSV format is "ncbi taxid \t scientific name"
# 9606 is the taxid for humans
curl -sS https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv \
    | awk -F'\t' '$8=="9606"{print $1"\t"$2}' \
    | sort -n \
    > human-viruses-raw.tsv
./expand-human-viruses.py human-viruses-raw.tsv human-viruses.tsv
