#!/usr/bin/env bash

# ncbi taxid TAB scientific name
curl -sS https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv \
    | awk -F'\t' '$8=="9606"{print $1"\t"$2}' \
    | sort -n \
    > human-viruses.tsv
