#!/usr/bin/bash

# usage: import-sratools.sh bioproject

BIOPROJECT="$1"

cat bioprojects/$BIOPROJECT/metadata/metadata.tsv | \
    awk -F '\t' '{print $1}' | \
    xargs -P 4 -I {} ./import-sratools-single.sh $BIOPROJECT {}

rm -r fasterq.tmp.assembly.*
