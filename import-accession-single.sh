#!/usr/bin/env bash

BIOPROJECT="$1"
URL="$2"
BASE=$(basename "$URL")

if aws s3 ls "s3://nao-mgs/$BIOPROJECT/raw/$BASE" | \
        awk '$3>0{print}' | \
        grep "$BASE" > /dev/null; then
    echo "$BASE already present"
else
    echo "Downloading $BASE"
    curl -sS "$URL" | aws s3 cp - "s3://nao-mgs/$BIOPROJECT/raw/$BASE"
fi
