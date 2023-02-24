#!/usr/bin/env bash

if [ $# -ne 1 ]; then
    echo "Usage: ./run.sh <study>" > /dev/stderr
    exit 1
fi

STUDY="$1"

