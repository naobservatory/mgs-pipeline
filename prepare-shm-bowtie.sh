#!/usr/bin/env bash

set -e
set -u

if [ -d /dev/shm/bowtie-db/ ]; then
    echo "DB already loaded"
    exit 0
fi

mkdir /dev/shm/bowtie-db/
cp bowtie/chm13.draft_v1.0_plusY* /dev/shm/bowtie-db/
cp bowtie/human-viruses* /dev/shm/bowtie-db/
