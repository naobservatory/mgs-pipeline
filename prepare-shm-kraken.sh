#!/usr/bin/env bash

set -e
set -u

if [ -d /dev/shm/kraken-db/ ]; then
    echo "DB already loaded"
    exit 0
fi

MEMORY_SIZE=$(cat /proc/meminfo | \
                  grep MemTotal | \
                  awk '{print  $2}')

if [[ $MEMORY_SIZE -lt 128000000 ]]; then
    echo "Insufficient memory; need c6a.16xlarge or bigger"
    exit 1
fi

mkdir /dev/shm/kraken-db/
aws s3 cp s3://genome-idx/kraken/k2_standard_20240605.tar.gz - | \
    tar -xzvv -C /dev/shm/kraken-db
