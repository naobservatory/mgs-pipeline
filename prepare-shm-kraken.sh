#!/usr/bin/env bash

set -e
set -u

if [ -d /dev/shm/kraken-db/ ]; then
    echo "DB already loaded"
    exit 0
fi

mkdir /dev/shm/kraken-db/
aws s3 cp s3://genome-idx/kraken/k2_standard_20240605.tar.gz - | \
    tar -xzv -C /dev/shm/kraken-db
