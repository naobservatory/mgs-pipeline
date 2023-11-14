#!/usr/bin/env bash

set -e # exit on error 

ROOT_DIR="$PWD"
DATA_DIR="$ROOT_DIR/validation-test"
BOWTIE_DIR="$ROOT_DIR/.." #or just add bowtie2 to your path

if [ ! -e $DATA_DIR/all_references.fna ]; then
    find raw-genomes | grep .fna$ | xargs cat > $DATA_DIR/all_references.fna
fi


if [ ! -e $ROOT_DIR/human-viruses.1.bt2 ]; then
        $BOWTIE_DIR/bowtie2-2.5.1-macos-arm64/bowtie2-build \
	-f \
        --threads 32 \
        --verbose \
        $DATA_DIR/all_references.fna \
	$ROOT_DIR/human-viruses
fi

if [ ! -e $DATA_DIR/hvsams ]; then
    mkdir $DATA_DIR/hvsams
    $ROOT_DIR/run_bowtie2.py \
	$BOWTIE_DIR
fi

if [ ! -e alignment_scores.tsv ]; then
    $ROOT_DIR/alignment_scores.py \
	$ROOT_DIR/hvsams \
	alignment_scores.tsv
fi


