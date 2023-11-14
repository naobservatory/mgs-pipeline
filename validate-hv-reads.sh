#!/usr/bin/env bash

set -e # exit on error 

ROOT_DIR="$PWD"
DATA_DIR="$ROOT_DIR/validation-test"
BOWTIE_DIR="$ROOT_DIR/.." #or just add bowtie2 to your path


if [ ! -d $DATA_DIR/hvfastqs ]; then
    mkdir $DATA_DIR/hvfastqs
    ls $DATA_DIR/hvreads | \
        xargs -P 32 -I {} \
	$ROOT_DIR/json_to_fasta.py $DATA_DIR/hvreads $DATA_DIR/hvfastqs {}
fi

if [ ! -e observed-human-virus-taxids.txt ]; then
    $ROOT_DIR/determine_hv_taxids.py \
        $DATA_DIR/hvreads/ \
        $ROOT_DIR/human-viruses.tsv \
        $ROOT_DIR/observed-human-virus-taxids.txt
fi

$ROOT_DIR/get_genomes.py

if [ ! -d $ROOT_DIR/raw-genomes ]; then
    mkdir $ROOT_DIR/raw-genomes
    for x in $(find $ROOT_DIR/refseq/ | grep gz$); do
        gunzip -c "$x" > $ROOT_DIR/raw-genomes/$(basename ${x/.fna.gz/.fna})
    done
fi

if [ ! -e $ROOT_DIR/all_references.fna ]; then
    find raw-genomes | grep .fna$ | xargs cat > $ROOT_DIR/all_references.fna
fi


if [ ! -e $ROOT_DIR/human-viruses.1.bt2 ]; then
        $BOWTIE_DIR/bowtie2-2.5.1-macos-arm64/bowtie2-build \
	-f \
        --threads 32 \
        --verbose \
        $ROOT_DIR/all_references.fna \
	$ROOT_DIR/human-viruses
fi

if [ ! -e $DATA_DIR/hvsams ]; then
    mkdir $DATA_DIR/hvsams
    $ROOT_DIR/run_bowtie2.py \
	$BOWTIE_DIR \
	$DATA_DIR
fi

if [ ! -e alignment_scores.tsv ]; then
    $ROOT_DIR/return_alignment_scores.py \
	$DATA_DIR/hvsams \
	alignment_scores.tsv
fi


