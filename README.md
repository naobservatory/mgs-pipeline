# Metagenomic Sequencing Pipeline

When we work with metagenomic sequencing data we have some common needs and
questions.  This pipeline will let us run all our incoming data through a
consistent set of steps, answer common questions, and prepare cleaned and
processed data for further analysis.

Questions we'll want it to answer include:

* Data quality:
  * What quality scores do we see?
  * Is there evidence of PCR duplicate issues?
  * What is the distribution of insert lengths?
  * Does the level of quality dropoff along reads make sense?
* Data contents:
  * What are the relative abundances of species of interest?
  * What sort of cross-sample variation is there?

## Status

As of 2023-02-27 handles the Rothman 2021 data up through cleaning.  Todo:

* Quality control
* Species classification
* Other data sources

## Design

### Data storage

At each stage data will be stored in an S3 bucket.  The structure is:

    s3://nao-mgs/
      [studyId]/
         raw/
         cleaned/
         processed/

In cases where the data comes from the SRA, the study ID is the SRA accession.
For example, "PRJNA729801" for the Rothman 2021 data.

Files under `raw/` have the contents as we received them, but have been renamed
to `[sampleID].fasta.gz`.  For paired-end data there will be two files,
`[sampleID]_1.fasta.gz` and `[sampleID]_2.fasta.gz`.

Files under `cleaned/` are the output of AdapterRemoval2; see below.

Files under `processed/` include our quality control results, species
classification, and any other processing we run.

### Metadata

Metadata goes in this repo under `studies/[accession]/metadata/`.  This
includes both the metadata file and the scripts that prepare it.

Each study has a `studies/[accession]/metadata/metadata.tsv` where the first
column is the sample ID and the remaining columns are study-specific.  For data
we downloaded from the SRA the sample ID is the SRA accession.  For example,
"SRR14530767" for the 2020-08-11 HTP sample in the Rothman 2021 data.

### Input

FASTQ files, in pairs if we have paired-end sequencing.  If we get data in
other formats, however, we'll include the conversion as a pre-processing step.

### Cleaning

AdapterRemoval2 to trim adapters, remove low-quality bases from the
ends of reads, and collapse overlapping paired-end reads.

### Quality Control

FastQC, plus additional custom checks for the issues we're specifically
concerned about.

### Species Classification

Kraken to assign taxonomic identifiers to reads, and Bracken to estimate
per-species relative abundance.  We'll may need to build our own database or
extend an existing one to ensure we have the level of coverage we need for
human viruses.

When we do want to build our own DB kraken2-build has good docs.

## Dependencies

### AdapterRemoval2

Visit https://github.com/MikkelSchubert/adapterremoval/ to get the latest
version number, and then:

```
wget -O adapterremoval-2.3.1.tar.gz \
     https://github.com/MikkelSchubert/adapterremoval/archive/v2.3.1.tar.gz
tar xvzf adapterremoval-2.3.1.tar.gz
cd adapterremoval-2.3.1
make
sudo make install
```

### Kraken

#### Install

```
git clone git@github.com:DerrickWood/kraken2.git
cd kraken2/
./install_kraken2.sh ~/kraken2-install
```

#### Set up DB

For now I'm using the Standard pre-built DB, with a 16GB cap.  Longer term what
we want depends on how much memory we're ok using; full Standard doesn't quite
fit on a machine with 64GB of RAM.

```
mkdir ~/kraken-db
cd ~/kraken-db
# avoids needing 2x the storage
aws s3 cp s3://genome-idx/kraken/k2_standard_16gb_20221209.tar.gz - | tar -xvz
```
