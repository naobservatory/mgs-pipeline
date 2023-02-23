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

## Design

### Data storage

At each stage data will be stored in an S3 bucket.  The structure is:

    [bucket]/
      [study id]/
         metadata.tsv
         raw/
         cleaned/
         processed/

In cases where the data comes from the SRA, the study ID is the SRA accession.
For example, "PRJNA729801" for the Rothman 2021 data.

The `metadata.tsv` file is tab-separated, where the first column is the sample
id and the remaining columns are study-specific.  For data we downloaded from
the SRA the sample ID is the SRA accession.  For example, "SRR14530767" for the
2020-08-11 HTP sample in the Rothman 2021 data.  The scripts used to create the
metadata file should be checked into this repo under `studies/[accession]/`.

Files under `raw/` have the contents as we received them, but have been renamed
to `[sampleID].fasta.gz`.  For paired-end data there will be two files,
`[sampleID]_1.fasta.gz` and `[sampleID]_2.fasta.gz`.

Files under `cleaned/` are the output of AdapterRemoval2; see below.

Files under `processed/` include our quality control results, species
classification, and any other processing we run.

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
