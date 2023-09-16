#!/bin/bash

# Set constants
FASTQ_PATH="s3://nao-mgs/PRJNA681566/cleaned/SRR13167434.collapsed.gz"
RRNA_OUTPUT_PATH="s3://lenni-analysis/ribodetector/SRR13167434.collapsed.rrna.fq"
NONRRNA_OUTPUT_PATH="s3://lenni-analysis/ribodetector/SRR13167434.collapsed.nonrrna.fq"

TEMP_FASTQ="/tmp/SRR13167434.collapsed.fq"
TEMP_RRNA_OUTPUT="/tmp/SRR13167434.collapsed.rrna.fq"
TEMP_NONRRNA_OUTPUT="/tmp/SRR13167434.collapsed.nonrrna.fq"
TEMP_RRNA_TEST_OUTPUT="/tmp/SRR13167434.collapsed.rrna.test.fq"

# Download files from S3
aws s3 cp $FASTQ_PATH - | gunzip -c > $TEMP_FASTQ
aws s3 cp $RRNA_OUTPUT_PATH $TEMP_RRNA_OUTPUT
aws s3 cp $NONRRNA_OUTPUT_PATH $TEMP_NONRRNA_OUTPUT

# Test 1: Check if total reads match
ORIGINAL_READS=$(cat $TEMP_FASTQ | wc -l | awk '{print $1/4}')
RRNA_READS=$(cat $TEMP_RRNA_OUTPUT | wc -l | awk '{print $1/4}')
NONRRNA_READS=$(cat $TEMP_NONRRNA_OUTPUT | wc -l | awk '{print $1/4}')

TOTAL_OUTPUT_READS=$(($RRNA_READS + $NONRRNA_READS))

if [ $ORIGINAL_READS -eq $TOTAL_OUTPUT_READS ]; then
    echo "Test 1 Passed: Total reads match"
else
    echo "Test 1 Failed: Total reads do not match"
    echo "Original Reads: $ORIGINAL_READS, Total Output Reads: $TOTAL_OUTPUT_READS"
fi

# Test 2: Check if there is no overlap in reads
OVERLAP=$(grep "^@" $TEMP_RRNA_OUTPUT | sort | comm -12 - <(grep "^@" $TEMP_NONRRNA_OUTPUT | sort) | wc -l)

if [ $OVERLAP -eq 0 ]; then
    echo "Test 2 Passed: No overlapping reads found"
else
    echo "Test 2 Failed: Overlapping reads found"
fi

# Test 3: Run the rRNA reads back through ribodetector and verify classification

# Calculate the average read length from the original FASTQ
AVG_READ_LENGTH=$(awk 'NR%4==2 {sum+=length; count++} END {print int(sum/count)}' $TEMP_FASTQ)

ribodetector_cpu -t 20 \
  -l $AVG_READ_LENGTH \
  -i $TEMP_RRNA_OUTPUT \
  -e rrna \
  --chunk_size 256 \
  -o $TEMP_RRNA_TEST_OUTPUT

# Check if the rRNA test output file is empty
RRNA_TEST_READS=$(cat $TEMP_RRNA_TEST_OUTPUT | wc -l | awk '{print $1/4}')

if [ $RRNA_TEST_READS -eq 0 ]; then
    echo "Test 3 Passed: All rRNA reads were correctly classified as rRNA by ribodetector in the second pass"
else
    echo "Test 3 Failed: Some rRNA reads were not classified as rRNA by ribodetector in the second pass"
    echo "Number of Reads Incorrectly Classified: $RRNA_TEST_READS"
fi

# Cleanup temp files
rm $TEMP_FASTQ
rm $TEMP_RRNA_OUTPUT
rm $TEMP_NONRRNA_OUTPUT
rm $TEMP_RRNA_TEST_OUTPUT


