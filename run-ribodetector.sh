#!/bin/bash

# Set constants
# If two input/output files are given, RiboDetector assumes they are paired-end format
FASTQ_PATH="s3://nao-mgs/PRJNA681566/cleaned/SRR13167434.collapsed.gz"
OUTPUT_PATH="s3://lenni-analysis/ribodetector"
TEMP_FASTQ_GZ="/tmp/SRR13167434.collapsed.gz"
TEMP_FASTQ="/tmp/SRR13167434.collapsed.fq"
TEMP_OUTPUT="/tmp/SRR13167434.collapsed.nonrrna.fq"
TEMP_RRNA_OUTPUT="/tmp/SRR13167434.collapsed.rrna.fq"
TEMP_NON_RRNA_IDS="/tmp/non_rrna_ids.txt"

# Stream the FASTQ data from S3 to a temp gzipped file
aws s3 cp $FASTQ_PATH $TEMP_FASTQ_GZ

# Decompress the FASTQ file
gunzip -c $TEMP_FASTQ_GZ > $TEMP_FASTQ

# Replace spaces with underscores in FASTQ headers
sed -i 's/ /_/' $TEMP_FASTQ

# Calculate the average read length
AVG_READ_LENGTH=$(awk 'NR%4==2 {sum+=length; count++} END {print int(sum/count)}' $TEMP_FASTQ)

# Execute ribodetector_cpu
ribodetector_cpu -t 20 \
  -l $AVG_READ_LENGTH \
  -i $TEMP_FASTQ \
  -e rrna \
  --chunk_size 256 \
  -o $TEMP_OUTPUT

# Extract non-rRNA IDs
grep "^@" $TEMP_OUTPUT | sed 's/^@//' | cut -d' ' -f1 > $TEMP_NON_RRNA_IDS

# Extract all IDs
grep "^@" $TEMP_FASTQ | sed 's/^@//' > /tmp/all_ids.txt

# Get the rRNA IDs (those not in the non-rRNA list)
comm -23 <(sort /tmp/all_ids.txt) <(sort $TEMP_NON_RRNA_IDS) > /tmp/rrna_ids.txt

# Extract rRNA reads using seqtk
/home/ec2-user/seqtk/seqtk subseq $TEMP_FASTQ /tmp/rrna_ids.txt > $TEMP_RRNA_OUTPUT 2>/dev/null

# Send the outputs to S3
aws s3 cp $TEMP_OUTPUT $OUTPUT_PATH/SRR13167434.collapsed.nonrrna.fq
aws s3 cp $TEMP_RRNA_OUTPUT $OUTPUT_PATH/SRR13167434.collapsed.rrna.fq

# Count the number of reads in both files and print the difference
ORIGINAL_READS=$(cat $TEMP_FASTQ | wc -l | awk '{print $1/4}')
FILTERED_READS=$(cat $TEMP_OUTPUT | wc -l | awk '{print $1/4}')

DIFFERENCE=$(($ORIGINAL_READS - $FILTERED_READS))
PERCENTAGE_REMOVED=$(echo "scale=2; ($DIFFERENCE / $ORIGINAL_READS) * 100" | bc)

echo "Original number of reads: $ORIGINAL_READS"
echo "Number of reads after ribodetector: $FILTERED_READS"
echo "Difference (number of reads removed): $DIFFERENCE ($PERCENTAGE_REMOVED% removed)"

# Optional: Cleanup temp files (if desired)
rm $TEMP_FASTQ_GZ
rm $TEMP_FASTQ
rm $TEMP_OUTPUT
rm $TEMP_RRNA_OUTPUT
rm $TEMP_NON_RRNA_IDS

