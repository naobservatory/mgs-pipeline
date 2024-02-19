#!/usr/bin/env bash

# I wish I could do this in pure python, which would have more robust error
# handling, but I'm not sure how to write this kind of branching command tree
# with the subprocess module.

set -e
set -u

if [ $# != 4 ] && [ $# != 5 ]; then
    >&2 echo "usage: $0 <disk|mmap> <'human'|'hv'> <sam_out> <fastq_in> [fastq2_in]"
    exit 1
fi

DISK_MEMORY="$1"
shift

DB="$1"
shift

OUT="$1"
shift

COPY_LOCAL=true

CMD="/home/ec2-user/bowtie2-2.5.2-linux-x86_64/bowtie2"
if [ "$DISK_MEMORY" = "disk" ]; then
    DB_DIR="/home/ec2-user/mgs-pipeline/bowtie"
    CMD+=" --threads 8"
elif [ "$DISK_MEMORY" = "mmap" ]; then
    DB_DIR="/dev/shm/bowtie-db"
    CMD+=" --threads 4"
    CMD+=" --mm"
else
    echo "Bad argument $DISK_MEMORY; expected 'disk' or 'mmap'"
    exit 1
fi
CMD+=" --no-unal"
CMD+=" --no-sq"
CMD+=" -S $OUT"

if [ "$DB" = "human" ]; then
    # Pre-built human DB
    CMD+=" -x $DB_DIR/chm13.draft_v1.0_plusY"
    # When identifying human reads use default tuning settings.
elif [ "$DB" = "hv" ]; then
    # Custom-built HV DB
    CMD+=" -x $DB_DIR/human-viruses"
    # When identifying HV reads use looser settings and filter more later.
    CMD+=" --local"
    CMD+=" --very-sensitive-local"
    CMD+=" --score-min G,1,0"
    CMD+=" --mp 4,1"
else
    >&2 echo "DB must be either 'human' or 'hv'; got '$DB'"
    exit 1
fi

if [ $# = 1 ]; then
    if $COPY_LOCAL; then
        aws s3 cp $1 tmp.1.fastq.gz
        $CMD -U <(cat tmp.1.fastq.gz | gunzip)
        rm tmp.1.fastq.gz
    else
        $CMD -U <(aws s3 cp $1 - | gunzip)
    fi
else
    if true; then
        aws s3 cp $1 tmp.1.fastq.gz
        aws s3 cp $2 tmp.2.fastq.gz
        $CMD -1 <(cat tmp.1.fastq.gz | gunzip) -2 <(cat tmp.2.fastq.gz | gunzip)
        rm tmp.1.fastq.gz
        rm tmp.2.fastq.gz
    else
        $CMD -1 <(aws s3 cp $1 - | gunzip) -2 <(aws s3 cp $2 - | gunzip)
    fi
fi
