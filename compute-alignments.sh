#!/usr/bin/env bash

# I wish I could do this in pure python, which would have more robust error
# handling, but I'm not sure how to write this kind of branching command tree
# with the subprocess module.

set -e
set +x

if [ $# != 3 ] && [ $# != 4 ]; then
    >&2 echo "usage: $0 <human|hv|hb> <sam_out> <fastq_in> [fastq2_in]"
    exit 1
fi

DB="$1"
shift

OUT="$1"
shift

CMD="/home/ec2-user/bowtie2-2.5.2-linux-x86_64/bowtie2"
CMD+=" --threads 8"
CMD+=" --no-unal"
CMD+=" --no-sq"
CMD+=" -S $OUT"

if [ "$DB" = "human" ]; then
    # Pre-built human DB
    CMD+=" -x /home/ec2-user/mgs-pipeline/bowtie/chm13.draft_v1.0_plusY"
    # When identifying human reads use default tuning settings.
elif [ "$DB" = "hv" ] || [ "$DB" = "hb" ]; then
    # Custom-built HV DB
    if [ "$DB" = "hv" ]; then
        CMD+=" -x /home/ec2-user/mgs-pipeline/bowtie/human-viruses"
    else
        CMD+=" -x /home/ec2-user/mgs-pipeline/bowtie/human-respiratory-bacteria"
    fi
    # When identifying HV reads use looser settings and filter more later.
    CMD+=" --local"
    CMD+=" --very-sensitive-local"
    CMD+=" --score-min G,1,0"
    CMD+=" --mp 4,1"
else
    >&2 echo "DB must be 'human', 'hv', or 'hb'; got '$DB'"
    exit 1
fi

if [ "$DB" = "hb" ]; then
    # The hb data is stored locally, since we only run this on a small portion
    # of reads.
    if [ $# = 1 ]; then
        $CMD -U $1
    else
        $CMD -1 $1 -2 $2
    fi
else
    # Otherwise, we run this on all reads as we stream them from S3.
    if [ $# = 1 ]; then
        $CMD -U <(aws s3 cp $1 - | gunzip)
    else
        $CMD -1 <(aws s3 cp $1 - | gunzip) -2 <(aws s3 cp $2 - | gunzip)
    fi
fi
