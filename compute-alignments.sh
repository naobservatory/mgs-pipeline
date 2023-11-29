#!/usr/bin/env bash

# I wish I could do this in pure python, which would have more robust error
# handling, but I'm not sure how to write this kind of branching command tree
# with the subprocess module.

set -e
set +x

CMD="/home/ec2-user/bowtie2-2.5.2-linux-x86_64/bowtie2"
CMD+=" --local"
CMD+=" -x /home/ec2-user/mgs-pipeline/bowtie/human-viruses"
CMD+=" --threads 28"
CMD+=" --very-sensitive-local"
CMD+=" --score-min G,1,0"
CMD+=" --mp 4,1"
CMD+=" --no-unal"
CMD+=" --no-sq"
CMD+=" -S $1"

if [ $# = 2 ]; then
    $CMD -U <(aws s3 cp $2 - | gunzip)
elif [ $# = 3 ]; then
    $CMD -1 <(aws s3 cp $2 - | gunzip) -2 <(aws s3 cp $3 - | gunzip)
else
    >&2 echo "usage: $0 <sam_out> <fastq_in> [fastq2_in]"
fi
