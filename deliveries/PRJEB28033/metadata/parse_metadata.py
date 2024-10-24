#!/usr/bin/env python3

# Run prepare_metadata.sh, which calls this script.

import sys


def start(raw_metadata_in, common_metadata_out):
    data = []
    with open(raw_metadata_in) as inf:
        for line in inf:
            line = line.strip()

            accession, fastq_ftps, sampleid = line.split("\t")

            if fastq_ftps == "fastq_ftp":
                continue  # skip header line

            sampleid = "NYC-%s" % (sampleid.split("-")[-1].zfill(2))

            data.append([accession, sampleid])

    with open(common_metadata_out, "w") as outf:
        for accession, sampleid in data:
            outf.write("\t".join([accession, sampleid]) + "\n")


if __name__ == "__main__":
    start(*sys.argv[1:])
