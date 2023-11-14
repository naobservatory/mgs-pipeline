#!/usr/bin/env python
import pysam
import sys
import csv
import os
from math import sqrt, log, e

hvsams_dir, alignment_scores_out = sys.argv[1:]


def sam_records():
    with open(alignment_scores_out, "w") as out_f:
        out_f = csv.writer(out_f, delimiter="\t")
        headers = [
            "read_id",
            "sequence",
            "alignment_score",
            "cigar",
            "reference_start",
            "sqrt_length_adj_score",
            "log_length_adj_score",
            "flag",
            "read_type",
            "read_length",
            "reference_id",
        ]
        out_f.writerow((headers))

        for filename in os.listdir(hvsams_dir):
            sam_file = pysam.AlignmentFile(hvsams_dir + "/" + filename, "r")
            for read in sam_file.fetch():
                read_id = read.query_name
                sequence = read.query_sequence
                cigar = read.cigarstring
                reference_start = read.reference_start
                reference_id = sam_file.get_reference_name(read.reference_id)
                flag = (
                    read.flag
                )  # https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#:~:text=all%20applicable%20flags.-,Flags,-relevant%20to%20Bowtie
                if flag < 64:
                    read_type = "combined"
                elif flag > 64 and flag < 128:
                    read_type = "read_1"
                elif flag > 128:
                    read_type = "read_2"
                read_length = int(read.query_length)
                log_read_length = log(read_length, e)
                sqrt_read_length = sqrt(read_length)
                try:
                    alignment_score = read.get_tag("AS")
                    sqrt_length_adj_score = alignment_score / sqrt_read_length
                    log_length_adj_score = alignment_score / log_read_length
                except:
                    alignment_score = 0
                    sqrt_length_adj_score = 0
                    log_length_adj_score = 0
                row = [
                    str(read_id),
                    str(sequence),
                    int(alignment_score),
                    str(cigar),
                    int(reference_start),
                    str(reference_id),
                    float(sqrt_length_adj_score),
                    float(log_length_adj_score),
                    str(flag),
                    str(read_type),
                    int(read_length),
                ]
                out_f.writerow(row)


def start():
    sam_records()


if __name__ == "__main__":
    start()
