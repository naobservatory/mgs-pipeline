#!/usr/bin/env python3
import os
import sys
import json

in_dir, out_dir, in_fname = sys.argv[1:]
sample = os.path.basename(in_fname).replace(".hvreads.json", "")

pair1 = []
pair2 = []
combined = []

with open(os.path.join(in_dir, in_fname)) as inf:
    for seq_id, details in sorted(json.load(inf).items()):
        kraken_assignment, kraken_info, *reads = details
        assert len(reads) in [1, 2]
        if len(reads) == 1:
            (read,) = reads
            combined.append((seq_id, read))
        else:
            read1, read2 = reads
            pair1.append((seq_id, read1))
            pair2.append((seq_id, read2))

for label, reads in [
    ("pair1", pair1),
    ("pair2", pair2),
    ("combined", combined),
]:
    with open(os.path.join(out_dir, sample + "." + label + ".fastq"), "w") as outf:
        for seq_id, read in reads:
            seq, quality = read
            outf.write("@%s\n%s\n+\n%s\n" % (seq_id, seq, quality))
