#!/usr/bin/env python3

import sys
import glob
import os.path
import subprocess
from collections import defaultdict

bowtie_dir, data_dir = sys.argv[1:]

samples = defaultdict(dict)
for fastq in glob.glob(f"{data_dir}/hvfastqs/*.fastq"):
    sample, category, _ = os.path.basename(fastq).rsplit(".", 2)
    samples[sample][category] = fastq
for sample, fastqs in sorted(samples.items()):
    out = f"{data_dir}/hvsams/%s.sam" % sample
    if os.path.exists(out):
        print("Skipping %s" % sample)
        continue

    print("Running bowtie2 on %s" % sample)
    cmd = [
        f"{bowtie_dir}/bowtie2-2.5.1-macos-arm64/bowtie2",
        "--local",
        "-x",
        "human-viruses",
        "--very-sensitive-local",
        "--score-min",
        "G,1,0",
        "--mp",
        "2,0",
        "--threads",
        "24",
        "-S",
        out,
    ]
    assert ("pair1" in fastqs) == ("pair2" in fastqs)
    if "pair1" in fastqs:
        cmd.extend(["-1", fastqs["pair1"], "-2", fastqs["pair2"]])
    if "combined" in fastqs:
        cmd.extend(["-U", fastqs["combined"]])

    try:
        subprocess.check_call(cmd)
    except Exception:
        print(" ".join(cmd))
        raise
