#!/usr/bin/env python3

import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator

reads = 0
for title, sequence, quality in FastqGeneralIterator(sys.stdin):
    reads += 1

print(reads)
