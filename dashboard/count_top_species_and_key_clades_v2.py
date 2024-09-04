#!/usr/bin/env python3

import sys
import json
from collections import defaultdict, Counter

taxids = set()
with open("key_clade_taxids.txt") as inf:
    for line in inf:
        line = line.strip()
        if line:
            taxids.add(int(line))

with open("comparison_taxids_v2.txt") as inf:
    for line in inf:
        line = line.strip()
        if line:
            taxids.add(int(line))

# sample -> taxid -> count
comparisons = defaultdict(Counter)
col = None
for line in sys.stdin:
    row = line.rstrip("\n").split("\t")
    if not col:
        col = row
        continue

    sample = row[col.index("sample")]
    taxid = int(row[col.index("taxid")])
    clade_assignments = int(row[col.index("n_reads_clade")])

    if taxid in taxids:
        comparisons[sample][taxid] += clade_assignments

print(json.dumps(comparisons))
