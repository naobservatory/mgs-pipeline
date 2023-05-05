#!/usr/bin/env python3

import sys
import json

taxids = set()
with open("key_clade_taxids.txt") as inf:
    for line in inf:
        line = line.strip()
        if line:
            taxids.add(int(line))
            
with open("comparison_taxids.txt") as inf:
    for line in inf:
        line = line.strip()
        if line:
            taxids.add(int(line))

comparisons = {}
for line in sys.stdin:
    taxid, direct_assignments, direct_hits, \
        clade_assignments, clade_hits = line.strip().split("\t")
    taxid = int(taxid)
    clade_assignments = int(clade_assignments)

    if taxid in taxids:
        comparisons[taxid] = clade_assignments

print(json.dumps(comparisons))
