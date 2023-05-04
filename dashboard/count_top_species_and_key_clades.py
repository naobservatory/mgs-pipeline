#!/usr/bin/env python3

import sys
import json

key_clades = {}
with open("key_clade_taxids.txt") as inf:
    for line in inf:
        line = line.strip()
        if line:
            key_clades[int(line)] = 0
            
comparisons = {}
with open("comparison_taxids.txt") as inf:
    for line in inf:
        line = line.strip()
        if line:
            comparisons[int(line)] = 0

for line in sys.stdin:
    taxid, direct_assignments, direct_hits, \
        clade_assignments, clade_hits = line.strip().split("\t")
    taxid = int(taxid)
    direct_assignments = int(direct_assignments)
    clade_assignments = int(clade_assignments)

    if taxid in key_clades:
        key_clades[taxid] = clade_assignments
    if taxid in comparisons:
        comparisons[taxid] = direct_assignments

print(json.dumps([key_clades, comparisons]))
