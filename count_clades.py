#!/usr/bin/env python3

# input: kraken output
# output: tsv of taxid, assignments, hits

import re
import sys
from collections import Counter

parents = {}  # child_taxid -> parent_taxid
with open("dashboard/nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)
        parents[child_taxid] = parent_taxid

direct_assignments = Counter() # taxid -> direct assignments
direct_hits = Counter()        # taxid -> direct hits
clade_assignments = Counter()  # taxid -> clade assignments
clade_hits = Counter()         # taxid -> clade hits

for line in sys.stdin:
    line = line.strip()
    if not line: continue
    _, _, name_and_taxid, _, encoded_hits = line.split("\t")

    taxid, = re.findall("^.*[(]taxid ([0-9]+)[)]$", name_and_taxid)
    taxid = int(taxid)
    direct_assignments[taxid] += 1
    while True:
        clade_assignments[taxid] += 1
        if taxid in [0, 1]:
            break
        taxid = parents[taxid]

    direct_incremented = set()
    clade_incremented = set()
    for hit in re.findall("([0-9]+):", encoded_hits):
        hit = int(hit)
        if hit not in direct_incremented:
            direct_hits[hit] += 1
            direct_incremented.add(hit)

        while hit not in clade_incremented:
            clade_hits[hit] += 1
            clade_incremented.add(hit)

            if hit in [0, 1]:
                break
            hit = parents[hit]

for taxid in sorted(clade_hits):
    print("%s\t%s\t%s\t%s\t%s" % (
        taxid,
        direct_assignments[taxid], direct_hits[taxid],
        clade_assignments[taxid], clade_hits[taxid]))

    
