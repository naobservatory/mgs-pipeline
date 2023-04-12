#!/usr/bin/env python3

# input: kraken output
# output: tsv of taxid, assignments, hits

import re
import sys
from collections import Counter

assignments = Counter() # taxid -> assignments
hits = Counter() # taxid -> hits
for line in sys.stdin:
    line = line.strip()
    if not line: continue
    _, _, name_and_taxid, _, encoded_hits = line.split("\t")

    taxid, = re.findall("^.*[(]taxid ([0-9]+)[)]$", name_and_taxid)
    taxid = int(taxid)
    assignments[taxid] += 1
    
    for hit in set(int(x) for x in re.findall("([0-9]+):", encoded_hits)):
        hits[hit] += 1

all_taxids = set(hits).union(assignments)
        
for taxid in sorted(all_taxids):
    print("%s\t%s\t%s" % (
        taxid, assignments[taxid], hits[taxid]))
        
        
