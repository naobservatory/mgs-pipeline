#!/usr/bin/env python3

# input: tsv of taxid, direct assignments, direct hits
# output: tsv of taxid, direct assignments, direct hits, all assignments, all hits

import sys

parents = {}  # child_taxid -> parent_taxid
children = {} # parent_taxid -> [child_taxid]
with open("dashboard/nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)

        if parent_taxid not in children:
            children[parent_taxid] = []
        if parent_taxid != child_taxid:
            children[parent_taxid].append(child_taxid)
        parents[child_taxid] = parent_taxid

direct_assignments = {} # taxid -> assignments
direct_hits = {} # taxid -> hits

for line in sys.stdin:
    line = line.strip()
    if not line: continue

    taxid, n_direct_assignments, n_direct_hits = [
        int(x) for x in line.split("\t")]

    direct_assignments[taxid] = n_direct_assignments
    direct_hits[taxid] = n_direct_hits

all_assignments = {}
all_hits = {}

# There are smart algorithms you could use here, which elegantly propagate
# counts up the tree.  We're not goint to do that.  Just:
#   1. Determine all nodes that have counts or whose children have any counts
#   2. For each of these nodes, count it and its children recursively

all_counted_taxids = set()
for taxid in direct_hits:
    while True:
        all_counted_taxids.add(taxid)
        if taxid in [0, 1]:
            break
        taxid = parents[taxid]

def count_children(target_taxid, current_taxid):
    all_assignments[target_taxid] += direct_assignments.get(current_taxid, 0)
    all_hits[target_taxid] += direct_hits.get(current_taxid, 0)

    for taxid in children.get(current_taxid, []):
        count_children(target_taxid, taxid)            
        
for target_taxid in all_counted_taxids:
    all_assignments[target_taxid] = all_hits[target_taxid] = 0
    count_children(target_taxid, target_taxid)

for taxid in sorted(all_counted_taxids):
    print("%s\t%s\t%s\t%s\t%s" % (
        taxid,
        direct_assignments.get(taxid, 0),
        direct_hits.get(taxid, 0),
        all_assignments.get(taxid, 0),
        all_hits.get(taxid, 0)))
