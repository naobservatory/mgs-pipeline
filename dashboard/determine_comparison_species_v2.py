#!/usr/bin/env python3

import sys
from collections import Counter, defaultdict

VIRUS = 10239

children = defaultdict(list)  # parent_taxid -> [children]
with open("nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = line.replace("\t|\n", "").split(
            "\t|\t"
        )
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)
        if child_taxid != parent_taxid:
            children[parent_taxid].append(child_taxid)


def populate_set(taxid, s):
    s.add(taxid)
    for child in children[taxid]:
        populate_set(child, s)


viruses = set()
populate_set(VIRUS, viruses)

all_counts = Counter()
viral_counts = Counter()

col = None
for line in sys.stdin:
    row = line.rstrip("\n").split("\t")
    if not col:
        col = row
        continue

    taxid = int(row[col.index("taxid")])
    direct_assignments = int(row[col.index("n_reads_direct")])

    all_counts[taxid] += direct_assignments
    if taxid in viruses:
        viral_counts[taxid] += direct_assignments

for count, taxid in sorted(
    [
        (count, taxid)
        for (taxid, count) in (
            set(all_counts.most_common(10)) | set(viral_counts.most_common(10))
        )
    ],
    reverse=True,
):
    # We take the top 30 now, because with the v2 pipeline it's per delivery
    # instead of per sample.
    if count > 30:
        print("%s\t%s" % (count, taxid))
