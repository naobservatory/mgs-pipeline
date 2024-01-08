#!/usr/bin/env python3

from collections import defaultdict

UNASSIGNED = 0
ROOT = 1
BACTERIA = 2
VIRUS = 10239
# Per
# https://wwwnc.cdc.gov/travel/yellowbook/2024/posttravel-evaluation/respiratory-infections,
# plus TB
RESPIRATORY_BACTERIA = [
    28450, # Burkholderia pseudomallei
    520,   # Bordetella pertussis
    83558, # Chlamydia pneumoniae
    1717,  # Corynebacterium diphtheriae
    727,   # Haemophilus influenzae
    2104,  # Mycoplasmoides pneumoniaea
    1313,  # Streptococcus pneumoniae
    1773,  # Mycobacterium tuberculosis
]

children = defaultdict(set)  # parent_taxid -> [children]
with open("nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = line.replace("\t|\n", "").split(
            "\t|\t"
        )
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)
        if child_taxid != parent_taxid:
            children[parent_taxid].add(child_taxid)

key_clades = set([UNASSIGNED, ROOT, BACTERIA, VIRUS,
                  *RESPIRATORY_BACTERIA])

key_clades |= children[ROOT]
key_clades |= children[BACTERIA]
key_clades |= children[VIRUS]

for key_clade in sorted(key_clades):
    print(key_clade)
