#!/usr/bin/env python3

import sys
import glob
import json

hvreads_dir, human_viruses, taxids_out = sys.argv[1:]

hv_taxids = set()  # strings
seen_hv_taxids = set()  # ints

with open(human_viruses) as inf:
    for line in inf:
        taxid = line.split("\t")[0]
        hv_taxids.add(taxid)
print(hvreads_dir)
print(human_viruses)
for fname in glob.glob("%s/*.json" % hvreads_dir):
    print(fname)    
    with open(fname) as inf:
        for _, kraken_hits, *_ in json.load(inf).values():
            for token in kraken_hits.split():
                taxid, count = token.split(":")
                if taxid in hv_taxids:
                    seen_hv_taxids.add(int(taxid))

with open(taxids_out, "w") as outf:
    for taxid in sorted(seen_hv_taxids):
        outf.write("%s\n" % taxid)
