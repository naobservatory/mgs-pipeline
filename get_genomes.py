#!/usr/bin/env python3

import os
import json
import subprocess
from collections import defaultdict

print(
    "Script version that creates reference database based on kraken identified human hits."
)
# Two kinds of taxid:
#   hv_taxid: any human virus kraken assigned a read to
#   detailed_taxid: a hv_taxid, or any of its descendents

hv_taxids = set()
print("Reading observed-human-virus-taxids.txt...")
with open("observed-human-virus-taxids.txt") as inf:
    for line in inf:
        hv_taxids.add(int(line.strip()))

detailed_taxids_fname = "detailed-taxids.txt"
hv_taxid_to_detailed_fname = "hv_taxid_to_detailed.json"
hv_taxid_to_detailed = defaultdict(list)

if not os.path.exists(detailed_taxids_fname):
    fetch = []
    for hv_taxid in hv_taxids:
        fetch.append(hv_taxid)
        hv_taxid_to_detailed[hv_taxid].append(hv_taxid)

        for line in (
            subprocess.check_output(["gimme_taxa.py", str(hv_taxid)])
            .decode("utf-8")
            .split("\n")
        ):
            line = line.strip()
            if line.startswith("parent_taxid") or not line:
                continue
            _, descendent_taxid, descendent_name = line.split("\t")
            descendent_taxid = int(descendent_taxid)
            fetch.append(descendent_taxid)
            hv_taxid_to_detailed[hv_taxid].append(descendent_taxid)
    with open(detailed_taxids_fname, "w") as outf:
        for detailed_taxid in fetch:
            outf.write("%s\n" % detailed_taxid)

    with open(hv_taxid_to_detailed_fname, "w") as outf:
        json.dump(hv_taxid_to_detailed, outf)

metadata_fname = "ncbi-fetch-metadata.txt"

if not os.path.exists(metadata_fname):
    print("ncbi-fetch-metadata.txt not found. Downloading from NCBI...")
    subprocess.check_call(
        [
            "ncbi-genome-download",
            "--section",
            "refseq",
            "--verbose",
            "--taxids",
            detailed_taxids_fname,
            "--formats",
            "fasta",
            "--metadata-table",
            metadata_fname,
            "viral",
        ]
    )
