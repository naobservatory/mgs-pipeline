#!/usr/bin/env python3

import csv
import os
import json
import subprocess
from collections import defaultdict

# Two kinds of taxid:
#   hv_taxid: any human virus kraken assigned a read to
#   detailed_taxid: a hv_taxid, or any of its descendents

hv_taxids = set()
with open("human-viruses.tsv") as inf:
    for line in inf:
        taxid = line.split("\t")[0]
        hv_taxids.add(int(taxid))


detailed_taxids_fname = "detailed-taxids.txt"
hv_taxid_to_detailed_fname = "hv_taxid_to_detailed.json"
hv_taxid_to_detailed = defaultdict(list)
if not os.path.exists(detailed_taxids_fname):
    print("Fetching detailed taxids...")
    fetch = []
    for hv_taxid in hv_taxids:
        print("Fetching descendant taxids for hv_taxid %s" % hv_taxid)
        fetch.append(hv_taxid)
        hv_taxid_to_detailed[hv_taxid].append(hv_taxid)
        # check if the script gimme_taxa.py is in the same directory:
        if not os.path.exists("gimme_taxa.py"):
            raise Exception(
                "gimme_taxa.py is required for execution, please download it from: \nhttps://github.com/kblin/ncbi-genome-download/blob/master/contrib/gimme_taxa.py"
            )

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
try:
    if not os.path.exists(metadata_fname):
        print("Fetching viral refseq genomes...")
        subprocess.check_call(
            [
                "ncbi-genome-download",
                "--section",
                "genbank",
                "--taxids",
                detailed_taxids_fname,
                "--formats",
                "fasta",
                "--metadata-table",
                metadata_fname,
                "viral",
            ]
        )
except subprocess.CalledProcessError as e:
    print("The subprocess failed with exit code {}".format(e.returncode))
    print(
        "Maybe ncbi-genome-download is not installed. If so, please download it from:\nhttps://github.com/kblin/ncbi-genome-download"
    )
    raise
