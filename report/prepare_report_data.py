#!/usr/bin/env python3

import sys
import gzip

ROOT_DIR, MGS_PIPELINE_DIR = sys.argv[1:]
REPORT_DIR = ROOT_DIR + "/report/"

sys.path.insert(0, REPORT_DIR)

def count_reads_by_taxonomic_category(cladecounts_file):

    # Taxids and names
    tax_data = { 
        0: {"name": "unmatched", "parent": 0},
        1: {"name": "root", "parent": 1},
        10239: {"name": "Viruses", "parent": 1},
        131567: {"name": "cellular organisms", "parent": 1},
        2787854: {"name": "other entries", "parent": 1},
        2787823: {"name": "unclassified entries", "parent": 1},
        2: {"name": "Bacteria", "parent": 131567},
        2759: {"name": "Eukaryota", "parent": 131567},
        2157: {"name": "Archaea", "parent": 131567}
    }
    
    counts = {i: 0 for i in tax_data.keys()}
    
    # Read the clade file and count reads by search taxid
    with gzip.open(cladecounts_file, 'rt') as f:
        for line in f:
            columns = line.strip().split("\t")
            taxid = int(columns[0])
            clade_assignments = int(columns[3])
            if taxid in tax_data:
                counts[taxid] += clade_assignments

    # Write the results to a TSV file
    with open("tax_category_counts.tsv", "w") as out:
        out.write("taxid\ttaxname\tparent taxid\tcount\n")
        for taxid, data in tax_data.items():
            out.write(f"{taxid}\t{data['name']}\t{data['parent']}\t{counts[taxid]}\n")

if __name__ == "__main__":
    count_reads_by_taxonomic_category("tax_category_counts.tsv")  # You should provide the path to your clade file here

