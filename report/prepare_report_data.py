#!/usr/bin/env python3

import sys
import gzip
import os

ROOT_DIR, MGS_PIPELINE_DIR, BIOPROJECT, SAMPLES, REPORT_NAME = sys.argv[1:]
SAMPLE_LIST = SAMPLES.split(',')
REPORT_DIR = os.path.join(ROOT_DIR, f"bioprojects/{BIOPROJECT}/reports/{REPORT_NAME}")

sys.path.insert(0, REPORT_DIR)

def count_reads_by_taxonomic_category(sample):
    cladecounts_dir = os.path.join(REPORT_DIR, "cladecounts")
    cladecounts_file = os.path.join(cladecounts_dir, f"{sample}.tsv.gz")


    output_dir = os.path.join(REPORT_DIR, "tax-categories")
    os.makedirs(output_dir, exist_ok=True)
    output_file_path = os.path.join(output_dir, f"{sample}.tax-category.tsv")

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
    with open(output_file_path, "w") as out:
        out.write("taxid\ttaxname\tparent taxid\tcount\n")
        for taxid, data in tax_data.items():
            out.write(f"{taxid}\t{data['name']}\t{data['parent']}\t{counts[taxid]}\n")

if __name__ == "__main__":
    for sample in SAMPLE_LIST:
        count_reads_by_taxonomic_category(sample)
