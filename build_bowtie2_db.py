#!/usr/bin/env python3
import subprocess
import os
import json
from collections import defaultdict
import glob
import gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser

THISDIR = os.path.abspath(os.path.dirname(__file__))


def cd_bowtie_dir():
    bowtie_dir = os.path.join(THISDIR, "bowtie")
    if not os.path.exists(bowtie_dir):
        os.mkdir(bowtie_dir)
    os.chdir(bowtie_dir)

def load_human_viruses():
    human_viruses = {}
    with open(os.path.join(THISDIR, "human-viruses.tsv")) as inf:
        for line in inf:
            taxid, name = line.strip().split("\t")
            human_viruses[int(taxid)] = name
    return human_viruses


def build_detailed_taxids(detailed_taxids_fname, human_viruses):
    hv_taxid_to_detailed_fname = "hv_taxid_to_detailed.json"
    if os.path.exists(detailed_taxids_fname):
        return
    print("Downloading all child taxids for each human virus...")

    hv_taxid_to_detailed = defaultdict(list)
    fetch = set()
    for hv_taxid in human_viruses:
        print("Fetching descendants of", hv_taxid, human_viruses[hv_taxid])
        fetch.add(hv_taxid)
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
            fetch.add(descendent_taxid)
            hv_taxid_to_detailed[hv_taxid].append(descendent_taxid)
    with open(detailed_taxids_fname, "w") as outf:
        for detailed_taxid in sorted(fetch):
            outf.write("%s\n" % detailed_taxid)

    with open(hv_taxid_to_detailed_fname, "w") as outf:
        json.dump(hv_taxid_to_detailed, outf)


def fetch_genomes(detailed_taxids_fname, metadata_fname, section, group):
    if os.path.exists(metadata_fname):
        return
    print("Fetching %s %s genomes..." % (section, group))

    subprocess.check_call(
        [
            "ncbi-genome-download",
            "--section",
            section,
            "--taxids",
            detailed_taxids_fname,
            "--formats",
            "fasta",
            "--metadata-table",
            metadata_fname,
            group,
        ]
    )


def create_genome_taxid_map(metadata_fname, genome_taxid_map_fname):
    if os.path.exists(genome_taxid_map_fname):
        return
    print("Creating map from genome to taxid...")

    genome_to_taxid = {}
    with open(metadata_fname) as inf:
        cols = None
        for line in inf:
            bits = line.rstrip("\n").split("\t")
            if not cols:
                cols = bits
                continue

            with gzip.open(bits[cols.index("local_filename")], "rt") as inf:
                for title, sequence in SimpleFastaParser(inf):
                    genome_id, name = title.split(" ", 1)
                    taxid = int(bits[cols.index("taxid")])
                    genome_to_taxid[genome_id] = taxid, name

    with open(genome_taxid_map_fname, "w") as outf:
        json.dump(genome_to_taxid, outf)

def combine_and_mask_genomes(prefix, masked_genomes_s3_path):
    subprocess.check_call([
        os.path.join(THISDIR, "combine_and_mask_genomes.sh"),
        prefix,
        masked_genomes_s3_path,
    ])

def build_db(bowtie_db_prefix, genomes_s3_path):
    if os.path.exists(bowtie_db_prefix + ".1.bt2"):
        return
    print("Building BowtieDB...")

    subprocess.check_call([
        os.path.join(THISDIR, "build_bowtie2_db.sh"),
        genomes_s3_path,
        bowtie_db_prefix,
    ])

def bowtie_db():
    cd_bowtie_dir()

    # Human Viruses
    human_viruses = load_human_viruses()
    detailed_taxids_fname = "detailed-taxids.txt"
    build_detailed_taxids(detailed_taxids_fname, human_viruses)
    metadata_fname = "ncbi-fetch-metadata.txt"
    fetch_genomes(detailed_taxids_fname, metadata_fname, "genbank", "viral")
    genome_taxid_map_fname = "genomeid-to-taxid.json"
    create_genome_taxid_map(metadata_fname, genome_taxid_map_fname)
    masked_genomes_fname = "masked_genomes.fna.gz"
    masked_genomes_s3_path = "s3://nao-mgs/db/%s" % masked_genomes_fname
    combine_and_mask_genomes("genbank/viral", masked_genomes_s3_path)
    bowtie_db_prefix = "human-viruses"
    build_db(bowtie_db_prefix, masked_genomes_s3_path)

    # Human Respiratory Bacteria
    metadata_fname = "ncbi-fetch-metadata-bacterial.txt"
    # This is already a fully-detailed expansion of the target taxids, so we
    # don't need gimme-taxa.
    fetch_genomes("../hb-taxids.txt", metadata_fname, "refseq", "bacteria")
    genome_taxid_map_fname = "genomeid-to-taxid-bacterial.json"
    create_genome_taxid_map(metadata_fname, genome_taxid_map_fname)
    masked_genomes_fname = "masked_genomes_bacterial.fna.gz"
    masked_genomes_s3_path = "s3://nao-mgs/db/%s" % masked_genomes_fname
    combine_and_mask_genomes("refseq/bacteria", masked_genomes_s3_path)
    bowtie_db_prefix = "human-respiratory-bacteria"
    build_db(bowtie_db_prefix, masked_genomes_s3_path)

if __name__ == "__main__":
    bowtie_db()
