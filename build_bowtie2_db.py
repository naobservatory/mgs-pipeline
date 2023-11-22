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


def fetch_genomes(detailed_taxids_fname, metadata_fname):
    if os.path.exists(metadata_fname):
        return
    print("Fetching viral GenBank genomes...")

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

def create_genome_taxid_map(metadata_fname):
    genome_taxid_map_fname = "genomeid-to-taxid.json"
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

def combine_genomes(combined_genomes_fname):
    if os.path.exists(combined_genomes_fname):
        return
    print("Combining genomes into one large FASTA file...")

    with open(combined_genomes_fname, "w") as outf:
        for fname in glob.glob("genbank/**/*.fna.gz", recursive=True):
            with gzip.open(fname, "rt") as inf:
                outf.writelines(inf.readlines())

def mask_low_complexity_sequences(
        combined_genomes_fname, masked_genomes_fname):
    if os.path.exists(masked_genomes_fname):
        return
    print("Masking low complexity sequences...")

    # If some input genome has a sequence like GGGGGG...GGGGGG then anytime we
    # evaluate a read with that we'll decide it came from the input genome.
    # But that's not actually helpful: low-complexity sequences are a mixture
    # of artifacts (two-color sequencers emit a lot of G when they get no
    # signal) and non-informative sections of genomes.
    #
    # Kraken handles this when building its DB with mask_low_complexity.sh
    # which uses NCBI dustmasker to mask regions of the input genomes that the
    # DUST algorithm identifies as low-complexity.  Do that here too.

    subprocess.check_call([
        "dustmasker",
        "-in", combined_genomes_fname,
        "-out", masked_genomes_fname,
        "-outfmt", "fasta"])

    # Dustmasker lowercases bases to mask them, but Bowtie needs them to be an
    # unknown character.  It doesn't matter which one, so copy Kraken and use
    # x.
    #
    # This regexp replaces all lowercase letters that aren't on lines beginning
    # with '>', which in FASTA means everywhere except in the sequence IDs.
    subprocess.check_call([
        "sed",
        "/^>/!s/[a-z]/x/g",
        "-i",
        masked_genomes_fname])

def build_db(bowtie_db_prefix, genomes_fname):
    if os.path.exists(bowtie_db_prefix + ".1.bt2"):
        return
    print("Building BowtieDB...")

    subprocess.check_call(
        [
            "/home/ec2-user/bowtie2-2.5.2-linux-x86_64/bowtie2-build",
            "-f",
            "--threads",
            "32",
            "--verbose",
            genomes_fname,
            bowtie_db_prefix,
        ]
    )


def bowtie_db():
    cd_bowtie_dir()
    human_viruses = load_human_viruses()
    detailed_taxids_fname = "detailed-taxids.txt"
    build_detailed_taxids(detailed_taxids_fname, human_viruses)
    metadata_fname = "ncbi-fetch-metadata.txt"
    fetch_genomes(detailed_taxids_fname, metadata_fname)
    create_genome_taxid_map(metadata_fname)
    combined_genomes_fname = "combined_genomes.fna"
    combine_genomes(combined_genomes_fname)
    masked_genomes_fname = "masked_genomes.fna"
    mask_low_complexity_sequences(combined_genomes_fname, masked_genomes_fname)
    bowtie_db_prefix = "human-viruses"
    build_db(bowtie_db_prefix, masked_genomes_fname)


if __name__ == "__main__":
    bowtie_db()
