#!/usr/bin/env python3
import subprocess
import os
import json
from collections import defaultdict
import glob
import gzip

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


def build_detailed_taxids(detailed_taxids_fname):
    hv_taxid_to_detailed_fname = "hv_taxid_to_detailed.json"
    if os.path.exists(detailed_taxids_fname):
        return
    hv_taxid_to_detailed = defaultdict(list)
    fetch = []
    for hv_taxid in human_viruses:
        print("Fetching descendants of", hv_taxid, human_viruses[hv_taxid])
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


def fetch_genomes(detailed_taxids_fname):
    metadata_fname = "ncbi-fetch-metadata.txt"
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


def combine_genomes(combined_genomes_fname):
    if os.path.exists(combined_genomes_fname):
        return
    with open(combined_genomes_fname, "w") as outf:
        for fname in glob.glob("genbank/**/*.fna.gz", recursive=True):
            with gzip.open(fname, "rt") as inf:
                outf.writelines(inf.readlines())


def build_db(bowtie_db_prefix, combined_genomes_fname):
    if os.path.exists(bowtie_db_prefix + ".1.bt2"):
        return

    subprocess.check_call(
        [
            "/home/ec2-user/bowtie2-2.5.2-linux-x86_64/bowtie2-build",
            "-f",
            "--threads",
            "32",
            "--verbose",
            combined_genomes_fname,
            bowtie_db_prefix,
        ]
    )


def bowtie_db():
    cd_bowtie_dir()
    human_viruses = load_human_viruses()
    detailed_taxids_fname = "detailed-taxids.txt"
    build_detailed_taxids(detailed_taxids_fname)
    fetch_genomes(detailed_taxids_fname)
    combined_genomes_fname = "combined_genomes.fna"
    combine_genomes(combined_genomes_fname)
    bowtie_db_prefix = "human-viruses"
    build_db(bowtie_db_prefix, combined_genomes_fname)


if __name__ == "__main__":
    bowtie_db()
