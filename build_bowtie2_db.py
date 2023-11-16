#!/usr/bin/env python3
import subprocess
import os
import json
from collections import defaultdict

THISDIR = os.path.abspath(os.path.dirname(__file__))


def bowtie_db():
    bowtie_dir = os.path.join(THISDIR, "bowtie")
    if not os.path.exists(bowtie_dir):
        os.mkdir(bowtie_dir)
    os.chdir(bowtie_dir)
    human_viruses = {}
    with open(os.path.join(THISDIR, "first-then-human-viruses.tsv")) as inf:
        for line in inf:
            taxid, name = line.strip().split("\t")
            human_viruses[int(taxid)] = name

    detailed_taxids_fname = "detailed-taxids.txt"
    hv_taxid_to_detailed_fname = "hv_taxid_to_detailed.json"
    if not os.path.exists(detailed_taxids_fname):
        hv_taxid_to_detailed = defaultdict(list)
        fetch = []
        for hv_taxid in human_viruses:
            print(hv_taxid, human_viruses[hv_taxid])
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
    try:
        if not os.path.exists(metadata_fname):
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
    except subprocess.CalledProcessError as e:
        print("The subprocess failed with exit code {}".format(e.returncode))
        print(
            "Maybe ncbi-genome-download is not installed. If so, please download it from: \nhttps://github.com/kblin/ncbi-genome-download"
        )
        raise
    combined_genomes_fname = "combined_genomes.fna"
    if not os.path.exists(combined_genomes_fname):
        with open(combined_genomes_fname, "w") as outf:
            for fname in glob.glob("genbank/**/*.fna.gz", recursive=True):
                with gzip.open(fname) as inf:
                    outf.writelines(inf.readlines())

    bowtie_db_prefix = "human-viruses"
    if not os.path.exists(bowtie_db_prefix + ".1.bt2"):
        subprocess.check.call(
            [
                "~/bowtie2-2.5.2-linux-x86_64/bowtie2-build",
                "-f",
                "--threads",
                "32",
                "verbose",
                combined_genomes_fname,
                bowtie_db_prefix,
            ]
        )


if __name__ == "__main__":
    bowtie_db()
