#!/usr/bin/env python3

import re
import os
import warnings
import glob
import gzip
import json
import time
import atexit
import argparse
import tempfile
import contextlib
import subprocess
import numpy as np
import random
from collections import Counter
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

S3_BUCKET = None
WORK_ROOT = None
THISDIR = os.path.abspath(os.path.dirname(__file__))

COLOR_RED = "\x1b[0;31m"
COLOR_GREEN = "\x1b[0;32m"
COLOR_CYAN = "\x1b[0;36m"
COLOR_END = "\x1b[0m"


def check_call_shell(cmd):
    # Unlike subprocess.check_call, if any member of the pipeline fails then
    # this fails too.
    subprocess.check_call(["bash", "-c", "set -o  pipefail; %s" % cmd])


def check_output_shell(cmd):
    # Unlike subprocess.check_output, if any member of the pipeline fails then
    # this fails too.
    return subprocess.check_output(
        ["bash", "-c", "set -o  pipefail; %s" % cmd]
    )


def work_fname(*fnames):
    return os.path.join(THISDIR, WORK_ROOT, *fnames)


def get_samples(args):
    if args.sample:
        return [args.sample]

    with open(
        work_fname("bioprojects", args.bioproject, "metadata", "metadata.tsv")
    ) as inf:
        return [line.strip().split("\t")[0] for line in inf]


@contextlib.contextmanager
def tempdir(stage, msg):
    olddir = os.getcwd()
    with tempfile.TemporaryDirectory() as workdir:
        os.chdir(workdir)
        try:
            print("%s: handling %s in %s" % (stage, msg, workdir))
            yield workdir
        finally:
            os.chdir(olddir)


def exists_s3_prefix(s3_path):
    try:
        subprocess.check_call(
            ["aws", "s3", "ls", s3_path], stdout=subprocess.DEVNULL
        )
        return True  # exit code 0 if present
    except subprocess.CalledProcessError as e:
        if e.returncode == 1:
            return False  # exit code 1 if absent

        raise  # any other exit code means something else is wrong


def ls_s3_dir(s3_dir, min_size=0, min_date=""):
    try:
        cmd_out = subprocess.check_output(["aws", "s3", "ls", s3_dir])
    except subprocess.CalledProcessError as e:
        if e.returncode == 1:
            return []  # exit code 1 if absent or empty
        raise  # any other exit code means something is wrong

    for line in cmd_out.split(b"\n"):
        if not line.strip():
            continue
        try:
            date, time, size, fname = line.split()
        except ValueError:
            print(line)
            print(s3_dir)
            raise

        if int(size) < min_size:
            continue
        if date.decode("utf-8") < min_date:
            continue

        yield fname.decode("utf-8")


def get_adapters(in1, in2, adapter1_fname, adapter2_fname):
    output = subprocess.check_output(
        [
            "AdapterRemoval",
            "--file1",
            in1,
            "--file2",
            in2,
            "--identify-adapters",
            "--qualitymax",
            "45",  # Aviti goes up to N
            "--threads",
            "4",
        ]
    )
    output = output.decode("utf-8")

    for line in output.split("\n"):
        if "--adapter1:" in line:
            adapter1 = line.replace("--adapter1:", "").strip()
        elif "--adapter2:" in line:
            adapter2 = line.replace("--adapter2:", "").strip()

    for adapter, fname in [
        [adapter1, adapter1_fname],
        [adapter2, adapter2_fname],
    ]:
        if not all(x in "ACTGN" for x in adapter) or len(adapter) < 20:
            print(output)
            raise Exception(
                "Invalid adapter %r for %r and %r" % (adapter, in1, in2)
            )
        with open(fname, "w") as outf:
            outf.write(adapter)


def adapter_removal(args, dirname, trim_quality, collapse):
    adapter_dir = work_fname("bioprojects", args.bioproject, "adapters")
    try:
        os.mkdir(adapter_dir)
    except FileExistsError:
        pass

    available_inputs = get_files(args, "raw")
    existing_outputs = get_files(args, "cleaned", min_size=100)

    for sample in get_samples(args):
        raw1 = "%s_1.fastq.gz" % sample
        raw2 = "%s_2.fastq.gz" % sample

        if raw1 not in available_inputs or raw2 not in available_inputs:
            print("Skipping %s" % sample)
            continue

        if (
            "%s.collapsed.gz" % sample in existing_outputs
            and "%s.pair1.truncated.gz" % sample in existing_outputs
            and "%s.pair2.truncated.gz" % sample in existing_outputs
        ):
            # Already done
            continue

        with tempdir("adapter_removal", sample) as workdir:
            in1 = "in1.fastq.gz"
            in2 = "in2.fastq.gz"

            subprocess.check_call(
                [
                    "aws",
                    "s3",
                    "cp",
                    "%s/%s/raw/%s" % (S3_BUCKET, args.bioproject, raw1),
                    in1,
                ]
            )
            subprocess.check_call(
                [
                    "aws",
                    "s3",
                    "cp",
                    "%s/%s/raw/%s" % (S3_BUCKET, args.bioproject, raw2),
                    in2,
                ]
            )

            adapter1_fname = os.path.join(adapter_dir, "%s.fwd" % sample)
            adapter2_fname = os.path.join(adapter_dir, "%s.rev" % sample)

            if not os.path.exists(adapter1_fname) or not os.path.exists(
                adapter2_fname
            ):
                get_adapters(in1, in2, adapter1_fname, adapter2_fname)

            with open(adapter1_fname) as inf:
                adapter1 = inf.read().strip()
            with open(adapter2_fname) as inf:
                adapter2 = inf.read().strip()

            cmd = [
                "AdapterRemoval",
                "--file1",
                in1,
                "--file2",
                in2,
                "--basename",
                sample,
                "--threads",
                "4",
                "--qualitymax",
                "45",  # Aviti goes up to N
                "--adapter1",
                adapter1,
                "--adapter2",
                adapter2,
                "--gzip",
            ]

            if trim_quality:
                cmd.extend(["--trimns", "--trimqualities"])
            if collapse:
                cmd.append("--collapse")

            subprocess.check_call(cmd)

            for output in glob.glob("%s.*" % sample):
                subprocess.check_call(
                    [
                        "aws",
                        "s3",
                        "cp",
                        output,
                        "%s/%s/%s/" % (S3_BUCKET, args.bioproject, dirname),
                    ]
                )


def clean(args):
    adapter_removal(args, "cleaned", trim_quality=True, collapse=True)


def rmadapter(args):
    adapter_removal(args, "noadapters", trim_quality=False, collapse=False)


def get_files(args, dirname, min_size=1, min_date=""):
    return set(
        ls_s3_dir(
            "%s/%s/%s/" % (S3_BUCKET, args.bioproject, dirname),
            min_size=min_size,
            min_date=min_date,
        )
    )


def ribofrac(args, subset_size=1000):
    """Fast algorithm to compute fraction of reads identified as rRNA by RiboDetector"""

    available_inputs = get_files(
        args,
        "cleaned",
        # tiny files are empty; ignore them
        min_size=100,
    )
    existing_outputs = get_files(args, "ribofrac", min_date="2023-10-12")

    def first_subset_fastq(file_paths, subset_size):
        """Selects the first subset of reads from gzipped fastq files"""
        print(
            f"Counting reads in input and selecting the first {subset_size}..."
        )
        output_files = []
        total_reads = 0
        # Count the total number of reads only for the first input file
        with gzip.open(file_paths[0], "rt") as f:
            total_reads = sum(1 for _ in FastqGeneralIterator(f))

        for fp in file_paths:
            # When reads in file < subset_size, return actual number of reads in subset
            actual_subset_size = 0

            with gzip.open(fp, "rt") as f:
                # Create an output file handle
                output_file = fp.replace(".fq.gz", ".subset.fq")
                with open(output_file, "w") as out_handle:
                    for index, (title, seq, qual) in enumerate(
                        FastqGeneralIterator(f)
                    ):
                        if index >= subset_size:
                            break
                        actual_subset_size += 1
                        out_handle.write(
                            "@%s\n%s\n+\n%s\n" % (title, seq, qual)
                        )

                output_files.append(output_file)

        return output_files, total_reads, actual_subset_size

    def file_integrity_check(filename):
        """
        Checks if the file exists and contains any reads.
        If either condition fails, it raises a warning.
        """
        if not os.path.exists(filename):
            warnings.warn(f"File {filename} does not exist!")
            return False

        with gzip.open(filename, "rt") as f:
            try:
                # Check for the first read using the iterator
                next(FastqGeneralIterator(f))
            except StopIteration:
                warnings.warn(f"File {filename} contains no reads!")
                return False

        return True

    for sample in get_samples(args):
        # Check for name of output file
        sample_output_file = sample + ".ribofrac.txt"
        if sample_output_file in existing_outputs:
            continue
        # Track for error handling
        total_files_in_sample = 0
        empty_files_in_sample = 0

        total_reads_dict = {}
        subset_reads_dict = {}
        rrna_reads_dict = {}
        file_present = None 
        for potential_input in available_inputs:
            if not potential_input.startswith(sample):
                continue
            if ".settings" in potential_input:
                continue
            if "discarded" in potential_input:
                continue
            file_present = True 
            total_files_in_sample += 1

            # Number of output and input files must match
            tmp_fq_output = potential_input.replace(
                ".gz", ".subset.nonrrna.fq"
            )
            tmp_fq_outputs = [tmp_fq_output]
            inputs = [potential_input]

            if ".pair1." in potential_input:
                tmp_fq_outputs.append(
                    tmp_fq_output.replace(".pair1.", ".pair2.")
                )
                inputs.append(potential_input.replace(".pair1.", ".pair2."))
            elif ".pair2." in potential_input:
                # Ribodetector handles pair1 and pair2 together.
                continue

            with tempdir("ribofrac", sample + " inputs") as workdir:
                for input_fname in inputs:
                    subprocess.check_call(
                        [
                            "aws",
                            "s3",
                            "cp",
                            "%s/%s/cleaned/%s"
                            % (S3_BUCKET, args.bioproject, input_fname),
                            input_fname,
                        ]
                    )

                    # Check file integrity
                    file_valid = file_integrity_check(input_fname)
                    if not file_valid:
                        empty_files_in_sample += 1

                    # Ribodetector gets angry if the .fq extension isn't in the filename
                    os.rename(
                        input_fname, input_fname.replace(".gz", ".fq.gz")
                    )

                if total_files_in_sample == empty_files_in_sample:
                    print(f"Skipping {sample}... all files are empty.")
                    continue

                # Add .fq extensions to input files
                inputs = [i.replace(".gz", ".fq.gz") for i in inputs]

                # Get subset of inputs
                subsets, total_reads, subset_reads = first_subset_fastq(
                    inputs, subset_size
                )
                subset_reads_dict[inputs[0]] = subset_reads
                total_reads_dict[inputs[0]] = total_reads

                # Compute average read lengths. For paired-end reads, average length is
                # computed only from pair1 reads.
                print("Calculating average read length...")

                def calculate_average_read_length(file_path):
                    total_len = 0
                    total_reads = 0
                    with open(file_path, "rt") as inf:
                        for title, sequence, quality in FastqGeneralIterator(
                            inf
                        ):
                            total_len += len(sequence)
                            total_reads += 1
                    return round(total_len / total_reads)

                avg_length = calculate_average_read_length(subsets[0])
                print("Done. Average read length is ", avg_length)

                ribodetector_cmd = [
                    "ribodetector_cpu",
                    "--ensure",
                    "rrna",
                    "--threads",
                    "28",
                ]
                ribodetector_cmd.extend(["--len", str(avg_length)])

                ribodetector_cmd.append("--input")
                ribodetector_cmd.extend(subsets)

                # RiboDetector outputs fastq files containing non-rRNA sequences
                # https://github.com/hzi-bifo/RiboDetector
                ribodetector_cmd.append("--output")
                ribodetector_cmd.extend(tmp_fq_outputs)

                subprocess.check_call(ribodetector_cmd)

                # Count number of rRNA reads in subset
                non_rrna_count = sum(
                    1
                    for _ in FastqGeneralIterator(
                        open(tmp_fq_outputs[0], "rt")
                    )
                )
                rrna_reads_dict[inputs[0]] = subset_reads - non_rrna_count
        if not file_present:
            print("%s wasn't processed by ribofrac because it's not present in %s/%s/cleaned" % (sample, S3_BUCKET, args.bioproject))
            continue


        # Calculate the weighted average fraction of rRNA reads across all inputs in sample using numpy
        # Extract the fractions of rRNA reads for each input
        fractions_rrna_in_subset = [
            rrna_reads_dict[input_filename] / subset_reads_dict[input_filename]
            for input_filename in total_reads_dict
        ]
        # Use the total number of reads for each input as weights
        weights = list(total_reads_dict.values())
        try:
            weighted_rrna_fraction = np.average(
                fractions_rrna_in_subset, weights=weights
            )
        except ZeroDivisionError:
            weighted_rrna_fraction = 0
            print("Sample doesn't contain ribosomal reads.")
        fraction_rrna = round(weighted_rrna_fraction, 4)

        print(
            f"Estimated fraction of rRNA reads in {sample} = {round(fraction_rrna*100, 2)}%"
        )

        # Save fraction of rRNA reads
        with tempdir("ribofrac", sample + "_output") as workdir:
            ribofrac_file = os.path.join(workdir, f"{sample}.ribofrac.txt")

            with open(ribofrac_file, "w") as txt_file:
                txt_file.write(str(fraction_rrna))

            subprocess.check_call(
                [
                    "aws",
                    "s3",
                    "cp",
                    ribofrac_file,
                    "%s/%s/ribofrac/" % (S3_BUCKET, args.bioproject),
                ]
            )


def riboreads(args):
    """Save title of reads identified as rRNA by RiboDetector to AWS"""

    available_inputs = get_files(
        args,
        "cleaned",
        # tiny files are empty; ignore them
        min_size=100,
    )
    existing_outputs = get_files(args, "riboreads", min_date="2023-10-10")

    for sample in get_samples(args):
        # Check for name of output file
        sample_output_file = sample + ".riboreads.txt"
        if sample_output_file in existing_outputs:
            continue

        sample_reads = []
        sample_nonrrna_reads = []
        sample_rrna_reads = []
        for potential_input in available_inputs:
            if not potential_input.startswith(sample):
                continue
            if ".settings" in potential_input:
                continue
            if "discarded" in potential_input:
                continue

            # Number of output and input files must match
            tmp_fq_output = potential_input.replace(".gz", ".nonrrna.fq")
            tmp_fq_outputs = [tmp_fq_output]
            inputs = [potential_input]

            if ".pair1." in potential_input:
                tmp_fq_outputs.append(
                    tmp_fq_output.replace(".pair1.", ".pair2.")
                )
                inputs.append(potential_input.replace(".pair1.", ".pair2."))
            elif ".pair2." in potential_input:
                # Ribodetector handles pair1 and pair2 together.
                continue

            with tempdir("riboreads", sample + " inputs") as workdir:
                for input_fname in inputs:
                    subprocess.check_call(
                        [
                            "aws",
                            "s3",
                            "cp",
                            "%s/%s/cleaned/%s"
                            % (S3_BUCKET, args.bioproject, input_fname),
                            input_fname,
                        ]
                    )

                    # Ribodetector gets angry if the .fq extension isn't in the filename
                    os.rename(
                        input_fname, input_fname.replace(".gz", ".fq.gz")
                    )

                # Add .fq extensions to input files
                inputs = [i.replace(".gz", ".fq.gz") for i in inputs]

                # Compute average read lengths. For paired-end reads, average length is
                # computed only from pair1 reads.
                print("Calculating average read length...")

                def calculate_average_read_length(file_path):
                    total_len = 0
                    total_reads = 0
                    with gzip.open(file_path, "rt") as inf:
                        for title, sequence, quality in FastqGeneralIterator(
                            inf
                        ):
                            total_len += len(sequence)
                            total_reads += 1
                    return round(total_len / total_reads)

                avg_length = calculate_average_read_length(inputs[0])
                print("Done. Average read length is ", avg_length)

                ribodetector_cmd = [
                    "ribodetector_cpu",
                    "--ensure",
                    "rrna",
                    "--threads",
                    "28",
                    "--chunk_size",
                    "512",
                ]
                ribodetector_cmd.extend(["--len", str(avg_length)])

                ribodetector_cmd.append("--input")
                ribodetector_cmd.extend(inputs)

                # RiboDetector outputs fastq files containing non-rRNA sequences
                # https://github.com/hzi-bifo/RiboDetector
                ribodetector_cmd.append("--output")
                ribodetector_cmd.extend(tmp_fq_outputs)

                subprocess.check_call(ribodetector_cmd)

                def parse_reads(file_path):
                    total_reads = 0
                    seq_titles = []

                    # Check if the file is gzipped
                    _, file_extension = os.path.splitext(file_path)

                    if file_extension == ".gz":
                        open_func = gzip.open
                        mode = "rt"
                    else:
                        open_func = open
                        mode = "r"

                    with open_func(file_path, mode) as inf:
                        for title, sequence, quality in FastqGeneralIterator(
                            inf
                        ):
                            seq_titles.append(title)
                    return seq_titles

                # Collect and count all input reads. Paired-end reads are counted once.
                reads = parse_reads(inputs[0])
                sample_reads.extend(reads)

                # Collect and count non-rRNA reads. Paired-end reads are counted once.
                non_rrna_reads = parse_reads(tmp_fq_outputs[0])
                sample_nonrrna_reads.extend(non_rrna_reads)

        # Collect rRNA reads
        sample_rrna_reads = list(set(sample_reads) - set(sample_nonrrna_reads))

        print(
            f"Number of rRNA reads in {sample} = {len(sample_rrna_reads)}/{len(sample_reads)} "
            f"({round(len(sample_rrna_reads)/len(sample_reads)*100)}%)"
        )

        # Save titles of rRNA reads
        # For paired-end reads, only the title of the first reads is saved
        with tempdir("riboreads", sample + "_output") as workdir:
            riboreads_file = os.path.join(workdir, f"{sample}.riboreads.txt")
            gzipped_file_path = riboreads_file + ".gz"

            # Write and gzip the text file
            with gzip.open(gzipped_file_path, "wb") as gzipped_file:
                for title in sample_rrna_reads:
                    gzipped_file.write(title + "\n")

            subprocess.check_call(
                [
                    "aws",
                    "s3",
                    "cp",
                    gzipped_file_path,
                    "%s/%s/riboreads/" % (S3_BUCKET, args.bioproject),
                ]
            )


def interpret(args):
    available_inputs = get_files(
        args,
        "cleaned",
        # tiny files are empty; ignore them
        min_size=100,
    )
    existing_outputs = get_files(args, "processed")

    for sample in get_samples(args):
        for potential_input in available_inputs:
            if not potential_input.startswith(sample):
                continue
            if ".settings" in potential_input:
                continue
            if "discarded" in potential_input:
                continue

            output = potential_input.replace(".gz", ".kraken2.tsv")
            inputs = [potential_input]
            if ".pair1." in output:
                output = output.replace(".pair1.", ".")
                inputs.append(potential_input.replace(".pair1.", ".pair2."))
            elif ".pair2" in output:
                # We handle pair1 and pair2 together.
                continue

            compressed_output = output + ".gz"
            if compressed_output in existing_outputs:
                continue

            with tempdir("interpret", ", ".join(inputs)) as workdir:
                for input_fname in inputs:
                    subprocess.check_call(
                        [
                            "aws",
                            "s3",
                            "cp",
                            "%s/%s/cleaned/%s"
                            % (S3_BUCKET, args.bioproject, input_fname),
                            input_fname,
                        ]
                    )

                kraken_cmd = [
                    "/home/ec2-user/kraken2-install/kraken2",
                    "--db",
                    "/home/ec2-user/kraken-db/",
                    "--use-names",
                    "--threads",
                    "8",
                    "--output",
                    output,
                ]
                if len(inputs) > 1:
                    kraken_cmd.append("--paired")
                kraken_cmd.extend(inputs)

                subprocess.check_call(kraken_cmd)
                subprocess.check_call(["gzip", output])
                subprocess.check_call(
                    [
                        "aws",
                        "s3",
                        "cp",
                        compressed_output,
                        "%s/%s/processed/" % (S3_BUCKET, args.bioproject),
                    ]
                )


def cladecounts(args):
    available_inputs = get_files(args, "processed")
    existing_outputs = get_files(
        args, "cladecounts", min_size=100, min_date="2023-05-19"
    )

    for sample in get_samples(args):
        output = "%s.tsv.gz" % sample
        if output in existing_outputs:
            continue

        if not any(x.startswith(sample) for x in available_inputs):
            continue

        subprocess.check_call(
            ["./count_clades.sh", S3_BUCKET, args.bioproject, sample]
        )


SAMPLE_READS_TARGET_LEN = 100_000


def samplereads(args):
    human_viruses = set()
    with open(os.path.join(THISDIR, "human-viruses.tsv")) as inf:
        for line in inf:
            taxid, _ = line.strip().split("\t")
            human_viruses.add(int(taxid))

    parents = {}
    with open(os.path.join(THISDIR, "dashboard", "nodes.dmp")) as inf:
        for line in inf:
            child_taxid, parent_taxid, *_ = line.replace("\t|\n", "").split(
                "\t|\t"
            )
            parents[int(child_taxid)] = int(parent_taxid)

    def taxid_under(clade, taxid):
        while taxid not in [0, 1]:
            if taxid == clade:
                return True
            taxid = parents[taxid]
        return False

    def taxid_matches(taxid, category):
        if category == "all":
            return True
        if category == "humanviral":
            return taxid in human_viruses
        return taxid_under(
            {
                "bacterial": 2,
                "viral": 10239,
            }[category],
            taxid,
        )

    available_inputs = get_files(args, "processed")
    existing_outputs = get_files(args, "samplereads", min_date="2023-11-03")
    for sample in get_samples(args):
        output = "%s.sr.tsv.gz" % sample
        if output in existing_outputs:
            continue

        inputs = [x for x in available_inputs if x.startswith(sample)]
        if not any(inputs):
            continue

        read_ids = {}
        full_counts = Counter()
        fname_counts = defaultdict(Counter)

        for fname in inputs:
            read_ids[fname] = {
                "all": [],
                "bacterial": [],
                "viral": [],
                "humanviral": [],
            }
            process = subprocess.Popen(
                [
                    "aws",
                    "s3",
                    "cp",
                    "%s/%s/processed/%s" % (S3_BUCKET, args.bioproject, fname),
                    "-",
                ],
                stdout=subprocess.PIPE,
                shell=False,
            )

            try:
                with gzip.open(process.stdout, "rt") as inf:
                    for line in inf:
                        bits = line.rstrip("\n").split("\t")
                        read_id = bits[1]
                        full_assignment = bits[2]

                        taxid = int(full_assignment.split()[-1].rstrip(")"))

                        for category in read_ids[fname]:
                            if taxid_matches(taxid, category):
                                full_counts[category] += 1
                                fname_counts[fname][category] += 1
                                if (
                                    len(read_ids[fname][category])
                                    < SAMPLE_READS_TARGET_LEN
                                ):
                                    read_ids[fname][category].append(read_id)
            finally:
                process.terminate()

        subsetted_ids = {}
        for category, full_count in full_counts.items():
            subsetted_ids[category] = []
            for fname in read_ids:
                if full_count <= SAMPLE_READS_TARGET_LEN:
                    subsetted_ids[category].extend(read_ids[fname][category])
                else:
                    target = (
                        SAMPLE_READS_TARGET_LEN
                        * fname_counts[fname][category]
                        // full_count
                    )

                    subsetted_ids[category].extend(
                        read_ids[fname][category][:target]
                    )

        with tempdir("samplereads", sample) as workdir:
            with gzip.open(output, "wt") as outf:
                for category in sorted(subsetted_ids):
                    for selected_read_id in sorted(subsetted_ids[category]):
                        outf.write(
                            "%s\t%s\n" % (category[0], selected_read_id)
                        )

            subprocess.check_call(
                [
                    "aws",
                    "s3",
                    "cp",
                    output,
                    "%s/%s/samplereads/" % (S3_BUCKET, args.bioproject),
                ]
            )


def readlengths(args):
    available_samplereads_inputs = get_files(args, "samplereads")
    available_cleaned_inputs = get_files(args, "cleaned")
    existing_outputs = get_files(args, "readlengths", min_date="2023-11-04")

    for sample in get_samples(args):
        output = "%s.rl.json.gz" % sample
        if output in existing_outputs:
            continue

        inputs = [
            x for x in available_samplereads_inputs if x.startswith(sample)
        ]
        if not any(inputs):
            continue
        (fname,) = inputs

        target_read_ids = defaultdict(set)
        process = subprocess.Popen(
            [
                "aws",
                "s3",
                "cp",
                "%s/%s/samplereads/%s" % (S3_BUCKET, args.bioproject, fname),
                "-",
            ],
            stdout=subprocess.PIPE,
            shell=False,
        )
        try:
            with gzip.open(process.stdout, "rt") as inf:
                for line in inf:
                    bits = line.rstrip("\n").split("\t")
                    category, read_id = bits
                    target_read_ids[read_id].add(category)
        finally:
            process.terminate()

        inputs = [x for x in available_cleaned_inputs if x.startswith(sample)]
        assert inputs

        lengths = {}
        for category in "abhv":
            lengths[category] = {"NC": 0}

        for fname in inputs:
            if ".collapsed." not in fname:
                # We can only get fragment lengths from cases where we could
                # collapse.  Fragments longer than fwd + rev - minoverlap could be
                # any length for all we know.
                #
                # (It would be possible to do better by aligning to genomes, but
                # that's a ton of work)
                continue

            process = subprocess.Popen(
                [
                    "aws",
                    "s3",
                    "cp",
                    "%s/%s/cleaned/%s" % (S3_BUCKET, args.bioproject, fname),
                    "-",
                ],
                stdout=subprocess.PIPE,
                shell=False,
            )
            try:
                with gzip.open(process.stdout, "rt") as inf:
                    for title, sequence, quality in FastqGeneralIterator(inf):
                        title = title.split()[0]
                        if title not in target_read_ids:
                            continue

                        for category in target_read_ids[title]:
                            seql = len(sequence)
                            if seql not in lengths[category]:
                                lengths[category][seql] = 1
                            else:
                                lengths[category][seql] += 1

                        del target_read_ids[title]
            finally:
                process.terminate()

        # We removed as we went, so any left here are non-collapsed
        for target_read_id, categories in target_read_ids.items():
            for category in categories:
                lengths[category]["NC"] += 1

        with tempdir("readlengths", sample) as workdir:
            with gzip.open(output, "wt") as outf:
                json.dump(lengths, outf)

            subprocess.check_call(
                [
                    "aws",
                    "s3",
                    "cp",
                    output,
                    "%s/%s/readlengths/" % (S3_BUCKET, args.bioproject),
                ]
            )


def humanviruses(args):
    human_viruses = {}
    with open(os.path.join(THISDIR, "human-viruses.tsv")) as inf:
        for line in inf:
            taxid, name = line.strip().split("\t")
            human_viruses[int(taxid)] = name

    available_inputs = get_files(args, "processed")
    existing_outputs = get_files(args, "humanviruses")

    for sample in get_samples(args):
        output = "%s.humanviruses.tsv" % sample
        if output in existing_outputs:
            continue

        inputs = [
            input_fname
            for input_fname in available_inputs
            if input_fname.startswith(sample)
        ]
        if not inputs:
            continue

        counts = Counter()

        for input_fname in inputs:
            with tempdir("humanviruses", sample) as workdir:
                subprocess.check_call(
                    [
                        "aws",
                        "s3",
                        "cp",
                        "%s/%s/processed/%s"
                        % (S3_BUCKET, args.bioproject, input_fname),
                        input_fname,
                    ]
                )

                with gzip.open(input_fname, "rt") as inf:
                    for line in inf:
                        (taxid,) = re.findall("[(]taxid ([0-9]+)[)]", line)
                        taxid = int(taxid)
                        if taxid in human_viruses:
                            counts[taxid] += 1

        with tempdir("humanviruses", sample) as workdir:
            with open(output, "w") as outf:
                for taxid, count in sorted(counts.items()):
                    outf.write(
                        "%s\t%s\t%s\n" % (taxid, count, human_viruses[taxid])
                    )

            subprocess.check_call(
                [
                    "aws",
                    "s3",
                    "cp",
                    output,
                    "%s/%s/humanviruses/%s"
                    % (S3_BUCKET, args.bioproject, output),
                ]
            )


def allmatches(args):
    human_viruses = {}
    with open(os.path.join(THISDIR, "human-viruses.tsv")) as inf:
        for line in inf:
            taxid, name = line.strip().split("\t")
            human_viruses[int(taxid)] = name

    available_inputs = get_files(args, "processed")
    existing_outputs = get_files(args, "allmatches")

    for sample in get_samples(args):
        output = "%s.allmatches.tsv" % sample
        if output in existing_outputs:
            continue

        inputs = [
            input_fname
            for input_fname in available_inputs
            if input_fname.startswith(sample)
        ]
        if not inputs:
            continue

        with tempdir("allmatches", sample) as workdir:
            kept = []
            for input_fname in inputs:
                subprocess.check_call(
                    [
                        "aws",
                        "s3",
                        "cp",
                        "%s/%s/processed/%s"
                        % (S3_BUCKET, args.bioproject, input_fname),
                        input_fname,
                    ]
                )

                with gzip.open(input_fname, "rt") as inf:
                    for line in inf:
                        keep = False
                        try:
                            taxid_matches = line.strip().split("\t")[4]
                            for taxid_match in taxid_matches.split(" "):
                                taxid, n_kmers = taxid_match.split(":")
                                if taxid == "A":
                                    continue  # ambiguous nucleotide
                                if taxid == "|":
                                    continue  # paired end transition
                                taxid = int(taxid)
                                if taxid in human_viruses:
                                    keep = True
                        except Exception:
                            print(line)
                            raise
                        if keep:
                            kept.append(line)

            with open(output, "w") as outf:
                for line in kept:
                    outf.write(line)

            subprocess.check_call(
                [
                    "aws",
                    "s3",
                    "cp",
                    output,
                    "%s/%s/allmatches/%s"
                    % (S3_BUCKET, args.bioproject, output),
                ]
            )


def hvreads(args):
    available_inputs = get_files(args, "allmatches")
    available_cleaned_inputs = get_files(
        args,
        "cleaned",
        # tiny files are empty; ignore them
        min_size=100,
    )

    existing_outputs = get_files(
        args,
        "hvreads",
        # date we added kraken assignments
        min_date="2023-10-24",
    )

    for sample in get_samples(args):
        output = "%s.hvreads.json" % sample
        if output in existing_outputs:
            continue

        input_fname = "%s.allmatches.tsv" % sample
        if input_fname not in available_inputs:
            continue

        all_matches = [
            x.strip().split("\t")
            for x in subprocess.check_output(
                [
                    "aws",
                    "s3",
                    "cp",
                    "%s/%s/allmatches/%s"
                    % (S3_BUCKET, args.bioproject, input_fname),
                    "-",
                ]
            )
            .decode("utf-8")
            .split("\n")
            if x.strip()
        ]

        seqs = {}  # seqid -> kraken_assignment, kraken_hits, fwd, rev
        for _, seq_id, kraken_assignment, _, kraken_details in all_matches:
            assignment_taxid = int(
                re.search(r"\(taxid (\d+)\)", kraken_assignment).group(1)
            )

            seqs[seq_id] = [assignment_taxid, kraken_details]
        for cleaned_input in sorted(available_cleaned_inputs):
            if not cleaned_input.startswith(sample):
                continue
            if ".settings" in cleaned_input:
                continue

            with tempdir("hvreads", cleaned_input) as workdir:
                subprocess.check_call(
                    [
                        "aws",
                        "s3",
                        "cp",
                        "%s/%s/cleaned/%s"
                        % (S3_BUCKET, args.bioproject, cleaned_input),
                        cleaned_input,
                    ]
                )

                with gzip.open(cleaned_input, "rt") as inf:
                    for title, sequence, quality in FastqGeneralIterator(inf):
                        seq_id = title.split()[0]
                        if seq_id in seqs:
                            seqs[seq_id].append([sequence, quality])

        with tempdir("hvreads", output) as workdir:
            with open(output, "w") as outf:
                json.dump(seqs, outf, sort_keys=True)
            subprocess.check_call(
                [
                    "aws",
                    "s3",
                    "cp",
                    output,
                    "%s/%s/hvreads/%s" % (S3_BUCKET, args.bioproject, output),
                ]
            )


def alignments(args):
    available_inputs = get_files(
        args,
        "cleaned",
        # tiny files are empty; ignore them
        min_size=100,
    )

    existing_outputs = get_files(args, "alignments", min_date="2023-11-28")

    with open(
        os.path.join(THISDIR, "bowtie", "genomeid-to-taxid.json")
    ) as inf:
        genomeid_to_taxid = json.load(inf)

    for sample in get_samples(args):
        for db in ["human", "hv"]:
            combined_output_compressed = "%s.%s.alignments.tsv.gz" % (
                sample,
                db,
            )
            if combined_output_compressed in existing_outputs:
                continue

            with tempdir("%s alignments" % db, sample) as workdir:
                tmp_outputs = []
                for potential_input in available_inputs:
                    if not potential_input.startswith(sample):
                        continue
                    if ".settings" in potential_input:
                        continue
                    if "discarded" in potential_input:
                        continue

                    tmp_output = potential_input.replace(
                        ".gz", ".alignments.tsv"
                    )
                    inputs = [potential_input]
                    if ".pair1." in tmp_output:
                        tmp_output = tmp_output.replace(".pair1.", ".")
                        inputs.append(
                            potential_input.replace(".pair1.", ".pair2.")
                        )
                    elif ".pair2" in tmp_output:
                        # We handle pair1 and pair2 together.
                        continue

                    tmp_outputs.append(tmp_output)

                    # TODO(jefftk): do we have a problem when read1 or read2 is
                    # too short?  I remember Bowtie choking on these before.

                    full_inputs = [
                        "%s/%s/cleaned/%s"
                        % (S3_BUCKET, args.bioproject, input_fname)
                        for input_fname in inputs
                    ]

                    subprocess.check_call(
                        [
                            os.path.join(THISDIR, "compute-alignments.sh"),
                            db,
                            tmp_output,
                            *full_inputs,
                        ]
                    )

                with gzip.open(combined_output_compressed, "wt") as outf:
                    for tmp_output in tmp_outputs:
                        with open(tmp_output) as inf:
                            for line in inf:
                                if line.startswith("@"):
                                    continue
                                bits = line.rstrip("\n").split("\t")

                                query_name = bits[0]
                                genomeid = bits[2]
                                ref_start = bits[3]
                                # The start position is 1-indexed, but we use
                                # 0-indexing.
                                ref_start = int(ref_start) - 1
                                cigarstring = bits[5]
                                query_len = len(bits[9])

                                as_val = None
                                for token in bits[11:]:
                                    if token.startswith("AS:i:"):
                                        as_val = token.replace("AS:i:", "")
                                assert as_val
                                as_val = int(as_val)

                                if db == "human":
                                    taxid = 9606
                                    genome_name = genomeid
                                else:
                                    taxid, genome_name = genomeid_to_taxid[
                                        genomeid
                                    ]
                                outf.write(
                                    "%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                                    % (
                                        query_name,
                                        genomeid,
                                        taxid,
                                        cigarstring,
                                        ref_start,
                                        as_val,
                                        query_len,
                                    )
                                )

                subprocess.check_call(
                    [
                        "aws",
                        "s3",
                        "cp",
                        combined_output_compressed,
                        "%s/%s/alignments/%s"
                        % (
                            S3_BUCKET,
                            args.bioproject,
                            combined_output_compressed,
                        ),
                    ]
                )

def phred_to_q(phred_score):
    return ord(phred_score) - ord("!")


def average_quality(phred_counts):
    xs = []
    weights = []

    for phred_score, count in phred_counts.items():
        xs.append(phred_to_q(phred_score))
        weights.append(count)

    return round(np.average(xs, weights=weights))


def print_status(args):
    if args.bioproject:
        bioprojects = [args.bioproject]
    else:
        bioprojects = [
            os.path.basename(x)
            for x in glob.glob(work_fname("bioprojects", "*"))
        ]

    running_processes = subprocess.check_output(["ps", "aux"]).decode("utf-8")

    stages = [
        "raw",
        "cleaned",
        "ribofrac",
        "riboreads",
        "processed",
        "cladecounts",
        "humanviruses",
        "allmatches",
        "hvreads",
        "samplereads",
        "readlengths",
        "alignments",
    ]
    short_stages = [
        "raw",
        "clean",
        "rf",
        "rr",
        "kraken",
        "cc",
        "hv",
        "am",
        "hvr",
        "sr",
        "rl",
        "al",
    ]

    papers_to_projects = defaultdict(list)  # paper -> [project]
    for bioproject in bioprojects:
        metadata_dir = work_fname("bioprojects", bioproject, "metadata")
        if not os.path.exists(metadata_dir):
            continue
        with open(os.path.join(metadata_dir, "name.txt")) as inf:
            name = inf.read().strip()
        papers_to_projects[name].append(bioproject)

    name_width = 21
    print(" " * name_width, end="\t")
    print(*short_stages, sep="\t")
    for paper, bioprojects in sorted(papers_to_projects.items()):
        print(paper)

        for bioproject in bioprojects:
            if bioproject in running_processes:
                color = COLOR_CYAN
            else:
                color = ""

            print(
                color
                + ("  " + bioproject).ljust(name_width)
                + (COLOR_END if color else ""),
                end="",
                flush=True,
            )

            fully_processed_fname = work_fname(
                "bioprojects", bioproject, "fully_processed"
            )
            fully_processed = os.path.exists(fully_processed_fname)
            if fully_processed:
                with open(fully_processed_fname) as inf:
                    n_raw = inf.read().strip()
                print(COLOR_GREEN, end="")
                for stage in stages:
                    print("\t", end="")
                    print(n_raw if stage == "raw" else "-", end="")
                print(COLOR_END)
                continue

            s3_bioproject_dir = "%s/%s" % (S3_BUCKET, bioproject)
            metadata_dir = work_fname("bioprojects", bioproject, "metadata")

            stage_counters = defaultdict(Counter)  # sample -> stage -> count

            with open(os.path.join(metadata_dir, "metadata.tsv")) as inf:
                samples = [x.strip().split("\t")[0] for x in inf]

            prev = None
            for stage in stages:
                print("\t", end="", flush=True)

                seen = set()
                for fname in ls_s3_dir("%s/%s/" % (s3_bioproject_dir, stage)):
                    for sample in samples:
                        if fname.startswith(sample):
                            seen.add(sample)

                missing = prev is not None and len(seen) < prev
                color = COLOR_RED if missing else ""

                print(
                    "%s%s%s" % (color, len(seen), COLOR_END if color else ""),
                    end="",
                    flush=True,
                )
                prev = len(seen)

            print()


STAGES_ORDERED = []
STAGE_FNS = {}
for stage_name, stage_fn in [
    ("clean", clean),
    ("riboreads", riboreads),
    ("ribofrac", ribofrac),
    ("interpret", interpret),
    ("cladecounts", cladecounts),
    ("humanviruses", humanviruses),
    ("allmatches", allmatches),
    ("hvreads", hvreads),
    ("samplereads", samplereads),
    ("readlengths", readlengths),
    ("alignments", alignments),
]:
    STAGES_ORDERED.append(stage_name)
    STAGE_FNS[stage_name] = stage_fn


def start():
    parser = argparse.ArgumentParser(
        description="Run the Metagenomic Sequencing Pipeline"
    )

    parser.add_argument(
        "--restricted",
        action="store_true",
        help="Whether to work on private data",
    )

    parser.add_argument(
        "--bioproject", help="The ID of the bioproject to process"
    )
    parser.add_argument(
        "--sample",
        default="",
        help="The SRA run accession of the sample to process.  Leave blank "
        "for all samples.",
    )

    parser.add_argument(
        "--status",
        action="store_true",
        help="Instead of running anything, just print the status of bioprojects",
    )

    parser.add_argument(
        "--stages",
        default=",".join(STAGES_ORDERED),
        help="Comma-separated list of stages to run.  Allowed stages: %s"
        % (", ".join(repr(x) for x in STAGES_ORDERED)),
    )

    parser.add_argument(
        "--skip-stages",
        default="",
        help="Comma-separated list of stages not to run.",
    )

    args = parser.parse_args()

    global S3_BUCKET
    global WORK_ROOT
    if args.restricted:
        S3_BUCKET = "s3://nao-restricted"
        WORK_ROOT = "../mgs-restricted/"
    else:
        S3_BUCKET = "s3://nao-mgs"
        WORK_ROOT = "./"

    if not args.status and not args.bioproject:
        parser.print_help()
        exit(1)

    if args.status:
        print_status(args)
        return

    if not os.path.isdir(work_fname("bioprojects", args.bioproject)):
        raise Exception(
            "Bioproject %s not found in %s/bioprojects"
            % (args.bioproject, WORK_ROOT)
        )

    selected_stages = args.stages.split(",")
    skipped_stages = args.skip_stages.split(",")
    for stage in selected_stages:
        if stage not in STAGE_FNS:
            raise Exception("Unknown stage %r" % stage)

    if skipped_stages != [""]:
        for stage in skipped_stages:
            if stage not in STAGE_FNS:
                raise Exception("Unknown stage %r" % stage)

    for stage in STAGES_ORDERED:
        if stage in selected_stages and stage not in skipped_stages:
            STAGE_FNS[stage](args)


if __name__ == "__main__":
    start()
