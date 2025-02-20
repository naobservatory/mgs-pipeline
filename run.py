#!/usr/bin/env python3

import re
import os
import warnings
import glob
import gzip
import json
import math
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

# Intended for https://github.com/naobservatory/mgs-pipeline/issues/65
# migration.
with open(os.path.join(THISDIR, "reference-suffix.txt")) as inf:
    REFERENCE_SUFFIX = inf.read().strip()

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

def get_sample_priority(sample):
    m = re.findall(r"L00\d$", sample)
    if not m:
        return "A"
    m, = m
    return m

def get_samples(args):
    if args.sample:
        return [args.sample]

    with open(
        work_fname("deliveries", args.delivery, "metadata", "metadata.tsv")
    ) as inf:
        samples = [line.strip().split("\t")[0] for line in inf]

    for sample in samples:
        if "." in sample:
            raise Exception("Bad sample for %s: %s" % (
                args.delivery, sample))

    decorated_samples = [
        (get_sample_priority(sample), sample)
        for sample in samples]
    decorated_samples.sort()

    return [sample for _, sample in decorated_samples]

@contextlib.contextmanager
def tempdir(stage, msg):
    olddir = os.getcwd()
    with tempfile.TemporaryDirectory(
            dir=os.path.expanduser("~/tmp/")) as workdir:
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
        if line.endswith(b"0 "):
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


def get_adapters(ins, adapter1_fname, adapter2_fname):
    cmd = [
        "AdapterRemoval",
    ]
    if len(ins) == 2:
        cmd.extend(
            ["--file1",
             ins[0],
             "--file2",
             ins[1],
             ])
    elif len(ins) == 1:
        cmd.extend(
            ["--interleaved-input",
             "--file1",
             ins[0]])
    else:
        assert False

    cmd.extend([
        "--identify-adapters",
        "--qualitymax",
        "45",  # Aviti goes up to N
        "--threads",
        "4",
    ])


    output = subprocess.check_output(cmd)
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
                "Invalid adapter %r for %r" % (adapter, ins)
            )
        with open(fname, "w") as outf:
            outf.write(adapter)

def is_nanopore(args):
    return args.delivery.startswith("NAO-ONT-")

def rm_human(args):
    return "-Zephyr" in args.delivery

def no_adapters_dirname(args):
    if is_nanopore(args):
        return "raw"
    return "cleaned"

def final_fastq_dirname(args):
    if rm_human(args):
        return "nonhuman"
    else:
        return no_adapters_dirname(args)

def adapter_removal(args, dirname, trim_quality, collapse):
    adapter_dir = work_fname("deliveries", args.delivery, "adapters")
    try:
        os.mkdir(adapter_dir)
    except FileExistsError:
        pass

    available_inputs = get_files(args, "raw")
    existing_outputs = get_files(args, "cleaned", min_size=100)

    for sample in get_samples(args):
        raw1 = "%s_1.fastq.gz" % sample
        raw2 = "%s_2.fastq.gz" % sample
        rawI = "%s.fastq.gz" % sample

        ins_raw = []
        if rawI in available_inputs:
            ins_raw.append(rawI)
        elif raw1 in available_inputs and raw2 in available_inputs:
            ins_raw.append(raw1)
            ins_raw.append(raw2)
        else:
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
            if len(ins_raw) == 2:
                ins = ["in1.fastq.gz",
                       "in2.fastq.gz"]
            elif len(ins_raw) == 1:
                ins = ["in.fastq.gz"]

            for remote_fname, local_fname in zip(
                    ins_raw, ins):
                s3_copy_down(args, "raw", remote_fname, local_fname=local_fname)

            adapter1_fname = os.path.join(adapter_dir, "%s.fwd" % sample)
            adapter2_fname = os.path.join(adapter_dir, "%s.rev" % sample)

            if not os.path.exists(adapter1_fname) or not os.path.exists(
                adapter2_fname
            ):
                get_adapters(ins, adapter1_fname, adapter2_fname)

            with open(adapter1_fname) as inf:
                adapter1 = inf.read().strip()
            with open(adapter2_fname) as inf:
                adapter2 = inf.read().strip()

            cmd = [
                "AdapterRemoval",
            ]

            if len(ins) == 2:
                cmd.extend([
                    "--file1",
                    ins[0],
                    "--file2",
                    ins[1],
                ])
            elif len(ins) == 1:
                cmd.extend([
                    "--interleaved-input",
                    "--file1",
                    ins[0]
                ])
            else:
                assert False

            cmd.extend([
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
            ])

            if trim_quality:
                cmd.extend(["--trimns", "--trimqualities"])
            if collapse:
                cmd.append("--collapse")

            subprocess.check_call(cmd)

            for output in glob.glob("%s.*" % sample):
                s3_copy_up(args, output, dirname)

def clean(args):
    if is_nanopore(args):
        return

    adapter_removal(args, "cleaned", trim_quality=True, collapse=True)

def full_s3_dirname(dirname):
    if dirname in ["raw", "cleaned", "ribofrac", "nonhuman"]:
        return dirname
    return "%s-%s" % (dirname, REFERENCE_SUFFIX)

def s3_dir(args, dirname):
    return "%s/%s/%s/" % (S3_BUCKET, args.delivery, full_s3_dirname(dirname))

def s3_file(args, dirname, fname):
    return "%s%s" % (s3_dir(args, dirname), fname)

def s3_copy_down(args, dirname, remote_fname, local_fname=None):
    if not local_fname:
        local_fname = os.path.basename(remote_fname)
    subprocess.check_call([
        "aws",
        "s3",
        "cp",
        s3_file(args, dirname, remote_fname),
        local_fname,
    ])

def s3_copy_up(args, local_fname, dirname, remote_fname=None):
    if not remote_fname:
        remote_fname = os.path.basename(local_fname)

    subprocess.check_call([
        "aws",
        "s3",
        "cp",
        local_fname,
        s3_file(args, dirname, remote_fname),
    ])


def get_files(args, dirname, min_size=1, min_date=""):
    return set(ls_s3_dir(s3_dir(args, dirname),
                         min_size=min_size,
                         min_date=min_date))



def ribofrac(args, subset_size=1000):
    """Fast algorithm to compute fraction of reads identified as rRNA by RiboDetector"""

    available_inputs = get_files(
        args,
        no_adapters_dirname(args),
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
        for potential_input in available_inputs:
            if not potential_input.startswith(sample):
                continue
            if ".settings" in potential_input:
                continue
            if "discarded" in potential_input:
                continue
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
                    s3_copy_down(args, no_adapters_dirname(args), input_fname)

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
        if not total_files_in_sample:
            print("%s wasn't processed by ribofrac because it's not present in cleaned" % sample)
            continue


        # Calculate the weighted average fraction of rRNA reads across all inputs in sample using numpy
        # Extract the fractions of rRNA reads for each input
        fractions_rrna_in_subset = [
            rrna_reads_dict[input_filename] / subset_reads_dict[input_filename]
            for input_filename in total_reads_dict
        ]

        # Use the total number of reads for each input as weights
        weights = list(total_reads_dict.values())
        weighted_rrna_fraction = np.average(
            fractions_rrna_in_subset, weights=weights
        )
        fraction_rrna = round(weighted_rrna_fraction, 4)

        print(
            f"Estimated fraction of rRNA reads in {sample} = {round(fraction_rrna*100, 2)}%"
        )

        # Save fraction of rRNA reads
        with tempdir("ribofrac", sample + "_output") as workdir:
            ribofrac_file = os.path.join(workdir, f"{sample}.ribofrac.txt")

            with open(ribofrac_file, "w") as txt_file:
                txt_file.write(str(fraction_rrna))

            s3_copy_up(args, ribofrac_file, "ribofrac")

def interpret(args):
    available_inputs = get_files(
        args,
        final_fastq_dirname(args),
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
                    s3_copy_down(args, final_fastq_dirname(args), input_fname)

                kraken_cmd = [
                    "/home/ec2-user/kraken2-install/kraken2",
                    "--use-names",
                    "--output",
                    output,
                ]

                db = "/dev/shm/kraken-db/"
                kraken_cmd.append("--memory-mapping")
                threads = "4"

                assert os.path.exists(db)
                kraken_cmd.append("--db")
                kraken_cmd.append(db)
                kraken_cmd.append("--threads")
                kraken_cmd.append(threads)

                if len(inputs) > 1:
                    kraken_cmd.append("--paired")
                kraken_cmd.extend(inputs)

                subprocess.check_call(kraken_cmd)
                subprocess.check_call(["gzip", output])
                s3_copy_up(args, compressed_output, "processed")

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
            ["./count_clades.sh",
             S3_BUCKET,
             args.delivery,
             sample,
             REFERENCE_SUFFIX]
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
                    s3_file(args, "processed", fname),
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

            s3_copy_up(args, output, "samplereads")

def readlengths(args):
    available_samplereads_inputs = get_files(args, "samplereads")
    available_cleaned_inputs = get_files(args, final_fastq_dirname(args))
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
                s3_file(args, "samplereads", fname),
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
            if ".collapsed." not in fname and not is_nanopore(args):
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
                    s3_file(args, final_fastq_dirname(args), fname),
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

            s3_copy_up(args, output, "readlengths")


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
                s3_copy_down(args, "processed", input_fname)

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

            s3_copy_up(args, output, "humanviruses")

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
                s3_copy_down(args, "processed", input_fname)

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

            s3_copy_up(args, output, "allmatches")

def valreads(args):
    # The subset of hvreads where that pass an alignment threshold.
    available_hvreads_inputs = get_files(args, "hvreads")
    available_alignments2_inputs = get_files(args, "alignments2")
    existing_outputs = get_files(args, "valreads")

    for sample in get_samples(args):
        output = "%s.valreads.json" % sample
        if output in existing_outputs:
            continue

        input_hvreads_fname = "%s.hvreads.json" % sample
        input_alignments2_fname = "%s.hv.alignments2.tsv.gz" % sample

        if input_hvreads_fname not in available_hvreads_inputs:
            continue
        if input_alignments2_fname not in available_alignments2_inputs:
            continue

        with tempdir("valreads", sample) as workdir:
            s3_copy_down(args, "hvreads", input_hvreads_fname)
            s3_copy_down(args, "alignments2", input_alignments2_fname)

            accepted_read_ids = set()
            with gzip.open(input_alignments2_fname, "rt") as inf:
                for line in inf:
                    (query_name, genomeid, taxid, cigarstring, ref_start,
                     as_val, query_len) = line.rstrip("\n").split("\t")

                    length_adjusted_score = int(as_val) / math.log(int(query_len))
                    if length_adjusted_score > 20:
                        accepted_read_ids.add(query_name)

            valreads_out = {}
            with open(input_hvreads_fname) as inf:
               for read_id, record in json.load(inf).items():
                   if read_id in accepted_read_ids:
                       valreads_out[read_id] = record

            with open(output, "w") as outf:
                json.dump(valreads_out, outf)

            s3_copy_up(args, output, "valreads")

def tmpvalreads(args):
    available_hvreads_inputs = get_files(args, "hvreads")
    available_alignments2_inputs = get_files(args, "alignments2")
    existing_outputs = get_files(args, "tmpvalreads")

    for sample in get_samples(args):
        output = "%s.tmpvalreads.json" % sample
        if output in existing_outputs:
            continue

        input_hvreads_fname = "%s.hvreads.json" % sample
        input_alignments2_fname = "%s.hv.alignments2.tsv.gz" % sample

        if input_hvreads_fname not in available_hvreads_inputs:
            continue
        if input_alignments2_fname not in available_alignments2_inputs:
            continue

        with tempdir("tmpvalreads", sample) as workdir:
            s3_copy_down(args, "hvreads", input_hvreads_fname)
            s3_copy_down(args, "alignments2", input_alignments2_fname)

            read_scores = defaultdict(float)
            with gzip.open(input_alignments2_fname, "rt") as inf:
                for line in inf:
                    (query_name, genomeid, taxid, cigarstring, ref_start,
                     as_val, query_len) = line.rstrip("\n").split("\t")

                    length_adjusted_score = int(as_val) / math.log(int(query_len))
                    read_scores[query_name] = max(
                        read_scores[query_name], length_adjusted_score)

            tmpvalreads_out = {}
            with open(input_hvreads_fname) as inf:
               for read_id, record in json.load(inf).items():
                   record.insert(0, read_scores[read_id])
                   tmpvalreads_out[read_id] = record

            with open(output, "w") as outf:
                json.dump(tmpvalreads_out, outf)

            s3_copy_up(args, output, "tmpvalreads")

def hvreads(args):
    available_inputs = get_files(args, "allmatches")
    available_cleaned_inputs = get_files(
        args,
        final_fastq_dirname(args),
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
                    s3_file(args, "allmatches", input_fname),
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
                s3_copy_down(args, final_fastq_dirname(args), cleaned_input)
                with gzip.open(cleaned_input, "rt") as inf:
                    for title, sequence, quality in FastqGeneralIterator(inf):
                        seq_id = title.split()[0]
                        if seq_id.endswith("/1") or seq_id.endswith("/2"):
                            seq_id = seq_id[:-2]
                        if seq_id in seqs:
                            seqs[seq_id].append([sequence, quality])

        with tempdir("hvreads", output) as workdir:
            with open(output, "w") as outf:
                json.dump(seqs, outf, sort_keys=True)
            s3_copy_up(args, output, "hvreads")

DB_DIR="/dev/shm/bowtie-db"
def nonhuman(args):
    if not rm_human(args):
        return

    available_inputs = get_files(
        args,
        no_adapters_dirname(args),
        # tiny files are empty; ignore them
        min_size=100,
    )

    existing_outputs = get_files(args, "nonhuman", min_size=100)

    for sample in get_samples(args):
        output= "%s.fastq.gz" % sample
        if output in existing_outputs:
            continue

        with tempdir("nonhuman", sample) as workdir:
            for potential_input in available_inputs:
                if not potential_input.startswith(sample):
                    continue
                s3_copy_down(args, no_adapters_dirname(args), potential_input)

                local_output="nonhuman.fastq.gz"
                subprocess.check_call([
                    "/home/ec2-user/bowtie2-2.5.2-linux-x86_64/bowtie2",
                    "-x", "%s/chm13.draft_v1.0_plusY" % DB_DIR,
                    "--threads", "4", "--mm",
                    "-U", potential_input,
                    "--un-gz", local_output,
                    "-S", "/dev/null",

                    # Tweak the settings because we're running Nanopore

                    # allow more mismatches
                    "--score-min", "L,0,-0.6",

                    # shorter seed length
                    "-L", "15",

                    # more frequent reseeding
                    "-i", "S,1,0.5",

                    # allow one mismatch in the seed alignment
                    "-N", "1",

                    # less stringent gap penalties
                    "--rdg", "5,3",  # read gap and open
                    "--rfg", "5,3",  # reference gap and open

                ])

                s3_copy_up(args, local_output, "nonhuman", remote_fname=output)

def alignments2(args):
    available_inputs = get_files(
        args,
        "hvreads",
        # tiny files are empty; ignore them
        min_size=100,
    )

    existing_outputs = get_files(args, "alignments2", min_size=100)

    with open(
        "/dev/shm/bowtie-db/v1-pipeline-bowtie-genomeid-to-taxid.json"
    ) as inf:
        genomeid_to_taxid = json.load(inf)

    for sample in get_samples(args):
        combined_output_compressed = "%s.hv.alignments2.tsv.gz" % sample
        if combined_output_compressed in existing_outputs:
            continue

        with tempdir("alignments2", sample) as workdir:
            tmp_outputs = []
            any_output = False
            for potential_input in available_inputs:
                if not potential_input.startswith(sample):
                    continue
                any_output = True

                tmp_output = potential_input.replace(
                    ".hvreads.json", ".alignments2.tsv"
                )

                any_paired = False
                any_collapsed = False

                s3_copy_down(args, "hvreads", potential_input)

                with open(potential_input) as inf, \
                     open("pair1.fastq", "w") as outf1, \
                     open("pair2.fastq", "w") as outf2, \
                     open("collapsed.fastq", "w") as outfC:

                    for title, record in json.load(inf).items():
                        taxid, kraken_info, *reads = record

                        reads = [
                            (s, q)
                            for (s, q) in reads
                            if len(s) >= 20
                        ]

                        if len(reads) == 1:
                            (s, q), = reads
                            outfC.write("@%s\n%s\n+\n%s\n" % (title, s, q))
                            any_collapsed = True
                        elif len(reads) == 2:
                            (s1, q1), (s2, q2) = reads
                            outf1.write("@%s/1\n%s\n+\n%s\n" % (
                                title, s1, q1))
                            outf2.write("@%s/2\n%s\n+\n%s\n" % (
                                title, s2, q2))
                            any_paired = True
                        else:
                            # Both reads were too short after trimming to be
                            # used.
                            continue

                if not any_paired and not any_collapsed:
                    continue

                cmd = [
                    "/home/ec2-user/bowtie2-2.5.2-linux-x86_64/bowtie2"
                ]
                cmd.extend(["--threads", "4", "--mm"])

                cmd.extend(["--no-unal",
                            "--no-sq",
                            "-S", tmp_output])

                # Custom-built HV DB
                cmd.extend(
                    ["-x", "%s/human-viruses" % DB_DIR])
                # When identifying HV reads use looser settings and
                # filter more later.
                cmd.extend(
                    ["--local", "--very-sensitive-local",
                     "--score-min", "G,1,0",
                     "--mp", "4,1"])

                if any_paired:
                    cmd.extend([
                        "-1", "pair1.fastq",
                        "-2", "pair2.fastq",
                    ])
                if any_collapsed:
                    cmd.extend(["-U", "collapsed.fastq"])

                subprocess.check_call(cmd)

                tmp_outputs.append(tmp_output)

            if not any_output:
                continue

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

            s3_copy_up(args, combined_output_compressed, "alignments2")

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
    if args.delivery:
        deliveries = [args.delivery]
    else:
        deliveries = [
            os.path.basename(os.path.dirname(x))
            for x in glob.glob(work_fname("deliveries", "*/"))
        ]

    running_processes = subprocess.check_output(["ps", "aux"]).decode("utf-8")

    stages = [
        "raw",
        "cleaned",
        "nonhuman",
        "ribofrac",
        "processed",
        "cladecounts",
        "humanviruses",
        "allmatches",
        "hvreads",
        "samplereads",
        "readlengths",
        "alignments2",
        "valreads",
        "tmpvalreads",
    ]
    short_stages = [
        "raw",
        "clean",
        "nh",
        "rf",
        "kraken",
        "cc",
        "hv",
        "am",
        "hvr",
        "sr",
        "rl",
        "al",
        "vr",
        "tvr",
    ]

    name_width = 21
    print(" " * name_width, end="\t")
    print(*short_stages, sep="\t")
    for delivery in sorted(deliveries):
        print(delivery)

        if True:
            if delivery in running_processes:
                color = COLOR_CYAN
            else:
                color = ""

            print(
                color
                + ("  " + delivery).ljust(name_width)
                + (COLOR_END if color else ""),
                end="",
                flush=True,
            )

            fully_processed_fname = work_fname(
                "deliveries", delivery, "fully_processed"
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

            s3_delivery_dir = "%s/%s" % (S3_BUCKET, delivery)
            metadata_dir = work_fname("deliveries", delivery, "metadata")

            stage_counters = defaultdict(Counter)  # sample -> stage -> count

            with open(os.path.join(metadata_dir, "metadata.tsv")) as inf:
                samples = [x.strip().split("\t")[0] for x in inf]

            prev = None
            for stage in stages:
                print("\t", end="", flush=True)

                if stage == "clean" and is_nanopore(args):
                    print("n/a", end="", flush=True)
                    continue

                seen = set()
                for fname in ls_s3_dir("%s/%s/" % (
                        s3_delivery_dir, full_s3_dirname(stage))):
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
    ("ribofrac", ribofrac),
    ("nonhuman", nonhuman),
    ("interpret", interpret),
    ("cladecounts", cladecounts),
    ("humanviruses", humanviruses),
    ("allmatches", allmatches),
    ("hvreads", hvreads),
    ("samplereads", samplereads),
    ("readlengths", readlengths),
    ("alignments2", alignments2),
    ("valreads", valreads),
    ("tmpvalreads", tmpvalreads),
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
        "--delivery", help="The ID of the delivery to process"
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
        help="Instead of running anything, just print the status of deliveries",
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

    if not args.status and not args.delivery:
        parser.print_help()
        exit(1)

    if args.status:
        print_status(args)
        return

    if not os.path.isdir(work_fname("deliveries", args.delivery)):
        raise Exception(
            "Delivery %s not found in %sdeliveries"
            % (args.delivery, WORK_ROOT)
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
