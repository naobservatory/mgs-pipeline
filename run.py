#!/usr/bin/env python3

import re
import os
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
from collections import Counter
from collections import defaultdict
from Bio.SeqIO.QualityIO import FastqGeneralIterator

S3_BUCKET="s3://nao-mgs"
THISDIR=os.path.abspath(os.path.dirname(__file__))

COLOR_RED = '\x1b[0;31m'
COLOR_GREEN = '\x1b[0;32m'
COLOR_CYAN = '\x1b[0;36m'
COLOR_END = '\x1b[0m'

def check_call_shell(cmd):
   # Unlike subprocess.check_call, if any member of the pipeline fails then
   # this fails too.
   subprocess.check_call([
      'bash', '-c', 'set -o  pipefail; %s' % cmd])

def check_output_shell(cmd):
   # Unlike subprocess.check_output, if any member of the pipeline fails then
   # this fails too.
   return subprocess.check_output([
      'bash', '-c', 'set -o  pipefail; %s' % cmd])

def relative_fname(fname):
   return os.path.join(THISDIR, fname)

def get_metadata_dir(args):
   return relative_fname("bioprojects/%s/metadata" % args.bioproject)

def get_samples(args):
   if args.sample:
      return [args.sample]

   with open(os.path.join(get_metadata_dir(args), "metadata.tsv")) as inf:
      return [line.strip().split("\t")[0]
              for line in inf]

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
      subprocess.check_call(["aws", "s3", "ls", s3_path],
                            stdout=subprocess.DEVNULL)
      return True  # exit code 0 if present
   except subprocess.CalledProcessError as e:
      if e.returncode == 1:
         return False # exit code 1 if absent

      raise # any other exit code means something else is wrong

def ls_s3_dir(s3_dir, min_size=0):
   try:
      cmd_out = subprocess.check_output(["aws", "s3", "ls", s3_dir])
   except subprocess.CalledProcessError as e:
      if e.returncode == 1:
         return [] # exit code 1 if absent or empty
      raise # any other exit code means something is wrong

   for line in cmd_out.split(b"\n"):
      if not line.strip():
         continue
      try:
         date, time, size, fname = line.split()
      except ValueError:
         print(line)
         print(s3_dir)
         raise

      if int(size) < min_size: continue
      yield fname.decode('utf-8')

def get_adapters(in1, in2, adapter1_fname, adapter2_fname):
   output = subprocess.check_output([
      "AdapterRemoval",
      "--file1", in1,
      "--file2", in2,
      "--identify-adapters",
      "--threads", "4"])
   output = output.decode("utf-8")

   for line in output.split("\n"):
      if "--adapter1:" in line:
         adapter1 = line.replace("--adapter1:", "").strip()
      elif "--adapter2:" in line:
         adapter2 = line.replace("--adapter2:", "").strip()

   for adapter, fname in [[adapter1, adapter1_fname],
                          [adapter2, adapter2_fname]]:
      if not all(x in 'ACTGN' for x in adapter) or len(adapter) < 20:
         print(output)
         raise Exception("Invalid adapter %r for %r and %r" % (
            adapter, in1, in2))
      with open(fname, 'w') as outf:
         outf.write(adapter)

def adapter_removal(args, dirname, trim_quality, collapse):
   adapter_dir = os.path.join(THISDIR, "bioprojects", args.bioproject, "adapters")
   try:
      os.mkdir(adapter_dir)
   except FileExistsError:
      pass

   for sample in get_samples(args):
      if exists_s3_prefix("%s/%s/%s/%s" % (
            S3_BUCKET, args.bioproject, dirname, sample)):
         continue

      with tempdir("adapter_removal", sample) as workdir:
         in1="in1.fastq.gz"
         in2="in2.fastq.gz"

         subprocess.check_call([
            "aws", "s3", "cp", "%s/%s/raw/%s_1.fastq.gz" % (
               S3_BUCKET, args.bioproject, sample), in1])
         subprocess.check_call([
            "aws", "s3", "cp", "%s/%s/raw/%s_2.fastq.gz" % (
               S3_BUCKET, args.bioproject, sample), in2])

         adapter1_fname = os.path.join(adapter_dir, "%s.fwd" % sample)
         adapter2_fname = os.path.join(adapter_dir, "%s.rev" % sample)

         if (not os.path.exists(adapter1_fname) or
             not os.path.exists(adapter2_fname)):
            get_adapters(in1, in2, adapter1_fname, adapter2_fname)

         with open(adapter1_fname) as inf:
            adapter1 = inf.read().strip()
         with open(adapter2_fname) as inf:
            adapter2 = inf.read().strip()

         cmd = ["AdapterRemoval",
                "--file1", in1,
                "--file2", in2,
                "--basename", sample,
                "--threads", "4",
                "--adapter1", adapter1,
                "--adapter2", adapter2,
                "--gzip"]

         if trim_quality:
            cmd.extend(["--trimns",
                        "--trimqualities"])
         if collapse:
            cmd.append("--collapse")

         subprocess.check_call(cmd)

         for output in glob.glob("%s.*" % sample):
            subprocess.check_call([
               "aws", "s3", "cp", output, "%s/%s/%s/" % (
                  S3_BUCKET, args.bioproject, dirname)])

def clean(args):
   adapter_removal(args, "cleaned", trim_quality=True, collapse=True)

def rmadapter(args):
   adapter_removal(args, "noadapters", trim_quality=False, collapse=False)

def get_files(args, dirname, min_size=1):
   return set(ls_s3_dir("%s/%s/%s/" % (S3_BUCKET, args.bioproject, dirname),
                        min_size=min_size))

def interpret(args):
   available_inputs = get_files(args, "cleaned",
                                # tiny files are empty; ignore them
                                min_size=100)
   existing_outputs = get_files(args, "processed")

   for sample in get_samples(args):
      for potential_input in available_inputs:
         if not potential_input.startswith(sample): continue
         if ".settings" in potential_input: continue

         output = potential_input.replace(".gz", ".kraken2.tsv")
         inputs = [potential_input]
         if ".pair1." in output:
            output = output.replace(".pair1.", ".")
            inputs.append(potential_input.replace(".pair1.", ".pair2."))
         elif ".pair2" in output:
            # We handle pair1 and pair2 together.
            continue

         compressed_output = output + ".gz"
         if compressed_output in existing_outputs: continue

         with tempdir("interpret", ", ".join(inputs)) as workdir:
            for input_fname in inputs:
               subprocess.check_call([
                  "aws", "s3", "cp", "%s/%s/cleaned/%s" % (
                     S3_BUCKET, args.bioproject, input_fname), input_fname])

            kraken_cmd = [
               "/home/ec2-user/kraken2-install/kraken2",
               "--db", "/home/ec2-user/kraken-db/",
               "--use-names",
               "--output", output]
            if len(inputs) > 1:
               kraken_cmd.append("--paired")
            kraken_cmd.extend(inputs)

            subprocess.check_call(kraken_cmd)
            subprocess.check_call(["gzip", output])
            subprocess.check_call([
               "aws", "s3", "cp", compressed_output, "%s/%s/processed/" % (
                  S3_BUCKET, args.bioproject)])

def cladecounts(args):
   available_inputs = get_files(args, "processed")
   existing_outputs = get_files(args, "cladecounts", min_size=100)

   for sample in get_samples(args):
      output = "%s.tsv.gz" % sample
      if output in existing_outputs: continue

      if not any(x.startswith(sample) for x in available_inputs):
         continue

      subprocess.check_call([
         "./count_clades.sh", args.bioproject, sample])

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
      if output in existing_outputs: continue

      inputs = [
         input_fname
         for input_fname in available_inputs
         if input_fname.startswith(sample)]
      if not inputs:
         continue

      counts = Counter()

      for input_fname in inputs:
         with tempdir("humanviruses", sample) as workdir:
            subprocess.check_call([
               "aws", "s3", "cp", "%s/%s/processed/%s" % (
                  S3_BUCKET, args.bioproject, input_fname), input_fname])

            with gzip.open(input_fname, "rt") as inf:
               for line in inf:
                  taxid, = re.findall("[(]taxid ([0-9]+)[)]", line)
                  taxid = int(taxid)
                  if taxid in human_viruses:
                     counts[taxid] += 1

      with tempdir("humanviruses", sample) as workdir:
         with open(output, "w") as outf:
            for taxid, count in sorted(counts.items()):
               outf.write("%s\t%s\t%s\n" % (taxid, count, human_viruses[taxid]))

         subprocess.check_call([
            "aws", "s3", "cp", output, "%s/%s/humanviruses/%s" % (
               S3_BUCKET, args.bioproject, output)])

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
      if output in existing_outputs: continue

      inputs = [
         input_fname
         for input_fname in available_inputs
         if input_fname.startswith(sample)]
      if not inputs:
         continue

      with tempdir("allmatches", sample) as workdir:
         kept = []
         for input_fname in inputs:
            subprocess.check_call([
               "aws", "s3", "cp", "%s/%s/processed/%s" % (
                  S3_BUCKET, args.bioproject, input_fname), input_fname])

            with gzip.open(input_fname, "rt") as inf:
               for line in inf:
                  keep = False
                  try:
                     taxid_matches = line.strip().split("\t")[4]
                     for taxid_match in taxid_matches.split(" "):
                        taxid, n_kmers = taxid_match.split(":")
                        if taxid == "A": continue # ambiguous nucleotide
                        if taxid == "|": continue # paired end transition
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

         subprocess.check_call([
            "aws", "s3", "cp", output, "%s/%s/allmatches/%s" % (
               S3_BUCKET, args.bioproject, output)])

def hvreads(args):
   available_inputs = get_files(args, "allmatches")
   available_cleaned_inputs = get_files(
      args, "cleaned",
      # tiny files are empty; ignore them
      min_size=100)

   existing_outputs = get_files(args, "hvreads")

   for sample in get_samples(args):
      output = "%s.hvreads.json" % sample
      if output in existing_outputs: continue

      input_fname = "%s.allmatches.tsv" % sample
      if input_fname not in available_inputs: continue

      all_matches = [
         x.strip().split("\t")
         for x in subprocess.check_output([
               "aws", "s3", "cp", "%s/%s/allmatches/%s" % (
                  S3_BUCKET, args.bioproject, input_fname), "-"]).decode(
                     'utf-8').split("\n")
         if x.strip()]

      seqs = {} # seqid -> kraken, fwd, rev
      for _, seq_id, _, _, kraken_details in all_matches:
         seqs[seq_id] = [kraken_details]

      for cleaned_input in sorted(available_cleaned_inputs):
         if not cleaned_input.startswith(sample): continue
         if ".settings" in cleaned_input: continue

         with tempdir("hvreads", cleaned_input) as workdir:
            subprocess.check_call([
               "aws", "s3", "cp", "%s/%s/cleaned/%s" % (
                  S3_BUCKET, args.bioproject,
                  cleaned_input), cleaned_input])

            with gzip.open(cleaned_input, "rt") as inf:
               for (title, sequence, quality) in FastqGeneralIterator(inf):
                  seq_id = title.split()[0]
                  if seq_id in seqs:
                     seqs[seq_id].append(sequence)

      with tempdir("hvreads", output) as workdir:
         with open(output, "w") as outf:
            json.dump(seqs, outf, sort_keys=True)
         subprocess.check_call([
            "aws", "s3", "cp", output, "%s/%s/hvreads/%s" % (
               S3_BUCKET, args.bioproject, output)])

def phred_to_q(phred_score):
   return ord(phred_score) - ord('!')

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
         for x in glob.glob(os.path.join(THISDIR, "bioprojects", "*"))]

   running_processes = subprocess.check_output(["ps", "aux"]).decode("utf-8")

   stages = ["raw", "cleaned", "processed", "cladecounts", "humanviruses",
             "allmatches", "hvreads"]
   short_stages = ["raw", "clean", "kraken", "cc", "hv", "am", "hvr"]

   papers_to_projects = defaultdict(list) # paper -> [project]
   for bioproject in bioprojects:
      metadata_dir = os.path.join(THISDIR, "bioprojects", bioproject, "metadata")
      if not os.path.exists(metadata_dir): continue
      with open(os.path.join(metadata_dir, "name.txt")) as inf:
         name = inf.read().strip()
      papers_to_projects[name].append(bioproject)

   name_width=21
   print(" "*name_width, end='\t')
   print(*short_stages, sep='\t')
   for paper, bioprojects in sorted(papers_to_projects.items()):
      print(paper)

      for bioproject in bioprojects:
         if bioproject in running_processes:
            color = COLOR_CYAN
         else:
            color = ""

         print(color +
               ("  " + bioproject).ljust(name_width) +
               (COLOR_END if color else ""), end="", flush=True)

         fully_processed = os.path.exists(os.path.join(
            THISDIR, "bioprojects", bioproject, "fully_processed"))

         s3_bioproject_dir = "%s/%s" % (S3_BUCKET, bioproject)
         metadata_dir = os.path.join(THISDIR, "bioprojects", bioproject,
                                     "metadata")

         stage_counters = defaultdict(Counter) # sample -> stage -> count

         with open(os.path.join(metadata_dir, "metadata.tsv")) as inf:
            samples = [x.strip().split("\t")[0] for x in inf]

         prev = None
         for stage in stages:
            print("\t", end="", flush=True)

            if fully_processed and stage != "raw":
               print("-", end="", flush=True)
               continue

            seen = set()
            for fname in ls_s3_dir("%s/%s/" % (s3_bioproject_dir, stage)):
               for sample in samples:
                  if fname.startswith(sample):
                     seen.add(sample)
            missing = prev is not None and len(seen) < prev

            color = ""
            if missing:
               color = COLOR_RED
            elif fully_processed:
               color = COLOR_GREEN

            print("%s%s%s" % (
               color,
               len(seen),
               COLOR_END if color else ""),
                  end="", flush=True)
            prev = len(seen)

         print()

STAGES_ORDERED = []
STAGE_FNS = {}
for stage_name, stage_fn in [("clean", clean),
                             ("interpret", interpret),
                             ("cladecounts", cladecounts),
                             ("humanviruses", humanviruses),
                             ("allmatches", allmatches),
                             ("hvreads", hvreads),
                             ]:
   STAGES_ORDERED.append(stage_name)
   STAGE_FNS[stage_name] = stage_fn

def start():
   parser = argparse.ArgumentParser(
      description='Run the Metagenomic Sequencing Pipeline')

   parser.add_argument(
      '--bioproject', help='The ID of the bioproject to process')
   parser.add_argument(
      '--sample', default='',
      help='The SRA run accession of the sample to process.  Leave blank '
      'for all samples.')

   parser.add_argument(
      "--status", action='store_true',
      help='Instead of running anything, just print the status of bioprojects')

   parser.add_argument(
      '--stages',
      default=",".join(STAGES_ORDERED),
      help='Comma-separated list of stages to run.  Allowed stages: %s' % (
         ", ".join(repr(x) for x in STAGES_ORDERED)))

   args = parser.parse_args()

   if not args.status and not args.bioproject:
      parser.print_help()
      exit(1)

   if args.status:
      print_status(args)
      return

   if not os.path.isdir("bioprojects/%s" % args.bioproject):
      raise Exception(
         "Bioproject %s not found in bioprojects/" % args.bioproject)

   selected_stages = args.stages.split(",")
   for selected_stage in selected_stages:
      if selected_stage not in STAGE_FNS:
         raise Exception("Unknown stage %r" % selected_stage)

   for stage in STAGES_ORDERED:
      if stage in selected_stages:
         STAGE_FNS[stage](args)

if __name__ == "__main__":
   start()
