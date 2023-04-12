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

def get_accessions(args):
   if args.accession:
      return [args.accession]

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

def bioproject_config(args):
   with open(os.path.join(get_metadata_dir(args), "study.json")) as inf:
      return json.load(inf)

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
   if not bioproject_config(args)["is_paired_end"]:
      raise Exception("Only paired end sequencing currently supported")

   adapter_dir = os.path.join(THISDIR, "bioprojects", args.bioproject, "adapters")
   try:
      os.mkdir(adapter_dir)
   except FileExistsError:
      pass

   for accession in get_accessions(args):
      if exists_s3_prefix("%s/%s/%s/%s" % (
            S3_BUCKET, args.bioproject, dirname, accession)):
         continue

      with tempdir("adapter_removal", accession) as workdir:
         in1="in1.fastq.gz"
         in2="in2.fastq.gz"

         subprocess.check_call([
            "aws", "s3", "cp", "%s/%s/raw/%s_1.fastq.gz" % (
               S3_BUCKET, args.bioproject, accession), in1])
         subprocess.check_call([
            "aws", "s3", "cp", "%s/%s/raw/%s_2.fastq.gz" % (
               S3_BUCKET, args.bioproject, accession), in2])

         adapter1_fname = os.path.join(adapter_dir, "%s.fwd" % accession)
         adapter2_fname = os.path.join(adapter_dir, "%s.rev" % accession)

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
                "--basename", accession,
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

         for output in glob.glob("%s.*" % accession):
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

   for accession in get_accessions(args):
      for potential_input in available_inputs:
         if not potential_input.startswith(accession): continue
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

def viruscount(args):
   available_inputs = get_files(args, "processed")
   existing_outputs = get_files(args, "viruscounts")

   for accession in get_accessions(args):
      output = "%s.viruscounts.tsv" % accession
      if output in existing_outputs: continue

      inputs = [
         input_fname
         for input_fname in available_inputs
         if input_fname.startswith(accession)]
      if not inputs:
         continue

      with tempdir("viruscount", accession) as workdir:
         for input_fname in inputs:
            subprocess.check_call([
               "aws", "s3", "cp", "%s/%s/processed/%s" % (
                  S3_BUCKET, args.bioproject, input_fname), input_fname])

         try:
            check_call_shell(
               "cat %s.*.kraken2.tsv.gz | "
               "gunzip | "
               "grep -i virus | "
               "awk -F'\t' '{print $3}' | "
               "sort | uniq -c | sort -n | "
               "while read n rest ; do echo -e \"$n\\t$rest\" ; done | "
               "aws s3 cp - %s/%s/viruscounts/%s.viruscounts.tsv" % (
                  accession, S3_BUCKET, args.bioproject, accession))
         except subprocess.CalledProcessError:
            # Probably no viruses in the file
            pass

def allcounts(args):
   available_inputs = get_files(args, "processed")
   existing_outputs = get_files(args, "allcounts", min_size=100)

   for accession in get_accessions(args):
      output = "%s.tsv.gz" % accession
      if output in existing_outputs: continue

      if not any(x.startswith(accession) for x in available_inputs):
         continue

      subprocess.check_call([
         "./count_taxonomic_ids.sh", args.bioproject, accession])   
         
def humanviruses(args):
   human_viruses = {}
   with open(os.path.join(THISDIR, "human-viruses.tsv")) as inf:
      for line in inf:
         taxid, name = line.strip().split("\t")
         human_viruses[int(taxid)] = name

   available_inputs = get_files(args, "processed")
   existing_outputs = get_files(args, "humanviruses")

   for accession in get_accessions(args):
      output = "%s.humanviruses.tsv" % accession
      if output in existing_outputs: continue

      inputs = [
         input_fname
         for input_fname in available_inputs
         if input_fname.startswith(accession)]
      if not inputs:
         continue

      counts = Counter()

      for input_fname in inputs:
         with tempdir("humanviruses", accession) as workdir:
            subprocess.check_call([
               "aws", "s3", "cp", "%s/%s/processed/%s" % (
                  S3_BUCKET, args.bioproject, input_fname), input_fname])

            with gzip.open(input_fname, "rt") as inf:
               for line in inf:
                  taxid, = re.findall("[(]taxid ([0-9]+)[)]", line)
                  taxid = int(taxid)
                  if taxid in human_viruses:
                     counts[taxid] += 1

      with tempdir("humanviruses", accession) as workdir:
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

   for accession in get_accessions(args):
      output = "%s.allmatches.tsv" % accession
      if output in existing_outputs: continue

      inputs = [
         input_fname
         for input_fname in available_inputs
         if input_fname.startswith(accession)]
      if not inputs:
         continue

      with tempdir("allmatches", accession) as workdir:
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

   for accession in get_accessions(args):
      output = "%s.hvreads.json" % accession
      if output in existing_outputs: continue

      input_fname = "%s.allmatches.tsv" % accession
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
         if not cleaned_input.startswith(accession): continue
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

def qc_reads(args, accession, qc_info,
             available_raw, available_cleaned,
             available_processed, available_viruscounts):
   in1_fname = accession + "_1.fastq.gz"
   if in1_fname not in available_raw:
      return None

   return int(check_output_shell(
      "aws s3 cp %s/%s/raw/%s - | gunzip | grep ^@ | wc -l" % (
         S3_BUCKET, args.bioproject, in1_fname)))

def qc_cleaned_reads(args, accession, qc_info,
                     available_raw, available_cleaned,
                     available_processed, available_viruscounts):
   counts = {}
   for fname in available_cleaned:
      if fname.startswith(accession):
         slug = fname.replace("%s." % accession, "").replace(".gz", "")
         if slug == "settings": continue
         counts[slug] = int(check_output_shell(
            "aws s3 cp %s/%s/cleaned/%s - | gunzip | grep ^@ | wc -l" % (
               S3_BUCKET, args.bioproject, fname)))
   if not counts:
      return None

   return counts

def qc_cleaning_summary(args, accession, qc_info,
                        available_raw, available_cleaned,
                        available_processed,
                        available_viruscounts):
   if "reads" not in qc_info or "cleaned_reads" not in qc_info:
      return None

   cleaned_count = sum(
      count
      for (slug, count) in qc_info["cleaned_reads"].items()
      if not slug.startswith("pair2"))

   collapsed_count = sum(
      count
      for (slug, count) in qc_info["cleaned_reads"].items()
      if slug.startswith("collapsed"))

   summary = {}
   summary["dropped"] = qc_info["reads"] - cleaned_count
   summary["collapsed_fraction"] = collapsed_count / qc_info["reads"]

   return summary

def phred_to_q(phred_score):
   return ord(phred_score) - ord('!')

def average_quality(phred_counts):
   xs = []
   weights = []

   for phred_score, count in phred_counts.items():
      xs.append(phred_to_q(phred_score))
      weights.append(count)

   return round(np.average(xs, weights=weights))

def qc_post_cleaning(args, accession, qc_info,
                     available_raw, available_cleaned,
                     available_processed,
                     available_viruscounts):
   if "cleaned_reads" not in qc_info:
      return None

   lengths = []
   qualities = defaultdict(list)

   with tempdir("qc_post_cleaning", accession) as workdir:
      for fname in available_cleaned:
         if fname.startswith(accession):
            slug = fname.replace("%s." % accession, "").replace(".gz", "")
            slug_raw_qualities = []
            if "settings" in slug: continue
            subprocess.check_call([
               "aws", "s3", "cp", "%s/%s/cleaned/%s" % (
                  S3_BUCKET, args.bioproject, fname),fname])
            with gzip.open(fname, "rt") as inf:
               for (title, sequence, quality) in FastqGeneralIterator(inf):
                  if "collapsed" in slug:
                     while len(sequence) > len(lengths) - 1:
                        lengths.append(0)
                     lengths[len(sequence)] += 1

                  for pos, qval in enumerate(quality):
                     if pos > len(slug_raw_qualities) - 1:
                        slug_raw_qualities.append(Counter())
                     slug_raw_qualities[pos][qval] += 1
            qualities[slug] = [
               average_quality(phred_counts)
               for phred_counts in slug_raw_qualities]

   return {
      "lengths": lengths,
      "qualities": qualities,
   }

def qc(args):
   available_raw = get_files(args, "raw")
   available_cleaned = get_files(args, "cleaned", min_size=100)
   available_processed = get_files(args, "processed")
   available_viruscounts = get_files(args, "viruscounts")

   qc_dir = os.path.join(THISDIR, "bioprojects", args.bioproject, "qc")
   if not os.path.exists(qc_dir):
      os.mkdir(qc_dir)

   for accession in get_accessions(args):
      qc_fname = os.path.join(qc_dir, "%s.json" % accession)
      qc_info = {}
      if os.path.exists(qc_fname):
         with open(qc_fname) as inf:
            qc_info = json.load(inf)

      for qc_key, qc_fn in [
            ("reads", qc_reads),
            ("cleaned_reads", qc_cleaned_reads),
            ("cleaning_summary", qc_cleaning_summary),
            ("post_cleaning", qc_post_cleaning),
      ]:
         if qc_key not in qc_info:
            result = qc_fn(
               args,
               accession,
               qc_info,
               available_raw,
               available_cleaned,
               available_processed,
               available_viruscounts)
            if result is not None:
               qc_info[qc_key] = result
               with open(qc_fname, "w") as outf:
                  json.dump(qc_info, outf, indent=2, sort_keys=True)
                  outf.write("\n")

def print_status(args):
   if args.bioproject:
      bioprojects = [args.bioproject]
   else:
      bioprojects = [
         os.path.basename(x)
         for x in glob.glob(os.path.join(THISDIR, "bioprojects", "*"))]

   # Name -> Bioproject Accession -> Stage -> N/M
   info = defaultdict(dict)

   stages = ["raw", "cleaned", "processed", "allcounts",
             "viruscounts", "humanviruses",
             "allmatches", "hvreads"]
   short_stages = ["raw", "clean", "kraken", "ac", "vc", "hv", "am", "hvr"]

   for n, bioproject in enumerate(bioprojects):
      print("\rgathering status information %s/%s..." % (
         n+1, len(bioprojects)), end='', flush=True)

      metadata_dir = os.path.join(THISDIR, "bioprojects", bioproject, "metadata")
      if not os.path.exists(metadata_dir):
         continue

      with open(os.path.join(metadata_dir, "name.txt")) as inf:
         name = inf.read().strip()

      info[name][bioproject] = {}

      with open(os.path.join(metadata_dir, "metadata.tsv")) as inf:
         accessions = [x.strip().split("\t")[0] for x in inf]

      for accession in accessions:
         info[name][bioproject][accession] = Counter()

      s3_bioproject_dir = "%s/%s" % (S3_BUCKET, bioproject)

      for stage in stages:
         for fname in ls_s3_dir("%s/%s/" % (s3_bioproject_dir, stage)):
            for accession in accessions:
               if fname.startswith(accession):
                  info[name][bioproject][accession][stage] += 1


   print("\n")

   name_width=21
   print(" "*name_width, end='\t')
   print(*short_stages, sep='\t')

   for name in sorted(info):
      print(name)
      for bioproject in sorted(info[name]):
         row = [("  " + bioproject).ljust(name_width)]
         totals = Counter()

         for accession in info[name][bioproject]:
            for stage in stages:
               if info[name][bioproject][accession][stage]:
                  totals[stage] += 1

         prev = None
         for stage in stages:
            missing = prev is not None and totals[stage] < prev
            row.append("%s%s%s" % (
               COLOR_RED if missing else "",
               totals[stage],
               COLOR_END if missing else ""))
            prev = totals[stage]

         print("\t".join(row))


STAGES_ORDERED = []
STAGE_FNS = {}
for stage_name, stage_fn in [("clean", clean),
                             ("interpret", interpret),
                             ("allcounts", allcounts),
                             ("viruscount", viruscount),
                             ("humanviruses", humanviruses),
                             ("allmatches", allmatches),
                             ("hvreads", hvreads),
                             #("qc", qc),
                             ]:
   STAGES_ORDERED.append(stage_name)
   STAGE_FNS[stage_name] = stage_fn

def start():
   parser = argparse.ArgumentParser(
      description='Run the Metagenomic Sequencing Pipeline')

   parser.add_argument(
      '--bioproject', help='The ID of the bioproject to process')
   parser.add_argument(
      '--accession', default='',
      help='The ID of the sample to process.  Leave blank for all samples.')

   parser.add_argument(
      "--status", action='store_true',
      help='Instead of running anything, just print the status of bioprojects')

   parser.add_argument(
      '--stages',
      default=",".join(STAGES_ORDERED),
      help='Comma-separated list of stages to run.  Allowed stages: %s' % (
         ", ".join(repr(x) for x in STAGES_ORDERED)))

   args = parser.parse_args()

   if args.status:
      print_status(args)
      return

   if not args.bioproject:
      raise Exception("Must specify --bioproject BIOPROJECT")

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
