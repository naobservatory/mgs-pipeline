#!/usr/bin/env python3

import os
import glob
import json
import time
import atexit
import argparse
import tempfile
import subprocess

S3_BUCKET="s3://nao-mgs"
THISDIR=os.path.abspath(os.path.dirname(__file__))

def relative_fname(fname):
   return os.path.join(THISDIR, fname)

def get_metadata_dir(args):
   return relative_fname("studies/%s/metadata" % args.study)

def get_accessions(args):
   if args.accession:
      return [args.accession]

   with open(os.path.join(get_metadata_dir(args), "metadata.tsv")) as inf:
      return [line.split("\t")[0]
              for line in inf]

def study_config(args):
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
      date, time, size, fname = line.split()
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

def clean(args):
   if not study_config(args)["is_paired_end"]:
      raise Exception("Only paired end sequencing currently supported")

   adapter_dir = os.path.join(THISDIR, "studies", args.study, "adapters")
   try:
      os.mkdir(adapter_dir)
   except FileExistsError:
      pass

   for accession in get_accessions(args):
      if exists_s3_prefix("%s/%s/cleaned/%s" % (
            S3_BUCKET, args.study, accession)):
         continue

      with tempfile.TemporaryDirectory() as workdir:
         print("Handling %s in %s" % (accession, workdir))
         os.chdir(workdir)

         in1="in1.fastq.gz"
         in2="in2.fastq.gz"

         subprocess.check_call([
            "aws", "s3", "cp", "%s/%s/raw/%s_1.fastq.gz" % (
               S3_BUCKET, args.study, accession), in1])
         subprocess.check_call([
            "aws", "s3", "cp", "%s/%s/raw/%s_2.fastq.gz" % (
               S3_BUCKET, args.study, accession), in2])

         adapter1_fname = os.path.join(adapter_dir, "%s.fwd" % accession)
         adapter2_fname = os.path.join(adapter_dir, "%s.rev" % accession)

         if (not os.path.exists(adapter1_fname) or
             not os.path.exists(adapter2_fname)):
            get_adapters(in1, in2, adapter1_fname, adapter2_fname)

         with open(adapter1_fname) as inf:
            adapter1 = inf.read().strip()
         with open(adapter2_fname) as inf:
            adapter2 = inf.read().strip()

         subprocess.check_call([
            "AdapterRemoval",
            "--file1", in1,
            "--file2", in2,
            "--basename", accession,
            "--trimns",
            "--trimqualities",
            "--collapse",
            "--threads", "4",
            "--adapter1", adapter1,
            "--adapter2", adapter2,
         ])

         for output in glob.glob("%s.*" % accession):
            # TODO: this could be streaming
            subprocess.check_call(["gzip", output])
            output = output + ".gz"

            subprocess.check_call([
               "aws", "s3", "cp", output, "%s/%s/cleaned/" % (
                  S3_BUCKET, args.study)])

def interpret(args):
   available_inputs = set(
      ls_s3_dir("%s/%s/cleaned/" % (S3_BUCKET, args.study),
                # tiny files are empty; ignore them
                min_size=100))
   existing_outputs = set(
      ls_s3_dir("%s/%s/processed/" % (S3_BUCKET, args.study)))

   for accession in get_accessions(args):
      for potential_input in available_inputs:
         if not potential_input.startswith(accession): continue
         if ".settings." in potential_input: continue

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

         with tempfile.TemporaryDirectory() as workdir:
            print("Handling %s in %s" % (", ".join(inputs), workdir))
            os.chdir(workdir)

            for input_fname in inputs:
               subprocess.check_call([
                  "aws", "s3", "cp", "%s/%s/cleaned/%s" % (
                     S3_BUCKET, args.study, input_fname), input_fname])

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
                  S3_BUCKET, args.study)])

def viruscount(args):
   available_inputs = set(
      ls_s3_dir("%s/%s/processed/" % (S3_BUCKET, args.study)))
   existing_outputs = set(
      ls_s3_dir("%s/%s/viruscounts/" % (S3_BUCKET, args.study),
                # empty file outputs mean something went wrong in an earlier
                # iteration, so let's try again.
                min_size=1))

   for accession in get_accessions(args):
      output = "%s.viruscounts.tsv" % accession
      if output in existing_outputs: continue

      with tempfile.TemporaryDirectory() as workdir:
         print("Handling %s in %s" % (accession, workdir))
         os.chdir(workdir)

         for input_fname in available_inputs:
            if not input_fname.startswith(accession): continue

            subprocess.check_call([
               "aws", "s3", "cp", "%s/%s/processed/%s" % (
                  S3_BUCKET, args.study, input_fname), input_fname])

         subprocess.check_call(
            "cat %s.*.kraken2.tsv.gz | "
            "gunzip | "
            "grep -i virus | "
            "awk -F'\t' '{print $3}' | "
            "sort | uniq -c | sort -n | "
            "while read n rest ; do echo -e \"$n\\t$rest\" ; done | "
            "aws s3 cp - %s/%s/viruscounts/%s.viruscounts.tsv" % (
               accession, S3_BUCKET, args.study, accession),
            shell=True)

STAGES_ORDERED = []
STAGE_FNS = {}
for stage_name, stage_fn in [("clean", clean),
                             ("interpret", interpret),
                             ("viruscount", viruscount)]:
   STAGES_ORDERED.append(stage_name)
   STAGE_FNS[stage_name] = stage_fn

def start():
   parser = argparse.ArgumentParser(
      description='Run the Metagenomic Sequencing Pipeline')
   parser.add_argument(
      '--study', help='The ID of the study to process')
   parser.add_argument(
      '--accession', default='',
      help='The ID of the sample to process.  Leave blank for all samples.')

   parser.add_argument(
      '--stages',
      default=",".join(STAGES_ORDERED),
      help='Comma-separated list of stages to run.  Allowed stages: %s' % (
         ", ".join(repr(x) for x in STAGES_ORDERED)))

   args = parser.parse_args()
   selected_stages = args.stages.split(",")
   for selected_stage in selected_stages:
      if selected_stage not in STAGE_FNS:
         raise Exception("Unknown stage %r" % selected_stage)

   for stage in STAGES_ORDERED:
      if stage in selected_stages:
         STAGE_FNS[stage](args)

if __name__ == "__main__":
   start()
