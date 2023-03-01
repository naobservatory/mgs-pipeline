#!/usr/bin/env python3

import os
import glob
import json
import time
import atexit
import argparse
import tempfile
import subprocess

MAX_PARALLEL=32

S3_BUCKET="s3://nao-mgs"
THISDIR=os.path.abspath(os.path.dirname(__file__))

running_subprocesses = set()
def kill_subprocesses():
   for process in running_subprocesses:
      process.kill()
atexit.register(kill_subprocesses)

def relative_fname(fname):
   return os.path.join(THISDIR, fname)

def get_metadata_dir(args):
   return relative_fname("studies/%s/metadata" % args.study)

def get_accessions(args):
   with open(os.path.join(get_metadata_dir(args), "metadata.tsv")) as inf:
      return [line.split("\t")[0]
              for line in inf]

def study_config(args):
   with open(os.path.join(get_metadata_dir(args), "study.json")) as inf:
      return json.load(inf)

def get_adapters(args, adapter_dir):
   try:
      os.mkdir(adapter_dir)
   except FileExistsError:
      pass

   for accession in get_accessions(args):
      for in_fname, direction in [("in1.fastq.gz", "fwd"),
                                  ("in2.fastq.gz", "rev")]:
         out_fname = os.path.join(
            adapter_dir, "%s.%s" % (accession, direction))
         if (os.path.exists(out_fname) and
             os.path.getsize(out_fname) > 0): continue

         # Run at most MAX_PARALLEL simultaneous processes.
         while len(running_subprocesses) >= MAX_PARALLEL:
            for process in running_subprocesses:
               if process.poll() is not None:
                  running_subprocesses.remove(process)
                  break
            else:
               time.sleep(0.1)

         running_subprocesses.add(subprocess.Popen(
            "aws s3 cp '%s/%s/raw/%s' - | "
            "%s/determine_adapter.py /dev/stdin %s > "
            "%s" % (
               S3_BUCKET, args.study, in_fname,
               THISDIR, direction,
               out_fname),
            # All of these arguments are trusted and doing this
            # manually with Popen is really painful.
            shell=True))

   # Wait for them all to finish.
   while running_subprocesses:
      for process in running_subprocesses:
         if process.poll() is not None:
            running_subprocesses.remove(process)
            break
      else:
         time.sleep(0.1)

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

def clean(args):
   if not study_config(args)["is_paired_end"]:
      raise Exception("Only paired end sequencing currently supported")

   adapter_dir = os.path.join(THISDIR, "studies", args.study, "adapters")
   get_adapters(args, adapter_dir)

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

         with open(os.path.join(adapter_dir, "%s.fwd" % accession)) as inf:
            adapter1 = inf.read().strip()
         with open(os.path.join(adapter_dir, "%s.rev" % accession)) as inf:
            adapter2 = inf.read().strip()

         subprocess.check_call([
            "AdapterRemoval",
            "--file1", in1,
            "--file2", in2,
            "--basename", accession,
            "--trimns",
            "--trimqualities",
            "--collapse",
            "--adapter1", adapter1,
            "--adapter2", adapter2,
         ])

         for output in glob.glob("%s.*" % accession):
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
      ls_s3_dir("%s/%s/viruscounts/" % (S3_BUCKET, args.study)))

   for accession in get_accessions(args):
      output = "%s.viruscounts.tsv" % accession
      compressed_output = output + ".gz"
      if compressed_output in existing_outputs: continue

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
