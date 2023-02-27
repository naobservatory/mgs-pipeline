#!/usr/bin/env python3

import os
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
      in1="%s_1.fastq.gz" % accession
      in2="%s_2.fastq.gz" % accession

      for in_fname, direction in [(in1, "fwd"), (in2, "rev")]:
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

         in1="%s_1.fastq.gz" % accession
         in2="%s_2.fastq.gz" % accession

         subprocess.check_call([
            "aws", "s3", "cp", "%s/%s/raw/%s" % (
               S3_BUCKET, args.study, in1), "."])
         subprocess.check_call([
            "aws", "s3", "cp", "%s/%s/raw/%s" % (
               S3_BUCKET, args.study, in2), "."])

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
            "--interleaved-output",
            "--adapter1", adapter1,
            "--adapter2", adapter2,
         ])

         for output_suffix  in [
               "collapsed", "collapsed.truncated", "discarded",
               "paired.truncated", "singleton.truncated"]:
            output = "%s.%s" % (accession, output_suffix)
            subprocess.check_call(["gzip", output])
            output = output + ".gz"

            subprocess.check_call([
               "aws", "s3", "cp", output, "%s/%s/cleaned/" % (
                  S3_BUCKET, args.study)])

def start():
   parser = argparse.ArgumentParser(
      description='Run the Metagenomic Sequencing Pipeline')
   parser.add_argument(
      'study', help='The ID of the study to process')

   args = parser.parse_args()

   clean(args)

if __name__ == "__main__":
   start()
