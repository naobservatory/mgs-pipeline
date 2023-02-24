#!/usr/bin/env python3

import os
import json
import argparse
import tempfile
import subprocess
import determine_adapter

S3_BUCKET="s3://nao-mgs"
THISDIR=os.path.dirname(__file__)
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
   
def clean(args):
   if not study_config(args)["is_paired_end"]:
      raise Exception("Only paired end sequencing currently supported")
   
   for accession in get_accessions(args):
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

         adapter1 = determine_adapter.get_adapter(in1, "fwd")
         adapter2 = determine_adapter.get_adapter(in2, "rev")

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
