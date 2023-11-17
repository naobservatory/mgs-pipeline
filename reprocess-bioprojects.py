#!/usr/bin/env python3

# Runs the pipeline with the provided arguments across all bioprojects.  Useful
# when there's a new stage, or a stage needs to be re-run.
#
# Usage: ./reprocess-bioprojects.sh <max_jobs> <log prefix> [./run.py arguments]
#
# Examples:
#   ./reprocess-bioprojects.sh 12 rl --stages readlengths
#   ./reprocess-bioprojects.sh 2 nr --skip-stages ribocounts,ribofrac
#

import os
import sys
import datetime
import subprocess
from concurrent.futures import ThreadPoolExecutor

log_date = datetime.datetime.now().date().isoformat()
log_dir = "log"
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

def prepare_job(bioproject, log_prefix, *args):
    logfile = "%s/%s.%s.%s" % (log_dir, log_date, log_prefix, bioproject)
    return logfile, ["./run.py", "--bioproject", bioproject, *args]

def run_job(job):
    logfile, cmd = job
    with open(logfile, "w") as outf:
        result = subprocess.run(cmd, stdout=outf, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            outf.write("ERROR: %s\n" % (result.returncode))

def start(max_jobs, log_prefix, *args):
    assert max_jobs.isdigit()
    max_jobs = int(max_jobs)
    assert not log_prefix.startswith("--")

    regular_bioprojects = os.listdir("bioprojects")

    restricted_bioprojects = []
    restricted_dir = os.path.join("..", "mgs-restricted")
    if os.path.exists(restricted_dir):
        restricted_bioprojects = os.listdir(
            os.path.join(restricted_dir, "bioprojects"))

    job_queue = []

    for bioproject in regular_bioprojects:
        job_queue.append(prepare_job(bioproject, log_prefix, *args))
    for bioproject in restricted_bioprojects:
        job_queue.append(prepare_job(
            bioproject, log_prefix, "--restricted", *args))

    with ThreadPoolExecutor(max_workers=max_jobs) as executor:
        for job in job_queue:
            executor.submit(run_job, job)

if __name__ == "__main__":
   start(*sys.argv[1:])
