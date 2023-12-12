#!/usr/bin/env python3

# Runs the pipeline with the provided arguments across all bioprojects.  Useful
# when there's a new stage, or a stage needs to be re-run.
#
# Usage: ./reprocess-bioprojects.sh \
#    --max-jobs <N> --log-prefix <LP> -- arguments-for-run
#
# Example: ./reprocess-bioprojects.sh \
#    --max-jobs 12 --log-prefix rl -- --stages readlengths
#
# You can pass --bioprojects A,B,C to run on only a subset of projects.

import os
import sys
import datetime
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor

log_date = datetime.datetime.now().date().isoformat()
log_dir = "log"
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

regular_bioprojects = os.listdir("bioprojects")

restricted_bioprojects = []
restricted_dir = os.path.join("..", "mgs-restricted")
if os.path.exists(restricted_dir):
    restricted_bioprojects = os.listdir(
        os.path.join(restricted_dir, "bioprojects")
    )


def prepare_job(bioproject, log_prefix, run_args):
    logfile = "%s/%s.%s.%s" % (log_dir, log_date, log_prefix, bioproject)
    return logfile, ["./run.py", "--bioproject", bioproject, *run_args]


def run_job(job):
    logfile, cmd = job
    with open(logfile, "w") as outf:
        result = subprocess.run(cmd, stdout=outf, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            outf.write("ERROR: %s\n" % (result.returncode))


def parallelize(max_jobs, log_prefix, bioprojects, run_args):
    job_queue = []

    for bioproject in bioprojects:
        args = run_args[:]
        if bioproject in restricted_bioprojects:
            args.append("--restricted")
        elif bioproject in regular_bioprojects:
            pass
        else:
            raise Exception("Unknown bioproject %r" % bioproject)

        job_queue.append(prepare_job(bioproject, log_prefix, args))

    with ThreadPoolExecutor(max_workers=max_jobs) as executor:
        for job in job_queue:
            executor.submit(run_job, job)


def start():
    argv = sys.argv[1:]
    if "--" not in argv:
        raise Exception("Use -- to separate arguments to ./run.py.")

    our_args = argv[: argv.index("--")]
    run_args = argv[argv.index("--") + 1 :]

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--max-jobs",
        metavar="N",
        type=int,
        required=True,
        help="maximum number of jobs to run at once",
    )
    parser.add_argument(
        "--log-prefix",
        required=True,
        help="Log prefix, for storing this run under log/",
    )
    parser.add_argument(
        "--bioprojects",
        help="The IDs of the bioproject to process, comma separated",
    )
    args = parser.parse_args(our_args)

    if args.bioprojects:
        bioprojects = args.bioprojects.split(",")
    else:
        bioprojects = regular_bioprojects + restricted_bioprojects

    parallelize(args.max_jobs, args.log_prefix, bioprojects, run_args)


if __name__ == "__main__":
    start()
