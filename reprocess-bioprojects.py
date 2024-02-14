#!/usr/bin/env python3

# Runs the pipeline with the provided arguments across all bioprojects.  Useful
# when there's a new stage, or a stage needs to be re-run.
#
# Usage: ./reprocess-bioprojects.sh \
#    --max-jobs <N> --log-prefix <LP> [--sample-level] -- arguments-for-run
#
# Examples:
#    ./reprocess-bioprojects.sh \
#        --max-jobs 12 --log-prefix rl -- --stages readlengths
#
# You can pass --bioprojects A,B,C to run on only a subset of projects.
#
#    ./reprocess-bioprojects.sh \
#        --bioprojects PRJNA729801 --max-jobs 12 --log-prefix rl \
#        -- --stages readlengths
#
# And if you run --sample-level it will invoke run.py for each sample instead
# of for each bioproject, allowing more parallelism but also more overhead.
#
#    ./reprocess-bioprojects.sh \
#        --bioprojects PRJNA729801 --sample-level --max-jobs 12 \
#        --log-prefix rl -- --stages readlengths

import os
import re
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


def prepare_job(bioproject, log_prefix, sample, run_args):
    logfile = "%s/%s.%s.%s" % (log_dir, log_date, log_prefix, bioproject)
    if sample:
        logfile = "%s.%s" % (logfile, sample)
        run_args = list(run_args) + ["--sample", sample]

    return logfile, ["./run.py", "--bioproject", bioproject, *run_args]


def run_job(job):
    logfile, cmd = job
    with open(logfile, "w") as outf:
        result = subprocess.run(cmd, stdout=outf, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            outf.write("ERROR: %s\n" % (result.returncode))

def get_sample_priority(sample):
    m = re.findall(r"L00\d$", sample)
    if not m:
        return "A"
    m, = m
    return m

def parallelize(config, bioprojects, run_args):
    job_queue = []

    for bioproject in bioprojects:
        args = run_args[:]
        if bioproject in restricted_bioprojects:
            args.append("--restricted")
            root_dir = restricted_dir
        elif bioproject in regular_bioprojects:
            root_dir = "."
        else:
            raise Exception("Unknown bioproject %r" % bioproject)

        if config.sample_level:
            prioritized_samples = []
            with open(os.path.join(root_dir, "bioprojects", bioproject,
                                   "metadata", "metadata.tsv")) as inf:
                for line in inf:
                    sample = line.strip().split()[0]
                    prioritized_samples.append(
                        (get_sample_priority(sample), sample))
            prioritized_samples.sort()
            for priority, sample in prioritized_samples:
                job_queue.append(prepare_job(
                    bioproject, config.log_prefix, sample, args))
        else:
            job_queue.append(prepare_job(
                bioproject, config.log_prefix, None, args))

    with ThreadPoolExecutor(max_workers=config.max_jobs) as executor:
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
    parser.add_argument(
        "--sample-level",
        action="store_true",
        help="Parallelize at the sample level instead of the bioproject level")

    config = parser.parse_args(our_args)

    if config.bioprojects:
        bioprojects = config.bioprojects.split(",")
    else:
        bioprojects = regular_bioprojects + restricted_bioprojects

    if "--memory-mapping" in run_args:
        subprocess.check_call(["./prepare-shm-kraken.sh"])

    parallelize(config, bioprojects, run_args)


if __name__ == "__main__":
    start()
