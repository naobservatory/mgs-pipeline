#!/usr/bin/env python3

import os
import sys
import gzip
import json
import importlib
from collections import Counter, defaultdict

MGS_PIPELINE_DIR = os.path.join(os.path.dirname(__file__), "..")

def collect_classification_data(read_ids):
    cdata = {}

    human_viruses = {}
    with open(os.path.join(MGS_PIPELINE_DIR, "human-viruses.tsv")) as inf:
        for line in inf:
            taxid, name = line.rstrip("\n").split("\t")
            human_viruses[int(taxid)] = name

    cdata["human_viruses"] = human_viruses
    cdata["dashboard_dir"] = os.path.join(
        MGS_PIPELINE_DIR, "dashboard")

    sample_ids = set(read_id_to_sample_id(read_id)
                     for read_id in read_ids)

    cdata["hvreads"] = {}
    cdata["alignments"] = {
        "hv": defaultdict(list),
        "human": defaultdict(list),
    }
    for sample in sample_ids:
        with open(os.path.join(cdata["dashboard_dir"],
                               "hvreads",
                               "%s.hvreads.json" % sample)) as inf:
            hvr = json.load(inf)
            for read_id in read_ids:
                if read_id in hvr:
                    cdata["hvreads"][read_id] = hvr[read_id]

        for bowtie_db in cdata["alignments"]:
            with gzip.open(os.path.join(cdata["dashboard_dir"],
                                        "alignments",
                                        "%s.%s.alignments.tsv.gz" % (
                                            sample, bowtie_db)),
                           "rt") as inf:
                for line in inf:
                    bits = line.removesuffix("\n").split("\t")
                    read_id = bits[0]
                    if read_id in read_ids:
                        cdata["alignments"][bowtie_db][read_id].append(bits)

    return cdata

def read_id_to_sample_id(read_id):
    return read_id.removeprefix("M_").split(".")[0]

def start(*classifier_names):
    classifiers = [
        importlib.import_module(classifier_name)
        for classifier_name in classifier_names]

    with open("classified-reads.json") as inf:
        ground_truth = json.load(inf)

    all_read_ids = set(ground_truth["yes"]) | set(ground_truth["no"])

    cdata = collect_classification_data(all_read_ids)
    
    metrics = [Counter() for _ in classifiers]

    print("read_id", *classifier_names, sep="\t")
    for read_id in sorted(all_read_ids):
        sample = read_id_to_sample_id(read_id)
        is_hv = read_id in ground_truth["yes"]

        results = [classifier.classify(sample, read_id, cdata)
                   for classifier in classifiers]

        for metric, result in zip(metrics, results):
            status = {
                (True, True): "truepos",
                (False, False): "trueneg",
                (False, True): "falsepos",
                (True, False): "falseneg",
            }[is_hv, result]
            metric[status] += 1

    classifier_column_width = max(
        len(classifier_name)
        for classifier_name in classifier_names)
    metric_names = list(sorted(metrics[0]))
    print(" "*classifier_column_width,
          *metric_names,
          sep="  ")
    for classifier_name, classifier_metrics in zip(
            classifier_names, metrics):
        print(classifier_name.ljust(classifier_column_width),
              *[str(classifier_metrics[metric_name]).rjust(
                  len(metric_name))
                for metric_name in metric_names],
              sep="  ")

if __name__ == "__main__":
    start(*sys.argv[1:])
