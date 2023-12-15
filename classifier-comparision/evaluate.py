#!/usr/bin/env python3

import os
import sys
import glob
import gzip
import json
import importlib
from collections import Counter, defaultdict

THIS_DIR = os.path.dirname(__file__)
MGS_PIPELINE_DIR = os.path.join(THIS_DIR, "..")
DASHBOARD_DIR = os.path.join(MGS_PIPELINE_DIR, "dashboard")

parents = {}
with open(os.path.join(DASHBOARD_DIR, "nodes.dmp")) as inf:
    for line in inf:
        child_taxid, parent_taxid, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        parents[int(child_taxid)] = int(parent_taxid)

human_viruses = {}
with open(os.path.join(MGS_PIPELINE_DIR, "human-viruses.tsv")) as inf:
    for line in inf:
        taxid, name = line.rstrip("\n").split("\t")
        human_viruses[int(taxid)] = name

def taxid_under(clade, taxid):
    while taxid not in [0, 1]:
        if taxid == clade:
            return True
        taxid = parents[taxid]
    return False

def is_bacterial(taxid):
    return taxid_under(2, taxid)

def load_data(sample):
    cdata = {}

    cdata["human_viruses"] = human_viruses
    cdata["dashboard_dir"] = DASHBOARD_DIR
    cdata["is_bacterial"] = is_bacterial

    cdata["hvreads"] = {}
    cdata["alignments"] = {
        "hv": defaultdict(list),
        "human": defaultdict(list),
    }
    cdata["processed"] = {}

    all_reads = set()

    with open(os.path.join(cdata["dashboard_dir"],
                           "hvreads",
                           "%s.hvreads.json" % sample)) as inf:
        hvr = json.load(inf)
        for read_id in hvr:
            cdata["hvreads"][read_id] = hvr[read_id]
            all_reads.add(read_id)

    for bowtie_db in cdata["alignments"]:
        with gzip.open(os.path.join(cdata["dashboard_dir"],
                                    "alignments",
                                    "%s.%s.alignments.tsv.gz" % (
                                        sample, bowtie_db)),
                       "rt") as inf:
            for line in inf:
                bits = line.removesuffix("\n").split("\t")
                read_id = bits[0]
                all_reads.add(read_id)
                cdata["alignments"][bowtie_db][read_id].append(bits)

    # This one is so big we process it streaming.
    for fname in glob.glob(
            os.path.join(cdata["dashboard_dir"],
                         "processed", "%s.*.kraken2.tsv.gz" % sample)):
        with gzip.open(fname, "rt") as inf:
            for line in inf:
                bits = line.removesuffix("\n").split("\t")
                read_id = bits[1]
                if read_id not in all_reads:
                    # All of the classifiers need at least an HVR or alignment
                    # hit, so we can speed things up a lot by ignoring any
                    # reads that don't have one of these.
                    continue
                cdata["processed"] = bits
                yield read_id, cdata

def read_id_to_sample_id(read_id):
    return read_id.removeprefix("M_").split(".")[0]

def start(*classifier_names):
    classifiers = [
        importlib.import_module(classifier_name)
        for classifier_name in classifier_names]

    with open("classified-reads.json") as inf:
        ground_truth = json.load(inf)

    gt_read_ids = set(ground_truth["yes"]) | set(ground_truth["no"])
    samples = set(read_id_to_sample_id(read_id)
                  for read_id in gt_read_ids)

    metrics = [Counter() for _ in classifiers]
    metric_names = [
        "totpos", "truepos", "trueneg", "falsepos", "falseneg"]

    all_notes = {}
    for sample in sorted(samples):
        #if sample != "SRR14530891":
        #    continue
        sample_notes = {}
        for read_id, cdata in load_data(sample):
            read_notes = {}
            if read_id in gt_read_ids:
                read_notes["hv"] = read_id in ground_truth["yes"]

            for classifier_name, classifier, metric in zip(
                    classifier_names, classifiers, metrics):
                out = classifier.classify(sample, read_id, cdata)
                if type(out) == type(True):
                    result = out
                    classifier_notes = {}
                else:
                    result, classifier_notes = out

                classifier_notes["hv"] = result

                if read_id in gt_read_ids:
                    is_hv = read_id in ground_truth["yes"]

                    status = {
                        (True, True): "truepos",
                        (False, False): "trueneg",
                        (False, True): "falsepos",
                        (True, False): "falseneg",
                    }[is_hv, result]
                    metric[status] += 1

                if result:
                    metric["totpos"] += 1

                read_notes[classifier_name] = classifier_notes

            log_notes = False
            if read_id in gt_read_ids:
                log_notes = True
            elif (read_notes.get("kraken-hits", {}).get("hv") !=
                  read_notes.get(
                      "bowtie-on-kraken-unclassified", {}).get("hv")):
                log_notes = True

            if read_notes and log_notes:
                sample_notes[read_id] = read_notes
        all_notes[sample] = sample_notes

    with open(os.path.join(THIS_DIR, "evaluation_notes.json"), "w") as outf:
        json.dump(all_notes, outf)

    classifier_column_width = max(
        len(classifier_name)
        for classifier_name in classifier_names)

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
