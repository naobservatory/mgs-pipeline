#!/usr/bin/env python3

# Run prepare_metadata.sh, which calls this script.

import sys


def start(raw_metadata_in, parsed_metadata_out):
    data = []
    with open(raw_metadata_in) as inf:
        for line in inf:
            line = line.strip()

            _, _, fastq_ftps, alias = line.split("\t")

            if fastq_ftps == "fastq_ftp":
                continue  # skip header line

            enriched = not alias.endswith("_unenriched")
            alias = alias.removesuffix("_unenriched")

            alias = alias.removesuffix("_INF")

            plant, month, day, year = alias.split("_")

            if len(year) == 2:
                year = "20%s" % year

            date = "%s-%s-%s" % (year, month.zfill(2), day.zfill(2))

            for fastq_ftp in fastq_ftps.split(";"):
                if not fastq_ftp:
                    continue

                filename = fastq_ftp.split("/")[-1]

                data.append([plant, date, filename, "1" if enriched else "0"])

    data.sort()

    with open(parsed_metadata_out, "w") as outf:
        outf.write("\t".join(["filename", "date", "plant", "is_enriched"]) + "\n")
        for plant, date, filename, is_enriched in data:
            outf.write("\t".join([filename, date, plant, is_enriched]) + "\n")


if __name__ == "__main__":
    start(*sys.argv[1:])
