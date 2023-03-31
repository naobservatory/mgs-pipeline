#!/usr/bin/env python3

import os
from collections import Counter

# samples

# paper, sample, country, location, date, reads
samples = []

def remove_newline(s):
    assert s.endswith("\n")
    return s[:-1]

for bioproject in os.listdir("bioprojects"):
    with open("bioprojects/%s/metadata/name.txt" % bioproject) as inf:
        paper = inf.read().strip()

    # not processed yet
    if paper == "Hjelmso 2019": continue

    # not actually all this country, since it's airplanes
    if paper == "Petersen 2015": continue    

    na, country, location, date = {
        # TODO: supplementary data has sampling dates and plane origins
        "Petersen 2015": ("DNA", "Denmark", "Airplanes", "Summer 2013"),
        "Munk 2022": ("DNA", None, None, None),
        # TODO: supplementary data has sampling dates and boroughs
        "Maritz 2019": ("DNA", "USA", "New York City", "2014-2015"),
        # TODO: integrate WTP level data and dates from metadata.txt
        "Rothman 2021": ("RNA", "USA", "Los Angeles", "Fall 2020"),
    }[paper]

    with open("bioprojects/%s/metadata/metadata.tsv" % bioproject) as inf:
        for line in inf:
            bits = remove_newline(line).split("\t")
            if paper == "Munk 2022":
                sample, country, location, date = bits
            elif paper == "Rothman 2021":
                sample, date, wtp = bits
            else:
                sample = bits[0]

            country = {
                "Macedonia, Republic of": "North Macedonia",
                "Macedonia": "North Macedonia",
                "United States of America": "USA",
                "Iran, Islamic Republic of": "Iran",
                "C?te d'Ivoire": "Cote d'Ivoire",
                "Czech Rep": "Czech Republic",
                "Kosova": "Kosovo",
            }.get(country, country)
                
            samples.append([paper, bioproject, sample, country,
                            location, date, na])

for record in samples:
    bioproject = record[1]
    sample = record[2]
    with open("bioprojects/%s/metadata/%s.n_reads" % (
            bioproject, sample)) as inf:
        n_reads = int(inf.read().strip())
    record.append(n_reads)

countries = set(record[3] for record in samples)

#for record in samples:
#    print(*record, sep="\t")

for target_country in sorted(countries):
    paper_count = Counter()
    dna_total = 0
    for paper, bioproject, sample, country, location, date, na, count in samples:
        if country != target_country: continue
        if na != "DNA": continue

        dna_total += count
        paper_count[paper] += count

    print(dna_total, target_country, ";".join(sorted(paper_count.keys())),
          sep="\t")
            
