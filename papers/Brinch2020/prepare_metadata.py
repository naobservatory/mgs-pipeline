#!/usr/bin/env python3


sample_to_details = {}
with open("table-s1.tsv") as inf:
    for line in inf:
        if not line.strip(): continue
        line = line[:-1]

        if line.startswith("seque"): continue

        bits = line.split("\t")

        sample_to_details[bits[7]] = bits[15], bits[17]

for project in ["PRJEB34633", "PRJEB13832"]:
    with open("../../bioprojects/%s/metadata/metadata_raw.tsv" %
              project) as inf:
        with open("../../bioprojects/%s/metadata/metadata.tsv" %
                  project, "w") as outf:
            for line in inf:
                sample = line.strip()
                if sample not in sample_to_details:
                    print("Unknown sample %s" % sample)
                    continue
                loc, date = sample_to_details[sample]
                outf.write("%s\t%s\t%s\n" % (sample, loc, date))
                
        
