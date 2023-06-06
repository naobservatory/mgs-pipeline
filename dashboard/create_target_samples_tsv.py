import json
import csv

# TODO: Make script target studies by providing a list of papers, querying
# "metadata_papers.json" to pull bioproject ids from, and then pull sample ids
# from metadata_bioprojects.json

# load the json data
with open("metadata_samples.json") as f:
    samples = json.load(f)

# with open("metadata_papers.json") as f:
#     papers = json.load(f)

# with open("metadata_bioprojects.json") as f:
#     bioprojects = json.load(f)

# initialize an empty set of headers
headers = set()
rows = []

# iterate through the json samples and prepare rows
for accession, values in samples.items():
    # add the new keys to the headers
    headers.update(values.keys())
    # create a row dictionary where missing values are represented as empty strings
    row = {header: values.get(header, "") for header in headers}
    # If the key "state" doesn't exist yet, or its value is "", skip this row
    if "state" not in row:
        continue
    if row["state"] == "":
        continue
    row["accession"] = accession
    rows.append(row)

# write to tsv file
with open("metadata_target_samples.tsv", "w", newline="") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=["accession"] + list(headers))
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
