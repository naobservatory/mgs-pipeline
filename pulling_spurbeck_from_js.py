import json
import re
import csv

# Open the file and load the data
with open('data.js', 'r') as file:
    data = file.read()

# Use a regex to find the sample_metadata object
match = re.search(r'sample_metadata\s*=\s*({.*?})\s*(;|$)', data, re.DOTALL)

if match:
    json_data = match.group(1)

    # Parse the JSON data
    sample_metadata = json.loads(json_data)

    # Prepare data for CSV
    csv_data = []
    for sample, properties in sample_metadata.items():
        if properties.get('state') == "Ohio":  # 'state' field considered
            csv_data.append([sample] + list(properties.values()))

    # Write data to CSV
    with open('output.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Adjust the header according to your data structure
        writer.writerow(['Sample', 'Country', 'County', 'Date', 'Enrichment', 'Fine Location', 'Location', 'Method', 'Reads', 'State'])
        writer.writerows(csv_data)


methods = {
    "AB" : "1",
    "C" : "2",
    "D" : "3",
    "EFGH": "4",
    "IJ": "J",
}