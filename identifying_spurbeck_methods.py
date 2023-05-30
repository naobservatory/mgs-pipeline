import csv

# Define the methods dictionary
methods = {
    "AB": "1",
    "C": "2",
    "D": "3",
    "EFGH": "4",
    "IJ": "J",
}

# Open the input and output CSV files
with open('input.csv', 'r') as input_file, open('output.csv', 'w', newline='') as output_file:
    # Create CSV reader and writer objects
    reader = csv.reader(input_file)
    writer = csv.writer(output_file)

    # Read the header row and add a new "Method Value" column
    header = next(reader)
    header.append("Method Value")
    writer.writerow(header)

    # Loop through each row in the input CSV file
    for row in reader:
        # Get the value in the "Method" column
        method = row[2]

        # Look up the corresponding value in the methods dictionary
        method_value = methods.get(method, "")

        # Add the method value to the row and write it to the output CSV file
        row.append(method_value)
        writer.writerow(row)

print("Output file created successfully!")