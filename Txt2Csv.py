import csv

# Open the text file
with open("RGCAULevels.txt", "r") as txt_file:
    # Read the lines from the text file
    lines = txt_file.readlines()

# Open a CSV file for writing
with open("RGCAULevels.csv", "w", newline="") as csv_file:
    # Create a CSV writer object
    writer = csv.writer(csv_file, delimiter=',')
    
    # Write the header
    writer.writerow(["Gene", "GC", "AU", "Length", "Polar%", "Non-Polar%"])
    
    # Write the data from the text file to the CSV file
    for line in lines[1:]:  # Skip the header line
        fields = line.strip().split("\t")
        writer.writerow(fields)

