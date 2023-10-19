import pandas as pd
import csv

def categorize_spots(num_spots):
    if num_spots <= 5:
        return 'low'
    elif 6 <= num_spots <= 15:
        return 'mid'
    else:
        return 'high'

def categorize_area(spot_area):
    if 0 <= spot_area < 0.1:
        return '0-0.1'
    elif 0.1 <= spot_area < 0.2:
        return '0.1-0.2'
    elif 0.2 <= spot_area < 0.3:
        return '0.2-0.3'
    elif 0.3 <= spot_area < 0.4:
        return '0.3-0.4'
    elif 0.4 <= spot_area < 0.5:
        return '0.4-0.5'
    else:
        return '0.5+'

def process_data(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames

        fieldnames.extend(['Spot Category', 'Area Category'])

        # Create a CSV writer for the output file
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            num_spots = int(row['num_spots'])
            spot_area = float(row['spot_area'])
            row['Spot Category'] = categorize_spots(num_spots)
            row['Area Category'] = categorize_area(spot_area)
            writer.writerow(row)

if __name__ == "__main__":
    input_file = 'test_data.csv'
    output_file = 'output.csv'
    process_data(input_file, output_file)

