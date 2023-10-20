import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv

# Define the data processing functions
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
    df = pd.read_csv(input_file)
    df['Spot Category'] = df['num_spots'].apply(categorize_spots)
    df['Area Category'] = df['spot_area'].apply(categorize_area)
    df.to_csv(output_file, index=False)

# Process the data
input_file = 'test_data.csv'
output_file = 'output.csv'
process_data(input_file, output_file)

# Create a plot using the processed data
df = pd.read_csv(output_file)
plt.figure(figsize=(10, 6))

cell_a_spot_counts = df[df['Cell'] == 'Cell A']['Spot Category'].value_counts()
cell_b_spot_counts = df[df['Cell'] == 'Cell B']['Spot Category'].value_counts()

categories = ['low', 'mid', 'high']

x = range(len(categories))
plt.bar(x, cell_a_spot_counts[categories], width=0.4, label='Cell A')
plt.bar([i + 0.4 for i in x], cell_b_spot_counts[categories], width=0.4, label='Cell B')

plt.xlabel('Spot Category')
plt.ylabel('Count')
plt.title('Spot Number Comparison between Cell A and Cell B')
plt.xticks([i + 0.2 for i in x], categories)

plt.legend()

# Save the plot as a PNG file
plt.savefig("pbody.png")

# Display the plot
plt.show()
