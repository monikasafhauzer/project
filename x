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

# Load the processed data from the output file
data = pd.read_csv('output.csv')

# Sort data for Cells A and B by num_spots
cell_a_data = data[data['Cell'] == 'Cell A'].sort_values('num_spots')
cell_b_data = data[data['Cell'] == 'Cell B'].sort_values('num_spots')

# Create a bar plot for num_spots
plt.figure(figsize=(10, 6))

plt.bar(cell_a_data['Cell'], cell_a_data['num_spots'], label='Cell A')
plt.bar(cell_b_data['Cell'], cell_b_data['num_spots'], label='Cell B')

plt.xlabel('Cell')
plt.ylabel('num_spots')
plt.title('Comparison of Cells A and B by num_spots')
plt.legend()
plt.show()

# Sort data for Cells A and B by spot_area
cell_a_data = data[data['Cell'] == 'Cell A'].sort_values('spot_area')
cell_b_data = data[data['Cell'] == 'Cell B'].sort_values('spot_area')

# Create a bar plot for spot_area
plt.figure(figsize=(10, 6))

plt.bar(cell_a_data['Cell'], cell_a_data['spot_area'], label='Cell A')
plt.bar(cell_b_data['Cell'], cell_b_data['spot_area'], label='Cell B')

plt.xlabel('Cell')
plt.ylabel('spot_area')
plt.title('Comparison of Cells A and B by spot_area')
plt.legend()
plt.show()
