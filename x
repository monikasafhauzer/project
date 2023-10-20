import matplotlib.pyplot as plt
import numpy as np

# Make a plot to visualize cells
# 0-5 (low), 6-15(mid), >15(high) cell A vs cell B
# areas cell A vs cell B

plt.figure(figsize=(10, 6))

# Group by 'Spot Category'
cell_a_spot_counts = cell_a_df['Spot Category'].value_counts()
cell_b_spot_counts = cell_b_df['Spot Category'].value_counts()

# Define the categories
categories = ['low', 'mid', 'high']

x = range(len(categories))
plt.bar(x, cell_a_spot_counts[categories], width=0.4, label='Cell A')
plt.bar([i + 0.4 for i in x], cell_b_spot_counts[categories], width=0.4, label='Cell B')

# Add labels to the plot
plt.xlabel('Spot Category')
plt.ylabel('Count')
plt.title('Spot Number Comparison between Cell A and Cell B')
plt.xticks([i + 0.2 for i in x], categories)

plt.legend()

# Save the plot as a PNG file
plt.savefig("pbody.png")

# Display the plot
plt.show()
