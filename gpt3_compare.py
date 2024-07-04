import pandas as pd
import glob
from collections import defaultdict
import sys

# Step 1: Read and parse each file
file_pattern = sys.argv[2] + '/*.alternative.calls.tab' # Update this with your file path pattern
files = glob.glob(file_pattern)

# Initialize a dictionary to store positions and their reference/alternative values
position_data = defaultdict(lambda: {'ref_values': set(), 'alt_values': set(), 'file_counts': defaultdict(int)})

for file in files:
    df = pd.read_csv(file, delimiter='\t', header=None, 
                     names=['reference', 'position', 'dot', 'reference_value', 'alternative_value', 'quality', 'filter', 'info', 'format', 'sample'])
    for index, row in df.iterrows():
        position = row['position']
        ref_value = row['reference_value']
        alt_value = row['alternative_value']
        
        position_data[position]['ref_values'].add(ref_value)
        position_data[position]['alt_values'].add(alt_value)
        position_data[position]['file_counts'][alt_value] += 1

# Function to identify positions with variance
def find_positions_with_variance(threshold):
    positions_with_variance = []
    for position, values in position_data.items():
        if len(values['ref_values']) == 1:
            for alt_value, count in values['file_counts'].items():
                if count >= threshold:
                    positions_with_variance.append((position, alt_value))
    return positions_with_variance

# Specify the threshold
try:
    threshold = int(sys.argv[1])  # Adjust this value as needed
except ValueError:
    threshold = 1
    
# Get positions with variance
positions_with_variance = find_positions_with_variance(threshold)

# Output the positions with variance
print(f"#Positions with variance across at least {threshold} files:")
for position, alt_value in positions_with_variance:
    print(f"Position: {position}, Alternative Value: {alt_value}")
