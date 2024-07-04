import pandas as pd
import glob
from collections import defaultdict
import sys
import os

def read_alternative_values(file_path):
    # Read positions and alternative values from the file
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            pos = int(parts[1])
            alt_value = parts[4]
            data.append((pos, alt_value))
    return data

def combine_files(file_paths):
    combined_data = {}
    all_positions = set()
    
    for file_path in file_paths:
        sample_name = os.path.basename(file_path).replace('.alternative.calls.tab', '')
        data = read_alternative_values(file_path)
        for pos, alt_value in data:
            if pos not in combined_data:
                combined_data[pos] = {}
            combined_data[pos][sample_name] = alt_value
            all_positions.add(pos)
    
    # Create a DataFrame with all positions and alternative values from each file
    sorted_positions = sorted(all_positions)
    columns = ['Position'] + [os.path.basename(f).replace('.alternative.calls.tab', '') for f in file_paths]
    merged_df = pd.DataFrame(columns=columns)
    
    for pos in sorted_positions:
        row = [pos]
        for file_path in file_paths:
            sample_name = os.path.basename(file_path).replace('.alternative.calls.tab', '')
            row.append(combined_data.get(pos, {}).get(sample_name, '-'))
        merged_df.loc[len(merged_df)] = row
    
    return merged_df

def main(input_directory):
    # List all files in the directory ending with .alternative.calls.tab
    file_paths = glob.glob(input_directory + '/*.alternative.calls.tab')
    merged_df = combine_files(file_paths)
    
    # Write the merged DataFrame to a new file
    output_file = 'merged_alternative_values.csv'
    merged_df.to_csv(output_file, index=False)

    print(f"Output written to {output_file}")

# Example usage:
# python merge_alternative_values.py /path/to/input_directory 
if __name__ == "__main__":
    input_directory = sys.argv[1]  # Input directory containing .alternative.calls.tab files
    main(input_directory)
