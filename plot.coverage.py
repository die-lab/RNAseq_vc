import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import sys

def list_coverage_files(directory):
    coverage_files = [f for f in os.listdir(directory) if f.endswith('coverage.txt')]
    return coverage_files

directory = sys.argv[1]
coverage_files = list_coverage_files(directory)
for file in coverage_files:
    print(file)

def plot_coverage(coverage_files):
    plt.figure(figsize=(15, 8))
    
    for coverage_file in coverage_files:
        data = pd.read_csv(coverage_file, sep='\t', header=None, names=['chrom', 'pos', 'coverage'])
        plt.plot(data['pos'], data['coverage'], label=os.path.basename(coverage_file).replace('.coverage.txt', ''))
        plt.ylim(0,10000)
        plt.xlabel('Position')
        plt.ylabel('Coverage')
        plt.legend()
        plt.title('Coverage along the Reference Genome')
        plt.savefig(os.path.basename(coverage_file).replace('.txt','.png'))
        plt.show()
        plt.clf()

def main(coverage_files):
    plot_coverage(coverage_files)

# Example usage
main(coverage_files)
