from Bio import SeqIO
from Bio.Seq import Seq
import sys

# Path to the file containing positions with variance
variance_file = sys.argv[2] + '/out_of_' + sys.argv[1] + '.output.vcalling.txt' # Update this with your file path
reference_fasta = sys.argv[2] + '/' + sys.argv[3]  # Update this with your reference FASTA file path
output_fasta =  sys.argv[2] + '/alt_' + sys.argv[1] + '.fasta' # Output file path

# Function to read positions and alternative values from the file
def read_positions_with_variance(file_path):
    positions_with_variance = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) == 2:
                pos_part, alt_part = parts
                position = int(pos_part.split(':')[1].strip())
                alt_value = alt_part.split(':')[1].strip()
                positions_with_variance.append((position, alt_value))
    return positions_with_variance

# Read positions with variance from the file
positions_with_variance = read_positions_with_variance(variance_file)

# Load the reference FASTA file
sequences = SeqIO.to_dict(SeqIO.parse(reference_fasta, "fasta"))

# Apply the changes to the sequences
for position, alt_value in positions_with_variance:
    for record_id, record in sequences.items():
        seq_list = list(record.seq)  # Convert to a list for mutability
        seq_list[position - 1] = alt_value  # Update the position with the alternative value (0-based index)
        record.seq = Seq(''.join(seq_list))  # Convert back to a Seq object and update the record's sequence

# Write the updated sequences to a new FASTA file
with open(output_fasta, 'w') as output_handle:
    SeqIO.write(sequences.values(), output_handle, "fasta")

print(f"Updated reference FASTA file saved to {output_fasta}")
