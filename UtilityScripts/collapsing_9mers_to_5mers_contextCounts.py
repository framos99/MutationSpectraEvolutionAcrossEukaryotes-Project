import pandas as pd
import argparse

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description='Collapse DNA sequences to middle 5 nucleotides and sum counts.')
parser.add_argument('input_file', help='Input text file with counts and sequences')
parser.add_argument('output_file', help='Output file to save the collapsed sequences')
args = parser.parse_args()

# Load the data
data = pd.read_csv(args.input_file, sep='\t', header=None, names=['count', 'sequence'])

# Extract the middle 5 nucleotides
data['context'] = data['sequence'].str[2:7]

# Group by the middle sequence and sum the counts
collapsed_data = data.groupby('context', as_index=False)['count'].sum()

# Save the collapsed data
collapsed_data.to_csv(args.output_file, sep='\t', index=False)
print(f"Collapsed data saved to {args.output_file}")
