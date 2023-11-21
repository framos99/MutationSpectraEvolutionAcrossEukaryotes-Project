#!/usr/bin/Python
from Bio import SeqIO
import os
import subprocess
import sys
import gzip
import getopt

def help(error_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -I => <FASTA> Assembly file REQUIRED
    -T => <string> NCBI taxonomic group for repetitive element families REQUIRED
    -W => <string> Working directory path REQUIRED
    -O_dir => <string> RMSK_output_dir (within working_directory) REQUIRED
ASSUMPTIONS
    * RepeatMasker/4.1.5 module is loaded
""")
    sys.exit(error_num)

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "I:T:W:O:")
    except getopt.GetoptError:
        print("Error: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)

    ## Required arguments
    try:
        input_fasta = str(options_dict['-I'])
        NCBI_TAXONOMY = str(options_dict['-T'])
        working_directory = str(options_dict['-W'])
        RMSK_output_directory = str(options_dict['-O'])

    except KeyError as E:
        missing_argument = E.args[0]
        print(f"Error: The required argument '{missing_argument}' is missing.")
        help()
    print("Acceptable Inputs Given")

    split_fasta_and_run_repeatmasker(input_fasta = input_fasta, NCBI_TAXONOMY = NCBI_TAXONOMY, working_directory = working_directory, RMSK_output_directory = RMSK_output_directory)

def split_fasta_and_run_repeatmasker(input_fasta, NCBI_TAXONOMY, working_directory, RMSK_output_directory):
    #list of sequence IDs already processed
    Already_processed_sequences = RMSK_do_not_run_list(working_directory, RMSK_output_directory)
    # RepeatMasker command
    repeatmasker_command = "bsub RepeatMasker -species {NCBI_TAXONOMY} -e hmmer -libdir /misc/appl/RepeatMasker-4.1.5/Libraries/ -hmmer_dir /misc/appl/hmmer-3.1b2/bin/ {individual_sequence_fasta}"
    # File paths
    input_fasta_path = os.path.join(working_directory, f"Assembly/{input_fasta}")
    # Create the output directory if it doesn't exist
    output_directory = os.path.join(working_directory, f"Assembly/{RMSK_output_directory}")
    os.makedirs(output_directory, exist_ok=True)
    #Assembly sequence IDs
    sequence_records_file = os.path.join(working_directory, f"Assembly/assembly_sequenceIDs.txt")

    # Parse the input FASTA file
    if input_fasta_path.endswith(".gz"):
        with gzip.open(input_fasta_path, "rt") as fasta_file:
            # iterate through each sequence
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id not in Already_processed_sequences:
                    # Create a separate FASTA file for each sequence
                    output_file = os.path.join(output_directory, f"{record.id}.fasta")
                    with open(output_file, "w") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")
                    #sequence records list
                    with open(sequence_records_file, 'a') as sequenceID_file:
                        sequenceID_file.write(f"{record.id}\n")
                    # Run RepeatMasker on the individual FASTA file
                    individual_sequence_fasta = os.path.abspath(output_file)
                    command = repeatmasker_command.format(NCBI_TAXONOMY=NCBI_TAXONOMY, individual_sequence_fasta=individual_sequence_fasta)
                    subprocess.run(command, shell=True, check=True)
    else:
        with open(input_fasta_path, "r") as fasta_file:
            #iterate through each sequence
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id not in Already_processed_sequences:
                    # Create a separate FASTA file for each sequence
                    output_file = os.path.join(output_directory, f"{record.id}.fasta")
                    with open(output_file, "w") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")
                    #sequence records list
                    with open(sequence_records_file, 'a') as sequenceID_file:
                        sequenceID_file.write(f"{record.id}\n")
                    # Run RepeatMasker on the individual FASTA file
                    individual_sequence_fasta = os.path.abspath(output_file)
                    command = repeatmasker_command.format(NCBI_TAXONOMY=NCBI_TAXONOMY, individual_sequence_fasta=individual_sequence_fasta)
                    subprocess.run(command, shell=True, check=True)
                    
def RMSK_do_not_run_list(working_directory, RMSK_output_directory):
    sequence_records_file = os.path.join(working_directory, f"Assembly/assembly_sequenceIDs.txt")
    Already_processed_sequences = []
    if os.path.isfile(sequence_records_file):
        with open(sequence_records_file, 'r') as file:
            records_list = [line.strip() for line in file.readlines()]
        for RECORD in records_list:
            output_directory = os.path.join(working_directory, f"Assembly/{RMSK_output_directory}")
            RMSK_out_file_path = os.path.join(output_directory, f"{RECORD}.fasta.out")
            if os.path.exists(RMSK_out_file_path):
                Already_processed_sequences.append(RECORD)
    return Already_processed_sequences

if __name__ == "__main__":
    main(sys.argv[1:])