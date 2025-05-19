import json
import os
import sys
import getopt

def help(error_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -S => <string> species name REQUIRED
    -L => <string> layer name REQUIRED
ASSUMPTIONS
    *phi summary jsons are located in "/project/voight_subrate2/fabramos/100_Species/2024_BAYMER_OUT_POST_SUMMARIES"
""")
    sys.exit(error_num)

def extract_ind_values(SPECIES, LAYER):
    json_file = os.path.join(f"/project/voight_subrate2/fabramos/100_Species/2024_BAYMER_OUT_POST_SUMMARIES/{SPECIES}.non_genic.ALL.final_posterior_summaries.json")
    output_file = f"/project/voight_subrate2/fabramos/100_Species/IndValues_for_BaymerLayers/{SPECIES}_Layer_{LAYER}_BaymerInds.txt"
    # Load the JSON data from the file
    with open(output_file, 'w') as output:
        with open(json_file, 'r') as file:
            data = json.load(file)

        try:
            # Extract all sequences in the specified layer
            sequences = data["ALL"][LAYER]["ind"].keys()

            # Iterate through all sequences and extract ind values
            for SEQUENCE in sequences:
                ind_values = data["ALL"][LAYER]["ind"][SEQUENCE]

                if len(SEQUENCE) % 2 == 0:
                    middle_index = (len(SEQUENCE)//2) - 1
                else:
                    middle_index = len(SEQUENCE)//2

                middle_nucleotide = SEQUENCE[middle_index]
                if middle_nucleotide == "A":
                    for i, nucleotide in enumerate(["C","G","T"]):            
                        modified_sequence = f"{SEQUENCE[0:middle_index]}[{middle_nucleotide}>{nucleotide}]{SEQUENCE[middle_index+1:]}"
                        ind_value = ind_values[i] if i < len(ind_values) else None
                        output.write(f"{modified_sequence}\t{ind_value}\n")
                elif middle_nucleotide == "C":
                    for i, nucleotide in enumerate(["A","G","T"]):        
                        modified_sequence = f"{SEQUENCE[0:middle_index]}[{middle_nucleotide}>{nucleotide}]{SEQUENCE[middle_index+1:]}"
                        ind_value = ind_values[i] if i < len(ind_values) else None
                        output.write(f"{modified_sequence}\t{ind_value}\n")

        except KeyError:
            print(f"Layer {LAYER} not found or no sequences in the JSON file.")

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "S:L:")
    except getopt.GetoptError:
        print("Error: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)

    ## Required arguments
    try:
        SPECIES = str(options_dict['-S'])
        LAYER = str(options_dict['-L'])

    except KeyError as E:
        missing_argument = E.args[0]
        print(f"Error: The required argument '{missing_argument}' is missing.")
        help()
    
    extract_ind_values(SPECIES=SPECIES, LAYER=LAYER)

if __name__ == "__main__":
    main(sys.argv[1:])
