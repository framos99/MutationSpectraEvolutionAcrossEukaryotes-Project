import json
import os
import sys
import getopt

def help(error_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -S => <string> species name REQUIRED
    -Q => <string> sequence to extract phis for REQUIRED
ASSUMPTIONS
    *phi summary jsons are located in "/project/voight_subrate2/fabramos/100_Species/2024_BAYMER_OUT_POST_SUMMARIES"
""")
    sys.exit(error_num)

def extract_phi_values(SPECIES, SEQUENCE):
    json_file = os.path.join(f"/project/voight_subrate2/fabramos/100_Species/2024_BAYMER_OUT_POST_SUMMARIES/{SPECIES}.non_genic.ALL.final_posterior_summaries.json")
    # Load the JSON data from the file
    with open(json_file, 'r') as file:
        data = json.load(file)

    try:
        # Extract phi values for the specified sequence
        layer = len(SEQUENCE) - 1
        prob_values = data["ALL"][str(layer)]["p_vec"][SEQUENCE]["mean"]
        formatted_prob_values = "\t".join(map(str, prob_values))
        print(f"{SPECIES}\t{SEQUENCE}\t{formatted_prob_values}")
    except KeyError:
        print(f"Sequence {SEQUENCE} not found in the JSON file.")

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "S:Q:")
    except getopt.GetoptError:
        print("Error: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)

    ## Required arguments
    try:
        SPECIES = str(options_dict['-S'])
        SEQUENCE = str(options_dict['-Q'])

    except KeyError as E:
        missing_argument = E.args[0]
        print(f"Error: The required argument '{missing_argument}' is missing.")
        help()
    extract_phi_values(SPECIES=SPECIES, SEQUENCE=SEQUENCE)

if __name__ == "__main__":
    main(sys.argv[1:])
