import json
import os
import sys
import getopt
import numpy as np

def help(error_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -S => <string> species name REQUIRED
ASSUMPTIONS
    *summary jsons are located in "/project/voight_subrate2/fabramos/100_Species/2024_BAYMER_OUT_POST_SUMMARIES"
""")
    sys.exit(error_num)

def extract_EVEN_prob_values(SPECIES, LAYER):
    json_file = os.path.join(f"/project/voight_subrate2/fabramos/100_Species/2024_BAYMER_OUT_POST_SUMMARIES/{SPECIES}.non_genic.EVEN.final_posterior_summaries.json")
    output_file = f"/project/voight_subrate2/fabramos/100_Species/Baymer_Layer_probabilities_inCosmicAnalysisFormat/{SPECIES}_Layer_{LAYER}_EVEN_BaymerProbabilities.txt"
    EVEN_fittedRates = []
    EVEN_CpG_to_T_rates = []
    EVEN_C_to_T_rates = []
    EVEN_CHG_to_T_rates = []
    # Load the JSON data from the file
    with open(output_file, 'w') as output:
        output.write("Type\tProbability_Value\n")
        with open(json_file, 'r') as file:
            data = json.load(file)

        try:
            # Extract all sequences in the specified layer
            sequences = data["EVEN"][LAYER]["p_vec"].keys()

            # Iterate through all sequences and extract phi values
            for SEQUENCE in sequences:
                phi_values = data["EVEN"][LAYER]["p_vec"][SEQUENCE]["mean"]

                if len(SEQUENCE) % 2 == 0:
                    middle_index = (len(SEQUENCE)//2) - 1
                else:
                    middle_index = len(SEQUENCE)//2

                middle_nucleotide = SEQUENCE[middle_index]
                if middle_nucleotide == "A":
                    for i, nucleotide in enumerate(["C","G","T"]):            
                        modified_sequence = f"{SEQUENCE[0:middle_index]}[{middle_nucleotide}>{nucleotide}]{SEQUENCE[middle_index+1:]}"
                        phi_value = phi_values[i] if i < len(phi_values) else None
                        EVEN_fittedRates = np.append(EVEN_fittedRates, phi_value)
                        output.write(f"{modified_sequence}\t{phi_value}\n")
                elif middle_nucleotide == "C":
                    for i, nucleotide in enumerate(["A","G","T"]):        
                        modified_sequence = f"{SEQUENCE[0:middle_index]}[{middle_nucleotide}>{nucleotide}]{SEQUENCE[middle_index+1:]}"
                        phi_value = phi_values[i] if i < len(phi_values) else None
                        if nucleotide == "T" and SEQUENCE[middle_index+1] == "G":
                            EVEN_CpG_to_T_rates = np.append(EVEN_CpG_to_T_rates, phi_value)
                        elif nucleotide =="T" and SEQUENCE[middle_index+1] != "G" and SEQUENCE[middle_index+2] == "G":
                            EVEN_CHG_to_T_rates = np.append(EVEN_CHG_to_T_rates, phi_value)
                        elif nucleotide == "T" and SEQUENCE[middle_index+1] != "G" and SEQUENCE[middle_index+2] != "G":
                            EVEN_C_to_T_rates = np.append(EVEN_C_to_T_rates, phi_value)
                        EVEN_fittedRates = np.append(EVEN_fittedRates, phi_value)
                        output.write(f"{modified_sequence}\t{phi_value}\n")
        except KeyError:
            print(f"Layer {LAYER} not found or no sequences in the JSON file.")

    return EVEN_fittedRates, EVEN_CpG_to_T_rates, EVEN_CHG_to_T_rates, EVEN_C_to_T_rates

def extract_ODD_prob_values(SPECIES, LAYER):
    json_file = os.path.join(f"/project/voight_subrate2/fabramos/100_Species/2024_BAYMER_OUT_POST_SUMMARIES/{SPECIES}.non_genic.ODD.final_posterior_summaries.json")
    output_file = f"/project/voight_subrate2/fabramos/100_Species/Baymer_Layer_probabilities_inCosmicAnalysisFormat/{SPECIES}_Layer_{LAYER}_ODD_BaymerProbabilities.txt"
    ODD_fittedRates = []
    ODD_CpG_to_T_rates = []
    ODD_C_to_T_rates = []
    ODD_CHG_to_T_rates = []
    # Load the JSON data from the file
    with open(output_file, 'w') as output:
        output.write("Type\tProbability_Value\n")
        with open(json_file, 'r') as file:
            data = json.load(file)

        try:
            # Extract all sequences in the specified layer
            sequences = data["ODD"][LAYER]["p_vec"].keys()

            # Iterate through all sequences and extract phi values
            for SEQUENCE in sequences:
                phi_values = data["ODD"][LAYER]["p_vec"][SEQUENCE]["mean"]

                if len(SEQUENCE) % 2 == 0:
                    middle_index = (len(SEQUENCE)//2) - 1
                else:
                    middle_index = len(SEQUENCE)//2

                middle_nucleotide = SEQUENCE[middle_index]
                if middle_nucleotide == "A":
                    for i, nucleotide in enumerate(["C","G","T"]):
                        modified_sequence = f"{SEQUENCE[0:middle_index]}[{middle_nucleotide}>{nucleotide}]{SEQUENCE[middle_index+1:]}"
                        phi_value = phi_values[i] if i < len(phi_values) else None
                        ODD_fittedRates = np.append(ODD_fittedRates, phi_value)
                        output.write(f"{modified_sequence}\t{phi_value}\n")
                elif middle_nucleotide == "C":
                    for i, nucleotide in enumerate(["A","G","T"]):
                        modified_sequence = f"{SEQUENCE[0:middle_index]}[{middle_nucleotide}>{nucleotide}]{SEQUENCE[middle_index+1:]}"
                        phi_value = phi_values[i] if i < len(phi_values) else None
                        if nucleotide == "T" and SEQUENCE[middle_index+1] == "G":
                            ODD_CpG_to_T_rates = np.append(ODD_CpG_to_T_rates, phi_value)
                        elif nucleotide =="T" and SEQUENCE[middle_index+1] != "G" and SEQUENCE[middle_index+2] == "G":
                            ODD_CHG_to_T_rates = np.append(ODD_CHG_to_T_rates, phi_value)
                        elif nucleotide == "T" and SEQUENCE[middle_index+1] != "G" and SEQUENCE[middle_index+2] != "G":
                            ODD_C_to_T_rates = np.append(ODD_C_to_T_rates, phi_value)
                        ODD_fittedRates = np.append(ODD_fittedRates, phi_value)
                        output.write(f"{modified_sequence}\t{phi_value}\n")
        except KeyError:
            print(f"Layer {LAYER} not found or no sequences in the JSON file.")

    return ODD_fittedRates, ODD_CpG_to_T_rates, ODD_CHG_to_T_rates, ODD_C_to_T_rates

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
        #LAYER = str(options_dict['-L'])
        LAYER = "4"

    except KeyError as E:
        missing_argument = E.args[0]
        print(f"Error: The required argument '{missing_argument}' is missing.")
        help()
    EVEN_fittedRates = []
    EVEN_CpG_to_T_rates = []
    EVEN_C_to_T_rates = []
    EVEN_CHG_to_T_rates = []
    ODD_fittedRates = []
    ODD_CpG_to_T_rates = []
    ODD_C_to_T_rates = []
    ODD_CHG_to_T_rates = []

    EVEN_fittedRates, EVEN_CpG_to_T_rates, EVEN_CHG_to_T_rates, EVEN_C_to_T_rates = extract_EVEN_prob_values(SPECIES=SPECIES, LAYER=LAYER)
    ODD_fittedRates, ODD_CpG_to_T_rates, ODD_CHG_to_T_rates, ODD_C_to_T_rates = extract_ODD_prob_values(SPECIES=SPECIES, LAYER=LAYER)

    orthogonal_error = 0
    for i, prob in enumerate(EVEN_fittedRates):
        orthogonal_dist = float((EVEN_fittedRates[i] - ODD_fittedRates[i]) / np.sqrt(2))
        orthogonal_error += orthogonal_dist ** 2
    RMSE = np.sqrt(orthogonal_error/len(EVEN_fittedRates))

    CpG_orthogonal_error = 0
    for i, prob in enumerate(EVEN_CpG_to_T_rates):
        orthogonal_dist = float((EVEN_CpG_to_T_rates[i] - ODD_CpG_to_T_rates[i]) / np.sqrt(2))
        CpG_orthogonal_error += orthogonal_dist ** 2
    CpG_RMSE = np.sqrt(CpG_orthogonal_error/len(EVEN_CpG_to_T_rates))

    CHG_orthogonal_error = 0
    for i, prob in enumerate(EVEN_CHG_to_T_rates):
        orthogonal_dist = float((EVEN_CHG_to_T_rates[i] - ODD_CHG_to_T_rates[i]) / np.sqrt(2))
        CHG_orthogonal_error += orthogonal_dist ** 2
    CHG_RMSE = np.sqrt(CHG_orthogonal_error/len(EVEN_CHG_to_T_rates))

    C_orthogonal_error = 0
    for i, prob in enumerate(EVEN_C_to_T_rates):
        orthogonal_dist = float((EVEN_C_to_T_rates[i] - ODD_C_to_T_rates[i]) / np.sqrt(2))
        C_orthogonal_error += orthogonal_dist ** 2
    C_RMSE = np.sqrt(C_orthogonal_error/len(EVEN_C_to_T_rates))

    #print(SPECIES,"\t",RMSE)
    print(SPECIES,"\t",f"{RMSE:.10f}")
    #print(CpG_RMSE)
    #print(CHG_RMSE)
    #print(C_RMSE)

if __name__ == "__main__":
    main(sys.argv[1:])
