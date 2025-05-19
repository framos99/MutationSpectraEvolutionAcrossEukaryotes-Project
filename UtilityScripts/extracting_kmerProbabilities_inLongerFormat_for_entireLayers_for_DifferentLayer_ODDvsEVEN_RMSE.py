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
    json_file2 = os.path.join(f"/project/voight_subrate2/fabramos/100_Species/2024_BAYMER_OUT_POST_SUMMARIES/{SPECIES}.non_genic.ODD.final_posterior_summaries.json")
    L8_EVEN_fittedRates = []
    L0_orthogonal_error = 0
    L1_orthogonal_error = 0
    L2_orthogonal_error = 0
    L3_orthogonal_error = 0
    L4_orthogonal_error = 0
    L5_orthogonal_error = 0
    L6_orthogonal_error = 0
    L7_orthogonal_error = 0
    L8_orthogonal_error = 0

    # Load the JSON data from the file
    with open(json_file, 'r') as file, open(json_file2, 'r') as file2:
        data = json.load(file)
        data2 = json.load(file2)
    try:
        # Extract all 9-mer sequences
        sequences = data["EVEN"][LAYER]["p_vec"].keys()
        # Iterate through all sequences and extract phi values
        for j, SEQUENCE in enumerate(sequences):
            phi_values = data["EVEN"][LAYER]["p_vec"][SEQUENCE]["mean"]
            if len(SEQUENCE) % 2 == 0:
                middle_index = (len(SEQUENCE)//2) - 1
            else:
                middle_index = len(SEQUENCE)//2
            # for each 9-mer access the folded k-mer fitted probability
            L0_SEQUENCE = f"{SEQUENCE[4]}"
            L1_SEQUENCE = f"{SEQUENCE[4:6]}"
            L2_SEQUENCE = f"{SEQUENCE[3:6]}"
            L3_SEQUENCE = f"{SEQUENCE[3:7]}"
            L4_SEQUENCE = f"{SEQUENCE[2:7]}"
            L5_SEQUENCE = f"{SEQUENCE[2:8]}"
            L6_SEQUENCE = f"{SEQUENCE[1:8]}"
            L7_SEQUENCE = f"{SEQUENCE[1:9]}"
            L8_SEQUENCE = f"{SEQUENCE[0:9]}"

            L0_PREDphi_values = data2["ODD"]["0"]["p_vec"][L0_SEQUENCE]["mean"]
            L1_PREDphi_values = data2["ODD"]["1"]["p_vec"][L1_SEQUENCE]["mean"]
            L2_PREDphi_values = data2["ODD"]["2"]["p_vec"][L2_SEQUENCE]["mean"]
            L3_PREDphi_values = data2["ODD"]["3"]["p_vec"][L3_SEQUENCE]["mean"]
            L4_PREDphi_values = data2["ODD"]["4"]["p_vec"][L4_SEQUENCE]["mean"]
            L5_PREDphi_values = data2["ODD"]["5"]["p_vec"][L5_SEQUENCE]["mean"]
            L6_PREDphi_values = data2["ODD"]["6"]["p_vec"][L6_SEQUENCE]["mean"]
            L7_PREDphi_values = data2["ODD"]["7"]["p_vec"][L7_SEQUENCE]["mean"]
            L8_PREDphi_values = data2["ODD"]["8"]["p_vec"][L8_SEQUENCE]["mean"]

            middle_nucleotide = SEQUENCE[middle_index]
            for i, nucleotide in enumerate(["C","G","T"]):            
                phi_value = phi_values[i] if i < len(phi_values) else None
                L8_EVEN_fittedRates = np.append(L8_EVEN_fittedRates, phi_value)

                # append the fitted rate for each folded layer in 9-mer space
                L0_PREDphi_value = L0_PREDphi_values[i] if i < len(L0_PREDphi_values) else None
                L0_orthogonal_dist = float((np.log(phi_value) - np.log(L0_PREDphi_value)) / np.sqrt(2))
                L0_orthogonal_error += L0_orthogonal_dist ** 2

                L1_PREDphi_value = L1_PREDphi_values[i] if i < len(L1_PREDphi_values) else None
                L1_orthogonal_dist = float((np.log(phi_value) - np.log(L1_PREDphi_value)) / np.sqrt(2))
                L1_orthogonal_error += L1_orthogonal_dist ** 2

                L2_PREDphi_value = L2_PREDphi_values[i] if i < len(L2_PREDphi_values) else None
                L2_orthogonal_dist = float((np.log(phi_value) - np.log(L2_PREDphi_value)) / np.sqrt(2))
                L2_orthogonal_error += L2_orthogonal_dist ** 2

                L3_PREDphi_value = L3_PREDphi_values[i] if i < len(L3_PREDphi_values) else None
                L3_orthogonal_dist = float((np.log(phi_value) - np.log(L3_PREDphi_value)) / np.sqrt(2))
                L3_orthogonal_error += L3_orthogonal_dist ** 2

                L4_PREDphi_value = L4_PREDphi_values[i] if i < len(L4_PREDphi_values) else None
                L4_orthogonal_dist = float((np.log(phi_value) - np.log(L4_PREDphi_value)) / np.sqrt(2))
                L4_orthogonal_error += L4_orthogonal_dist ** 2

                L5_PREDphi_value = L5_PREDphi_values[i] if i < len(L5_PREDphi_values) else None
                L5_orthogonal_dist = float((np.log(phi_value) - np.log(L5_PREDphi_value)) / np.sqrt(2))
                L5_orthogonal_error += L5_orthogonal_dist ** 2

                L6_PREDphi_value = L6_PREDphi_values[i] if i < len(L6_PREDphi_values) else None
                L6_orthogonal_dist = float((np.log(phi_value) - np.log(L6_PREDphi_value)) / np.sqrt(2))
                L6_orthogonal_error += L6_orthogonal_dist ** 2

                L7_PREDphi_value = L7_PREDphi_values[i] if i < len(L7_PREDphi_values) else None
                L7_orthogonal_dist = float((np.log(phi_value) - np.log(L7_PREDphi_value)) / np.sqrt(2))
                L7_orthogonal_error += L7_orthogonal_dist ** 2

                L8_PREDphi_value = L8_PREDphi_values[i] if i < len(L8_PREDphi_values) else None
                L8_orthogonal_dist = float((np.log(phi_value) - np.log(L8_PREDphi_value)) / np.sqrt(2))
                L8_orthogonal_error += L8_orthogonal_dist ** 2

        RMSE_L0 = np.sqrt(L0_orthogonal_error/len(L8_EVEN_fittedRates))
        RMSE_L1 = np.sqrt(L1_orthogonal_error/len(L8_EVEN_fittedRates))
        RMSE_L2 = np.sqrt(L2_orthogonal_error/len(L8_EVEN_fittedRates))
        RMSE_L3 = np.sqrt(L3_orthogonal_error/len(L8_EVEN_fittedRates))
        RMSE_L4 = np.sqrt(L4_orthogonal_error/len(L8_EVEN_fittedRates))
        RMSE_L5 = np.sqrt(L5_orthogonal_error/len(L8_EVEN_fittedRates))
        RMSE_L6 = np.sqrt(L6_orthogonal_error/len(L8_EVEN_fittedRates))
        RMSE_L7 = np.sqrt(L7_orthogonal_error/len(L8_EVEN_fittedRates))
        RMSE_L8 = np.sqrt(L8_orthogonal_error/len(L8_EVEN_fittedRates))

    except KeyError:
        print(f"Layer {LAYER} not found or no sequences in the JSON file.")

    return RMSE_L0, RMSE_L1, RMSE_L2, RMSE_L3, RMSE_L4, RMSE_L5, RMSE_L6, RMSE_L7, RMSE_L8

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "S:")
    except getopt.GetoptError:
        print("Error: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)

    ## Required arguments
    try:
        SPECIES = str(options_dict['-S'])
        LAYER = "8"

    except KeyError as E:
        missing_argument = E.args[0]
        print(f"Error: The required argument '{missing_argument}' is missing.")
        help()

    RMSE_L0 = 0
    RMSE_L1 = 0
    RMSE_L2 = 0
    RMSE_L3 = 0
    RMSE_L4 = 0
    RMSE_L5 = 0
    RMSE_L6 = 0
    RMSE_L7 = 0
    RMSE_L8 = 0

    RMSE_L0, RMSE_L1, RMSE_L2, RMSE_L3, RMSE_L4, RMSE_L5, RMSE_L6, RMSE_L7, RMSE_L8 = extract_EVEN_prob_values(SPECIES=SPECIES, LAYER=LAYER)

    print(SPECIES,"\t",f"{RMSE_L0:.10f}","\t",f"{RMSE_L1:.10f}","\t",f"{RMSE_L2:.10f}","\t",f"{RMSE_L3:.10f}","\t",f"{RMSE_L4:.10f}","\t",f"{RMSE_L5:.10f}","\t",f"{RMSE_L6:.10f}","\t",f"{RMSE_L7:.10f}","\t",f"{RMSE_L8:.10f}")

if __name__ == "__main__":
    main(sys.argv[1:])
