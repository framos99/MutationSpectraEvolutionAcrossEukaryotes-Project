#!/bin/bash

cd /project/voight_subrate2/fabramos/100_Species/SPECIES_NAME/

python /project/voight_subrate2/fabramos/Baymer/baymer/context_counter.py -c SPECIES_NAME_context_and_mutation_counts.yaml --feature non_genic -m 9 --co ./Baymer_Out/context_counts.tsv -h
python /project/voight_subrate2/fabramos/Baymer/baymer/testingname.mutation_counter.py -c SPECIES_NAME_context_and_mutation_counts.yaml --feature non_genic -m 9 -p SPECIES_NAME --mo ./Baymer_Out/mutation_counts.tsv -h --ac 1 --min --ac-syntax "AC" --an-syntax "AN"

python /project/voight_subrate2/adamschr/git_repositories/Mutation_Rate_Modeling/scripts/species_variant_counting/species_baymer_wrapper.v2.py -c SPECIES_NAME_baymer_config.yaml 
