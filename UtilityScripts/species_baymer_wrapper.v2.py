#!/usr/bin/Python3

# Created by: Christopher J Adams 12/6/2022
# 

# NOTE: Written in python3

###############################################################################
###
### This script will use the mutation and context count tsvs to generate a host of baymer models
###
###############################################################################

#import cProfile
import sys
import getopt
import os
import re
import yaml
import subprocess
import time
import numpy as np

#sys.path.append(os.path.join(os.path.dirname(sys.path[0]), 'modules'))
#import alpha_generator_class

def help():
    print("""-----------------------------------------------------------------
ARGUMENTS
    -c => <file> config file REQUIRED

ASSUMPTIONS
    * output directory is empty and has no current directory structure, as this
      script will generate these for you.
    * Assumes the base parameter bed file is to be used
""")
    sys.exit()
###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################


## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "c:")
    except getopt.GetoptError:
        print("Error: Incorrect usage of getopts flags!")
        help() 
    options_dict = dict(opts)

    try:
        config_file = options_dict["-c"]
    except KeyError:
        print("Error: One of your required arguments is missing.")
        help()

    print("Acceptable Inputs Given")

    driver(config_file)



###############################################################################
#############################  DRIVER  ########################################
###############################################################################


## drive the script ##
## ONE-TIME CALL -- called by main

def driver(config_file):

    ## Load yaml
    config_dict = yaml.load(open(config_file, 'r'), Loader=yaml.SafeLoader)
    
    base_output_directory = config_dict["base_output_location"]
    feature = config_dict["feature"]
    species = config_dict["species"]
    pop = species

    count_generation_config = config_dict["count_generation_config"]
    context_counts_tsv = config_dict["context_count_tsv"]
    mutation_counts_csv = config_dict["mutation_count_tsv"]
    # create base directory structure
    subprocess.run(["mkdir", "{}/{}/".format(base_output_directory, species)])
    subprocess.run(["mkdir", "{}/{}/{}/".format(base_output_directory, species, feature)])
    bsub_directory = "{}/{}/bsub_files/".format(base_output_directory, species)
    out_and_err_directory = "{}/out_and_err_files/".format(bsub_directory)
    subprocess.run(["mkdir", bsub_directory])
    subprocess.run(["mkdir", out_and_err_directory])
    
    # gather appropriate scripts:
    combine_tsvs_script = "/project/voight_subrate2/fabramos/Baymer/baymer/generate_count_json.py" 
    run_baymer_script = "/project/voight_subrate2/fabramos/Baymer/baymer/run_baymer.py"
    plot_baymer_script = "/project/voight_subrate2/fabramos/Baymer/baymer/plot_baymer_posterior_distributions.py"
    generate_rate_dict_script = "/project/voight_subrate2/fabramos/Baymer/baymer/collapse_model/convert_mcmc_summary_json_to_count_dicts.py"

    # iterate through conditions and generate count structures
    for condition_name in list(config_dict["conditions"].keys()):
        
        baymer_output_directory = "{}/{}/{}/{}_baymer_out_files/".format(base_output_directory, species, feature, condition_name) 
        baymer_count_directory = baymer_output_directory + condition_name + "_count_jsons/"
        subprocess.run(["mkdir", "-p", baymer_count_directory])
        datasets = config_dict["conditions"][condition_name]["datasets"]
        
        max_af = config_dict["conditions"][condition_name]["filters"]["max_AF"]
        min_ac = config_dict["conditions"][condition_name]["filters"]["min_AC"]
        min_quality = config_dict["conditions"][condition_name]["filters"]["min_quality"]

        # generate the bsub files for each condition
        for dataset in datasets:

            # first create json bsubs
            create_json_bsub_file = bsub_directory + "{}.{}.{}.{}.create_jsons.sh".format(species, feature, condition_name, dataset)
            
            create_json_bsub = open(create_json_bsub_file, "w")
            subprocess.run(["chmod", "770", create_json_bsub_file])
            
            script_flags = "-c {} -f {} -p {} --mc {} --cc {} -d {} -o {}".format( \
                           count_generation_config, feature, species, mutation_counts_csv, context_counts_tsv, dataset, baymer_count_directory)
            
            # next go through the filters
            for flag, filter_amount in zip(["--max-af", "--min-ac", "--quality"], [max_af, min_ac, min_quality]):
                try:
                    int(filter_amount)
                except ValueError:
                    continue
                else:
                    script_flags += " {} {}".format(flag, filter_amount)



            create_json_bsub.write("python {} {}".format(combine_tsvs_script, script_flags))
            
            # submit json jobs
            create_json_job_name = condition_name + "_" + dataset + "_bp"
            create_json_err = out_and_err_directory + "create_json." + create_json_job_name + ".err"
            create_json_out = out_and_err_directory + "create_json." + create_json_job_name + ".out"
            
            create_jsons = subprocess.run(["bsub", "-q", "math_normal", "-eo", create_json_err, "-oo", create_json_out, create_json_bsub_file], stdout=subprocess.PIPE)
            create_jsons_job_id = get_job_id(create_jsons)
            
            json_config = "{}1_9mer.{}.{}.{}.hardcoded_count_files.yaml".format(baymer_count_directory, dataset, pop, feature)

            # write baymer config file - note this assumes default settings
            baymer_config_file =  baymer_output_directory + "baymer_config." + create_json_job_name + ".yaml"
            write_baymer_config(baymer_config_file, baymer_output_directory, pop, feature, dataset)

            
            # create baymer bsubs
            baymer_job_ids = []
            for rs_index in [0, 1]:
                run_baymer_bsub_file = bsub_directory + "{}.{}.{}.{}.rs_index_{}.run_baymer.sh".format(species, feature, condition_name, dataset, rs_index)
 
                run_baymer_bsub = open(run_baymer_bsub_file, "w")
                subprocess.run(["chmod", "770", run_baymer_bsub_file])
                
                #run_baymer_bsub.write("python {} \
                #-c {} -p {} -r {} -t 8".format(run_baymer_script, json_config, baymer_config_file, rs_index))
                run_baymer_bsub.write("python {} \
                -c {} -p {} -r {}".format(run_baymer_script, json_config, baymer_config_file, rs_index))
                
                run_baymer_job_name = create_json_job_name + "_rs_index_" + str(rs_index)
                run_baymer_err = out_and_err_directory + "run_baymer." + run_baymer_job_name + ".err"
                run_baymer_out = out_and_err_directory + "run_baymer." + run_baymer_job_name + ".out"
                # submit run baymer!
                #run_baymer = subprocess.run(["bsub", "-q", "math_normal", "-m", "math1", '-R "rusage[mem=25GB]"', "-M", "26G", "-w", create_jsons_job_id, "-e", run_baymer_err, "-o", run_baymer_out, run_baymer_bsub_file], stdout=subprocess.PIPE)
                run_baymer = subprocess.run(["bsub", "-q", "math_normal", "-n 8", '-R "rusage[mem=25GB]"', "-M", "26G", "-w", create_jsons_job_id, "-e", run_baymer_err, "-o", run_baymer_out, run_baymer_bsub_file], stdout=subprocess.PIPE)
                run_baymer_job_id = get_job_id(run_baymer)
                baymer_job_ids.append(run_baymer_job_id)
            
            # finally launch the plotting script
            plot_baymer_bsub_file = bsub_directory + "{}.{}.{}.{}.plot_baymer.sh".format(species, feature, condition_name, dataset)
            
            plot_baymer_bsub = open(plot_baymer_bsub_file, "w")
            subprocess.run(["chmod", "770", plot_baymer_bsub_file])
            
            plot_baymer_bsub.write("python {} -c {}".format(plot_baymer_script, baymer_config_file))

            plot_baymer_err = out_and_err_directory + "plot_baymer." + create_json_job_name + ".err"
            plot_baymer_out = out_and_err_directory + "plot_baymer." + create_json_job_name + ".out"
            
            run_baymer_job_dependency = '\"(' + " && ".join([str(x) for x in baymer_job_ids]) + ')\"'
            #print(run_baymer_job_dependency)
            subprocess_string = " ".join([str(x) for x in ["bsub", "-q", "math_normal", "-w", run_baymer_job_dependency, "-eo", plot_baymer_err, "-oo", plot_baymer_out, plot_baymer_bsub_file]]) 
            #print(subprocess_string)
            plot_baymer = subprocess.run(subprocess_string, shell = True, stdout=subprocess.PIPE)

            ## next I can perform all the other plotting I might be interested in
            plot_baymer_job_id = get_job_id(plot_baymer)
            
            mcmc_summary_dict = "{0}/{1}/{2}/{3}_baymer_out_files/{4}_outplots/{1}.{2}.{4}.final_posterior_summaries.json".format(base_output_directory, species, feature, condition_name, dataset)

            # generate rate dicts
            rate_dict_output_directory = "{}/{}/{}/{}_baymer_out_files/{}_outplots/rate_dicts/".format(base_output_directory, species, feature, condition_name, dataset)
            rate_dict_baymer_bsub_file = bsub_directory + "{}.{}.{}.{}.baymer_map_rate_dict.sh".format(species, feature, condition_name, dataset)
            

            rate_dict_baymer_bsub = open(rate_dict_baymer_bsub_file, "w")
            subprocess.run(["chmod", "770", rate_dict_baymer_bsub_file])
            rate_dict_baymer_bsub.write("mkdir " + rate_dict_output_directory + "\n")
            rate_dict_baymer_bsub.write("python {} -a {} -p {} -f {} -o {}".format(generate_rate_dict_script, mcmc_summary_dict, pop, feature, rate_dict_output_directory))

            rate_dict_baymer_err = out_and_err_directory + "baymer_map_rate_dict." + create_json_job_name + ".err"
            rate_dict_baymer_out = out_and_err_directory + "baymer_map_rate_dict." + create_json_job_name + ".out"
            
            run_baymer = subprocess.run(["bsub", "-q", "math_normal", "-w", plot_baymer_job_id, "-eo", rate_dict_baymer_err, "-oo", rate_dict_baymer_out, rate_dict_baymer_bsub_file], stdout=subprocess.PIPE)




            

def get_job_id(std_out):
    
    job_id_list = str(std_out).strip().split('<')
    job_id = job_id_list[1].split('>')[0]

    return job_id


def write_baymer_config(baymer_config_file, baymer_output_directory, pop, feature, dataset):
    
    layer_parameters = {0: {
                            "iteration": 30000,
                            "burnin": 20000,
                            "thinning_parameter": 10,
                            "alpha": 0,
                            "slab_sigma": 0,
                            "alpha_sigma": 0,
                            "sigma_sigma": 0,
                            "phi_sigma": 0.00005,
                            "indicator_sampling": "gibbs"
                            },
                        1: {
                            "iteration": 30000,
                            "burnin": 20000,
                            "thinning_parameter": 10,
                            "alpha": 0.99,
                            "slab_sigma": 1.0,
                            "alpha_sigma": "beta",
                            "sigma_sigma": 0.2,
                            "phi_sigma": 0.001,
                            "indicator_sampling": "gibbs"
                            },
                        2: {
                            "iteration": 30000,
                            "burnin": 20000,
                            "thinning_parameter": 10,
                            "alpha": 0.99,
                            "slab_sigma": 0.2,
                            "alpha_sigma": "beta",
                            "sigma_sigma": 0.02,
                            "phi_sigma": 0.005,
                            "indicator_sampling": "gibbs"
                            },

                        3: {
                            "iteration": 30000,
                            "burnin": 20000,
                            "thinning_parameter": 10,
                            "alpha": 0.99,
                            "slab_sigma": 0.18,
                            "alpha_sigma": "beta",
                            "sigma_sigma": 0.006,
                            "phi_sigma": 0.005,
                            "indicator_sampling": "gibbs"
                            },

                        4: {
                            "iteration": 40000,
                            "burnin": 30000,
                            "thinning_parameter": 10,
                            "alpha": 0.6,
                            "slab_sigma": 0.25,
                            "alpha_sigma": "beta",
                            "sigma_sigma": 0.004,
                            "phi_sigma": 0.005,
                            "indicator_sampling": "gibbs"
                            },

                        5: {
                            "iteration": 40000,
                            "burnin": 30000,
                            "thinning_parameter": 10,
                            "alpha": 0.5,
                            "slab_sigma": 0.15,
                            "alpha_sigma": "beta",
                            "sigma_sigma": 0.005,
                            "phi_sigma": 0.005,
                            "indicator_sampling": "gibbs"
                            },

                        6: {
                            "iteration": 50000,
                            "burnin": 40000,
                            "thinning_parameter": 10,
                            "alpha": 0.4,
                            "slab_sigma": 0.2,
                            "alpha_sigma": "beta",
                            "sigma_sigma": 0.005,
                            "phi_sigma": 0.005,
                            "indicator_sampling": "gibbs"
                            },

                        7: {
                            "iteration": 50000,
                            "burnin": 40000,
                            "thinning_parameter": 10,
                            "alpha": 0.2,
                            "slab_sigma": 0.2,
                            "alpha_sigma": "beta",
                            "sigma_sigma": 0.0025,
                            "phi_sigma": 0.005,
                            "indicator_sampling": "gibbs"
                            },

                        8: {
                            "iteration": 50000,
                            "burnin": 40000,
                            "thinning_parameter": 10,
                            "alpha": 0.1,
                            "slab_sigma": 0.33,
                            "alpha_sigma": "beta",
                            "sigma_sigma": 0.001,
                            "phi_sigma": 0.005,
                            "indicator_sampling": "gibbs"
                            }
                        }


    random_seeds = np.random.randint(1, 999, size=2, dtype=int)

    

    with open(baymer_config_file, "w") as config:
        config.write("---\n")
        config.write("max_mer: 9\n")
        config.write("pop: " + pop + "\n")
        config.write("feature: " + feature + "\n")
        config.write("c: 0.04472136\n")
        config.write("set_sigma: 0.8540486\n")
        config.write("posterior_dir: " + baymer_output_directory + "\n")
        config.write("dataset: " + dataset + "\n")
        config.write('alternation_pattern: right' + "\n")
        config.write("random_seeds:\n")
        for rs in random_seeds:
            config.write("  - " + str(rs) + "\n")
        for layer in layer_parameters:
            config.write(str(layer) + ":\n")
            
            num_iterations = layer_parameters[layer]["iteration"]
            burnin = layer_parameters[layer]["burnin"]
            config.write("  iteration_burnin:\n")
            config.write("    - " + str(num_iterations) + "\n")
            config.write("    - " + str(burnin) + "\n")
            
            thinning_parameter = layer_parameters[layer]["thinning_parameter"]
            config.write("  thinning_parameter: " + str(thinning_parameter) + "\n")
            
            alpha = layer_parameters[layer]["alpha"]
            config.write("  alpha: " + str(alpha) + "\n")
            
            slab_sigma = layer_parameters[layer]["slab_sigma"]
            config.write("  slab_sigma: " + str(slab_sigma) + "\n")
        
            config.write("  proposal_sigmas:\n")
            alpha_sigma = layer_parameters[layer]["alpha_sigma"]
            config.write("    alpha_sigma: " + str(alpha_sigma) + "\n")
            
            sigma_sigma = layer_parameters[layer]["sigma_sigma"]
            config.write("    sigma_sigma: " + str(sigma_sigma) + "\n")
            
            phi_sigma = layer_parameters[layer]["phi_sigma"]
            config.write("    phi_sigma: " + str(phi_sigma) + "\n")

            i_s = layer_parameters[layer]["indicator_sampling"]
            config.write("    indicator_sampling: " + str(i_s) + "\n")
        
        config.write("...")



class VCFPrepper:

    def __init__(self, output_directory, pop_class_list, config_dict, ancestral, bed, delete):
        
        self.base_dir = output_directory
        
        self.pop_class_list = pop_class_list
        
        self.config_dict = config_dict
        
        self.ancestral_bool = ancestral
        self.delete = delete

        self.bed_string = ""
        if bed:
            self.bed_string = "-R {}".format(self.config_dict['bed'][bed])
        
        self.full_pop_list = config_dict['full_pop_lists']['ALL_admixed']
        self.bsub_files_dir = "{}/bsub_files/".format(output_directory)
        self.intermediate_file_dir = "{}/intermediate_files/".format(output_directory)
        self.final_output_dir = "{}/final_output/".format(output_directory)
        
    def create_directories(self):
        
        subprocess.run(["mkdir", self.bsub_files_dir])
        subprocess.run(["mkdir", "{}/out_and_err_files/".format(self.bsub_files_dir)])

        subprocess.run(["mkdir", self.intermediate_file_dir])

        subprocess.run(["mkdir", self.final_output_dir])
        
        if self.ancestral_bool:
            subprocess.run(["mkdir", "{}/ancestral_variant_vcfs/".format(self.intermediate_file_dir)])

        for pop in self.pop_class_list:
            subprocess.run(["mkdir", "{}/{}/".format(self.final_output_dir, pop.pop_string)])
        
    def write_bcftools_filter(self, config_file):

        for chrom in self.config_dict['input_vcfs']:
            
            chrom_init_filter_name = "{}/{}.bcftools_view.vcf.gz".format(self.intermediate_file_dir, chrom)
            
            chrom_file = self.config_dict['input_vcfs'][chrom]
             
            job_name = "bcftools_filter.{}.sh".format(chrom)
            job_file = open("{}/{}".format(self.bsub_files_dir, job_name), 'w')
            job_file.write("module load bcftools\nmodule load htslib\nmodule load perl\n")
            
            job_file.write("## {}\n".format(chrom))

            ## removes multi-allelic, non-PASS, and indel sites, and also culls to just bed region
            
            job_file.write("bcftools view -O z -o {} {} -m 2 -M 2 -v snps -f .,PASS " \
                           "--force-samples --samples-file {} --min-ac 2 {}\n".format( \
                            chrom_init_filter_name, self.bed_string, self.full_pop_list, chrom_file))

            if self.ancestral_bool:
                
                job_file.write("python3 /project/voight_subrate2/fabramos/Baymer/" \
                               "baymer/ancestral_vcf_generator.py " \
                               "-c {} -i {} -o {} -n {} -v \n".format(config_file, \
                               chrom_init_filter_name, self.intermediate_file_dir, chrom[3:]))
                if self.delete:
                    job_file.write("command rm {}\n".format(chrom_init_filter_name))
                job_file.write("bgzip {}/ancestral_variant_vcfs/{}.bcftools_view.ancestral_variants.vcf\n" \
                                .format(self.intermediate_file_dir, chrom))
                job_file.write("tabix {}/ancestral_variant_vcfs/{}.bcftools_view.ancestral_variants.vcf.gz\n" \
                                .format(self.intermediate_file_dir, chrom))

                intermediate_bcftools_out = "{}/ancestral_variant_vcfs/{}.bcftools_view.ancestral_variants.vcf.gz".format(self.intermediate_file_dir, chrom)
                
                chrom_init_filter_name = "{}/{}.bcftools_view_2.ancestral.vcf.gz".format(self.intermediate_file_dir, chrom)

                job_file.write("bcftools view -O z -o {} {} -m 2 -M 2 -v snps -f .,PASS " \
                               "--force-samples --samples-file {} --min-ac 2 {}\n".format( \
                               chrom_init_filter_name, self.bed_string, self.full_pop_list, intermediate_bcftools_out))
            
            for pop in self.pop_class_list:
                bcftools_out = "{0}/{1}/{2}.{1}.vcf.gz".format(self.final_output_dir, pop.pop_string, chrom)
                pop.chrom_file_dict[chrom][0] = bcftools_out

                job_file.write("bcftools view -O z -o {} --min-ac 1 --force-samples --samples-file {} {}\n".format(
                                bcftools_out, pop.pop_list_file, chrom_init_filter_name))
                job_file.write("tabix -p vcf {}\n".format(bcftools_out))
            
            if self.delete:
                job_file.write("command rm {}\n".format(chrom_init_filter_name))
                job_file.write("command rm {}\n".format(intermediate_bcftools_out))
            
            subprocess.run(["chmod", "770", "{}/{}".format(self.bsub_files_dir, job_name)])



    def submit_jobs(self):
        
        for chrom in self.config_dict['input_vcfs']:

            ## first submit bcftools job
            bcftools_job = "{}/bcftools_filter.{}.sh".format(self.bsub_files_dir, chrom)
            bcftools_err = "{}/out_and_err_files/bcftools_filter.{}.err".format(self.bsub_files_dir, chrom)
            bcftools_out = "{}/out_and_err_files/bcftools_filter.{}.out".format(self.bsub_files_dir, chrom)

            bcftools_filter = subprocess.run(["bsub", "-q", "math_normal", "-eo", bcftools_err, "-oo", bcftools_out, bcftools_job], stdout=subprocess.PIPE)    
    
    def get_job_id(self, std_out):
        
        job_id_list = str(std_out).strip().split('<')
        job_id = job_id_list[1].split('>')[0]
        
        return job_id


class Population:
    
    def __init__(self, pop_string, config_dict, cosmo = False):
        
        self.pop_string = pop_string
        
        if cosmo:
            self.pop_list_file = None
        else:
            self.pop_list_file = config_dict['pop_lists'][self.pop_string]

        self.chrom_file_dict = self.establish_chrom_dict(config_dict)

        self.recode_info = self.get_recode_code()

    def establish_chrom_dict(self, config_dict):
        
        chrom_list = config_dict["input_vcfs"].keys()
        value_list = [["","","",""] for n in range(len(chrom_list))]

        return dict(zip(chrom_list, value_list))


    def get_recode_code(self):
        
        return "{}_AF".format(self.pop_string[0:3])



#######################################################################################################################################################
#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])
