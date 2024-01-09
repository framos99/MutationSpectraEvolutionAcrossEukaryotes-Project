#!/usr/bin/Python
import subprocess
import os
import gzip
import shutil
import tempfile
import pandas as pd
import sys
import getopt

def help(error_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -P => <string> Species Name (should be the same name used for working directory) REQUIRED
    -A => <FASTA> Assembly FASTA (expected to be in 'Assembly' subdirectory of the working directory) REQUIRED
    -V => <VCF> VCF file (expected to be in 'VCFs' subdirectory of the working directory) REQUIRED
    -R => <RMSK> RepeatMasker file (expected to be in 'Regions_To_Remove' subdirectory of the working directory) REQUIRED
    -F => <GTF/GFF> CDS feature file (expected to be in 'Regions_To_Remove' subdirectory of the working directory) REQUIRED
    -W => <string> Path to working directory REQUIRED
    --ploidy => <string> Sample ploidy; options are haploid or diploid (default: diploid) Currently unused by script!
ASSUMPTIONS
    * VCF file is bgzipped
    * FASTA, GTF/GFF, and RMSK files are gzipped
    * working directory has the same name as the given species name (-P)
    * 5 subdirectories are found in the working directory: Assembly, VCFs, Regions_To_Remove, Baymer_Out, Files_for_Baymer)
""")
    sys.exit(error_num)
#check if necessary file exists to avoid unnecessary processing
def run_command(command):
    subprocess.run(command, shell=True, check=True)

#remove gz extension from file name
def remove_gz_extension(file_path):
    base_name, extension = os.path.splitext(file_path)
    if extension.lower() == ".gz":
        return base_name
    else:
        return file_path
#if necessary file does not exist, execute generating commands  
def check_and_run_command(command, output_file):
    if not os.path.isfile(output_file):
        run_command(command)


#this function is invoked when the produced non-genic vcf does not contain AC and AN INFO tags (necessary for mutation counter baymer script):
def calculate_ac_an(input_vcf, output_vcf, ploidy):
    # Read the input gzipped VCF file and open the output VCF file for writing
    with gzip.open(input_vcf, 'rt') as vcf_in, open(output_vcf, 'w') as vcf_out:
        # Initialize flag to update the VCF header
        update_header = True

        for line in vcf_in:
            if line.startswith("#"):
                # This is a header line
                vcf_out.write(line)
            else:
                # This is a variant data line
                fields = line.strip().split("\t")
                info = fields[7].split(";")
                
                # Initialize AC and AN counters
                AC = 0
                AN = 0
                
                # Get the genotype columns starting from index 9
                genotypes = fields[9:]
                for gt_info in genotypes:
                    if not "./." in gt_info or ".|." in gt_info:
                        # Split the genotype information on the colon to isolate the GT field
                        gt_field = gt_info.split(":")[0]
                        # Split the GT part on '|' or '/' to count alleles
                        alleles = gt_field.split("|") if "|" in gt_field else gt_field.split("/")
                        if "1" in alleles:
                            AC += alleles.count("1")  # Count alternate alleles
                        AN += len(alleles)
                
                # Update the INFO field with AC and AN
                if "." in info:
                    info.remove(".")
                    info.insert(0, f"AC={AC}")
                    info.insert(1, f"AN={AN}")
                else:
                    info.insert(0, f"AC={AC}")
                    info.insert(1, f"AN={AN}")
                
                # Reconstruct the INFO field
                fields[7] = ";".join(info)
                
                # Write the updated variant line to the output VCF file
                vcf_out.write("\t".join(fields) + "\n")
    #update vcf header:
    subprocess.run(['sed', '-i', '2i ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">', output_vcf])
    subprocess.run(['sed', '-i', '3i ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in genotypes">', output_vcf])
    subprocess.run(['bgzip', output_vcf])

def check_and_update_ac_an_tags(input_vcf, output_vcf, ploidy, working_directory):
    # Check if AC and AN tags are already present in the VCF
    with gzip.open(input_vcf, 'rt') as vcf_in:
        for line in vcf_in:
            if line.startswith("##INFO") and "ID=AC" in line:
                print("AC and AN tags are already present in the VCF.")
                return
        else:
            temp_dir = os.path.join(working_directory, f"VCFs/temporary_work")
            os.mkdir(temp_dir)
            temp_vcf = os.path.join(temp_dir, f"temp.vcf")
            try:
                calculate_ac_an(input_vcf, temp_vcf, ploidy)

                # Move the temporary VCF to the final output location
                final_temp_vcf = os.path.join(f"{temp_vcf}.gz")
                shutil.move(final_temp_vcf, output_vcf)
            except Exception as e:
                print(f"Error: {e}")
            finally:
                # Clean up temporary directory
                shutil.rmtree(temp_dir)

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "P:A:V:R:F:W:L:", ['ploidy'])
    except getopt.GetoptError:
        print("Error: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)

    #check if file is provided to read in arguments
    if '-L' in options_dict:
        input_file = str(options_dict['-L'])
        print(f"Reading arguments from file: {input_file}")
        with open(input_file, 'r') as ArgFile:
            lines = ArgFile.readlines()
            #take rows from file as arguments
            file_args = [arg.strip() for line in lines for arg in line.split()]

        #combine file-extracted arguments with command-line arguments when file is incomplete
        command_line_args = sys.argv[1:] + file_args
        opts, args = getopt.getopt(command_line_args, "P:A:V:R:F:W:L:", ['ploidy'])

        options_dict = dict(opts)

    ## Required arguments
    try:
        PREFIX = str(options_dict['-P'])
        ASSEMBLY = str(options_dict['-A'])
        VCF = str(options_dict['-V'])
        RMSK = str(options_dict['-R'])
        CDS = str(options_dict['-F'])
        working_directory = str(options_dict['-W'])
        PLOIDY = options_dict.get('--ploidy', "diploid")

    except KeyError as E:
        missing_argument = E.args[0]
        print(f"Error: The required argument '{missing_argument}' is missing.")
        help()
    print("Acceptable Inputs Given")

    #Input file paths:
    assembly_file = os.path.join(working_directory, f"Assembly/{ASSEMBLY}")
    source_file_for_rmsk_bedfile = os.path.join(working_directory, f"Regions_To_Remove/{RMSK}")
    source_file_for_exon_bedfile = os.path.join(working_directory, f"Regions_To_Remove/{CDS}")
    vcf_file = os.path.join(working_directory, f"VCFs/{VCF}")

    #intermediate and output file paths
    rmsk_bedfile = os.path.join(working_directory, f"Regions_To_Remove/{PREFIX}_rmsk.bed")
    exons_bedfile = os.path.join(working_directory, f"Regions_To_Remove/{PREFIX}_exons.bed")
    callability_mask = os.path.join(working_directory, f"Regions_To_Remove/{PREFIX}_callabilityMask.bed")
    merged_bedfile = os.path.join(working_directory, f"Regions_To_Remove/{PREFIX}_exons.and.rmsk.and.lowcomplex.bed")
    whole_genome_bedfile = os.path.join(working_directory, f"Assembly/{PREFIX}_wholeGenome.bed")
    non_genic_bedfile = os.path.join(working_directory, f"Files_for_Baymer/{PREFIX}_nonGenic.bed")
    non_genic_fasta = os.path.join(working_directory, f"Files_for_Baymer/{PREFIX}_nonGenic.fa")
    uncompressed_assembly_file = remove_gz_extension(assembly_file)
    non_genic_vcf = os.path.join(working_directory, f"Files_for_Baymer/non.genic_{VCF}")
    temp_bedfile = os.path.join(working_directory, f"Regions_To_Remove/temp_combined.bed")
    assembly_index = os.path.join(working_directory, f"Assembly/{PREFIX}.fai")
    non_genic_fasta_index = os.path.join(working_directory, f"Files_for_Baymer/{PREFIX}_nonGenic.fa.fai")
    vcf_file_index = os.path.join(working_directory, f"VCFs/{VCF}.tbi")

    # Command-line commands
    gzipped_rmsk_bed_generating_command = ("gunzip -c {source_file_for_rmsk_bedfile} | sed '/^$/d' | awk 'BEGIN{{OFS=\"\\t\"}}; {{print $5,$6,$7}}' | grep -v -e 'query' -e 'sequence' > {rmsk_bedfile}").format(source_file_for_rmsk_bedfile=source_file_for_rmsk_bedfile, rmsk_bedfile=rmsk_bedfile)
    
    rmsk_bed_generating_command = ("sed '/^$/d' {source_file_for_rmsk_bedfile} | awk 'BEGIN{{OFS=\"\\t\"}}; {{print $5,$6,$7}}' | grep -v -e 'query' -e 'sequence' > {rmsk_bedfile}").format(source_file_for_rmsk_bedfile=source_file_for_rmsk_bedfile, rmsk_bedfile=rmsk_bedfile)
    
    bedfile_merging_command = ("sort -k1,1V -k2,2n {bed_to_merge} | bedtools merge -i stdin > {merged_bed}").format(bed_to_merge=temp_bedfile, merged_bed=merged_bedfile)
    
    assembly_indexing_command = ("samtools faidx {uncompressed_assembly_file} -o {assembly_index}").format(uncompressed_assembly_file=uncompressed_assembly_file, assembly_index=assembly_index)
    
    generating_wholegenome_bed_command = ("awk 'BEGIN{{OFS=\"\t\"}}; {{print $1,\"0\",$2}}' {assembly_index} | sort -k1,1V > {whole_genome_bedfile}").format(assembly_index=assembly_index, whole_genome_bedfile=whole_genome_bedfile)
    
    making_nongenic_bed_command = ("bedtools subtract -a {whole_genome_bedfile} -b {merged_bedfile} > {non_genic_bedfile}").format(whole_genome_bedfile=whole_genome_bedfile, merged_bedfile=merged_bedfile, non_genic_bedfile=non_genic_bedfile)

    uncompressing_assembly_command = ("gunzip {assembly_file}").format(assembly_file=assembly_file)

    making_nongenic_fasta_command = ("bedtools getfasta -fi {uncompressed_assembly_file} -bed {non_genic_bedfile} -fo {non_genic_fasta}").format(uncompressed_assembly_file=uncompressed_assembly_file, non_genic_bedfile=non_genic_bedfile, non_genic_fasta=non_genic_fasta)

    nongenic_assembly_indexing_command = ("samtools faidx {non_genic_fasta} -o {non_genic_fasta_index}").format(non_genic_fasta=non_genic_fasta, non_genic_fasta_index=non_genic_fasta_index)

    generate_vcf_index_command = ("bcftools index {vcf_file}").format(vcf_file=vcf_file)

    generating_nongenic_vcf_command = ("bcftools view -O z -R {non_genic_bedfile} -v snps {vcf_file} | bcftools annotate -O z -x FORMAT,^INFO/AC,^INFO/AF,^INFO/AN - | bcftools norm -O z --multiallelics - --fasta-ref {uncompressed_assembly_file} - | bcftools view -O z -o {non_genic_vcf} -v snps -e 'ALT[0]==\"*\" || ALT[0]==\".\"' -").format(non_genic_bedfile=non_genic_bedfile, vcf_file=vcf_file, uncompressed_assembly_file=uncompressed_assembly_file, non_genic_vcf=non_genic_vcf)

    compress_assembly_command = ("gzip {uncompressed_assembly_file}").format(uncompressed_assembly_file=uncompressed_assembly_file)

    dustmasker_command = ("dustmasker -in {uncompressed_assembly_file} -window 32 -outfmt acclist | sed 's/>//g' | awk 'BEGIN{{OFS=\"\t\"}}; {{print $1,$(NF-1),$NF}}' > {callability_mask}").format(uncompressed_assembly_file=uncompressed_assembly_file, callability_mask=callability_mask)

    #check if assembly is compressed, if so uncompress for processing
    if assembly_file.endswith('.gz'):
        subprocess.run(uncompressing_assembly_command, shell=True, check=True)
    # Check for RepeatMasker bedfile
    if not os.path.isfile(rmsk_bedfile):
        if source_file_for_rmsk_bedfile.endswith('.gz'):
            subprocess.run(gzipped_rmsk_bed_generating_command, shell=True, check=True)
        else:
            subprocess.run(rmsk_bed_generating_command, shell=True, check=True)

    # Check for coding sequence bedfile
    if not os.path.isfile(exons_bedfile):
        if source_file_for_exon_bedfile.endswith('.gz'):
            with gzip.open(source_file_for_exon_bedfile, 'rt') as f:
                exon_df = pd.read_csv(f, sep='\t', comment='#', header=None)

        else:
            exon_df = pd.read_csv(source_file_for_exon_bedfile, sep='\t', comment='#', header=None)
        # Filter rows based on conditions
        exon_df = exon_df[(exon_df[2].str.contains("exon"))]
        exon_df = exon_df[(~exon_df[8].str.contains("pseudogene"))]
        #save to exon bedfile
        exon_df[[0,3,4]].sort_values(by=[0, 3]).to_csv(exons_bedfile, sep='\t', index=False, header=False)

    # Check for low complexity regions mask; DustMasker
    if not os.path.isfile(callability_mask):
        subprocess.run(dustmasker_command, shell=True, check=True)


    # Merge exons, rmsk, and low complexity bedfiles
    if not os.path.isfile(merged_bedfile):
        repetitive_data = pd.read_csv(rmsk_bedfile, sep='\t', header=None, names=["chrom", "start", "end"])
        exons_data = pd.read_csv(exons_bedfile, sep='\t', header=None, names=["chrom", "start", "end"])
        lowcomplx_data = pd.read_csv(callability_mask, sep='\t', header=None, names=["chrom", "start", "end"])

        combined_rmsk_exon_lowcomplx_df = repetitive_data.append(exons_data, ignore_index=True)
        combined_rmsk_exon_lowcomplx_df = combined_rmsk_exon_lowcomplx_df.append(lowcomplx_data, ignore_index=True)
        #concatenated_df = pd.concat([repetitive_data, exons_data], axis=0, ignore_index=True)
        combined_rmsk_exon_lowcomplx_df = combined_rmsk_exon_lowcomplx_df.sort_values(by=['chrom', 'start'])

        combined_rmsk_exon_lowcomplx_df.to_csv(temp_bedfile, sep='\t', header=False, index=False)
        #run bedtools merge
        subprocess.run(bedfile_merging_command, shell=True, check=True)

        #remove temporary bedfile
        os.remove(temp_bedfile)

    # Whole genome bedfile
    if not os.path.isfile(whole_genome_bedfile):
        subprocess.run(assembly_indexing_command, shell=True, check=True)
        
        subprocess.run(generating_wholegenome_bed_command, shell=True, check=True)

    # Non-genic bedfile
    if not os.path.isfile(non_genic_bedfile):
        subprocess.run(making_nongenic_bed_command, shell=True, check=True)

    # Non-genic FASTA file
    if not os.path.isfile(non_genic_fasta):
        #run bedtools getfasta
        subprocess.run(making_nongenic_fasta_command, shell=True, check=True)

    # Check non-genic fasta index
    if not os.path.isfile(non_genic_fasta_index):
        subprocess.run(nongenic_assembly_indexing_command, shell=True, check=True)

    # Non-genic VCF file
    if not os.path.isfile(non_genic_vcf):
        #check that vcf index exists
        if not os.path.isfile(vcf_file_index):
            subprocess.run(generate_vcf_index_command, shell=True, check=True)

        #temp_dir = os.path.join(working_directory, f"VCFs/temporary_work")
        #os.mkdir(temp_dir)
        #temp_vcf = os.path.join(temp_dir, f"temp.vcf.gz")
        #fixing_vcf_header_command = ("gatk FixVcfHeader -I {vcf_file} -O {temp_vcf}").format(vcf_file=vcf_file, temp_vcf=temp_vcf)
        
        subprocess.run(generating_nongenic_vcf_command, shell=True, check=True)
        check_and_update_ac_an_tags(input_vcf=non_genic_vcf, output_vcf=non_genic_vcf, ploidy=PLOIDY, working_directory=working_directory)

    #re-compress assembly
    if os.path.isfile(non_genic_vcf) and os.path.isfile(non_genic_fasta):
        subprocess.run(compress_assembly_command, shell=True, check=True)

    # Make baymer config files
    baymer_config_file = os.path.join(working_directory, f"{PREFIX}_baymer_config.yaml")
    check_and_run_command(f"sed 's/SPECIES_NAME/{PREFIX}/g' /project/voight_subrate2/fabramos/100_Species/General_baymer_config.yaml > {baymer_config_file}", baymer_config_file)

    # Make context and mutation config file
    context_and_mutation_counts_file = os.path.join(working_directory, f"{PREFIX}_context_and_mutation_counts.yaml")
    check_and_run_command(f"sed 's/SPECIES_NAME/{PREFIX}/g' /project/voight_subrate2/fabramos/100_Species/General_context_and_mutation_counts.yaml | sed 's/FINAL_VCF/non.genic_{VCF}/g' | sed 's/NON_GENIC_BED/{PREFIX}_nonGenic.bed/g' | sed 's/NON_GENIC_FASTA/{PREFIX}_nonGenic.fa/g' > {context_and_mutation_counts_file}", context_and_mutation_counts_file)

    # Make baymer job scripts
    baymer_scripts_file = os.path.join(working_directory, f"{PREFIX}_baymer_scripts.bat")
    check_and_run_command(f"sed 's/SPECIES_NAME/{PREFIX}/g' /project/voight_subrate2/fabramos/100_Species/General_all_baymer_scripts.bat > {baymer_scripts_file}", baymer_scripts_file)

if __name__ == "__main__":
    main(sys.argv[1:])