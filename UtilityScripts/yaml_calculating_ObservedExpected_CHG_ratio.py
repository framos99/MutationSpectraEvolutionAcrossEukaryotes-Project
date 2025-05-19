from Bio import SeqIO
import sys
import yaml

def calculate_OE_CpG_content(fasta_file):
    fasta_C_count = 0
    fasta_G_count = 0
    fasta_A_count = 0
    fasta_T_count = 0
    fasta_CHG_count = 0
    fasta_nucleotide_count = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        G_count = sequence.count('G') + sequence.count('g')
        C_count = sequence.count('C') + sequence.count('c')
        A_count = sequence.count('A') + sequence.count('a')
        T_count = sequence.count('T') + sequence.count('t')

        CHG_count = sequence.count('CAG') + sequence.count('CCG') + sequence.count('CTG') + sequence.count('cag') + sequence.count('ccg') + sequence.count('ctg') + sequence.count('Cag') + sequence.count('Ccg') + sequence.count('Ctg') + sequence.count('cAg') + sequence.count('cCg') + sequence.count('cTg') + sequence.count('caG') + sequence.count('ccG') + sequence.count('ctG') + sequence.count('CAg') + sequence.count('CCg') + sequence.count('CTg') + sequence.count('cAG') + sequence.count('cCG') + sequence.count('cTG') + sequence.count('CaG') + sequence.count('CcG') + sequence.count('CtG')

        fasta_C_count += C_count
        fasta_G_count += G_count
        fasta_A_count += A_count
        fasta_T_count += T_count
        fasta_CHG_count += CHG_count
        fasta_nucleotide_count += sequence.count('G') + sequence.count('C') + sequence.count('T') + sequence.count('A') + sequence.count('g') + sequence.count('c') + sequence.count('t') + sequence.count('a')

    return fasta_C_count, fasta_G_count, fasta_A_count, fasta_T_count, fasta_CHG_count, fasta_nucleotide_count


def main(yaml_file):
    total_C_count = 0
    total_G_count = 0
    total_A_count = 0
    total_T_count = 0
    total_CHG_count = 0
    total_nucleotide_count = 0

    with open(yaml_file, 'r') as file:
        try:
            data = yaml.safe_load(file)
            fasta_files = data['features']['non_genic']['fastas']
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit(1)
        except KeyError as exc:
            print(f"Key error: {exc}")
            sys.exit(1)
            
    for sequence_name, fasta_file in fasta_files.items():
        fasta_C_count = 0
        fasta_G_count = 0
        fasta_A_count = 0
        fasta_T_count = 0

        fasta_CHG_count = 0
        fasta_nucleotide_count = 0
        #calculate counts for individual fastas
        fasta_C_count, fasta_G_count, fasta_A_count, fasta_T_count, fasta_CHG_count, fasta_nucleotide_count = calculate_OE_CpG_content(fasta_file)
        #sum counts across fasta files
        total_C_count += fasta_C_count
        total_G_count += fasta_G_count
        total_A_count += fasta_A_count
        total_T_count += fasta_T_count

        total_CHG_count += fasta_CHG_count
        total_nucleotide_count += fasta_nucleotide_count
    #compute Observed/Expected CpG ratio
    Expected_CHG_perc = ((total_C_count / total_nucleotide_count) * (total_A_count / total_nucleotide_count) * (total_G_count / total_nucleotide_count)) + ((total_C_count / total_nucleotide_count) * (total_C_count / total_nucleotide_count) * (total_G_count / total_nucleotide_count)) + ((total_C_count / total_nucleotide_count) * (total_T_count / total_nucleotide_count) * (total_G_count / total_nucleotide_count))
    Observed_CHG_perc = total_CHG_count / total_nucleotide_count
    OE_CHG_ratio = Observed_CHG_perc / Expected_CHG_perc
    print(OE_CHG_ratio)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <YAML>")
        sys.exit(1)

    yaml_file = sys.argv[1]

    main(yaml_file)
