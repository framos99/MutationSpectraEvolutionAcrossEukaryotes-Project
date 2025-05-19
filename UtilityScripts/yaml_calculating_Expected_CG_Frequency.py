from Bio import SeqIO
import sys
import yaml

def calculate_OE_CpG_content(fasta_file):
    fasta_C_count = 0
    fasta_G_count = 0
    fasta_CG_count = 0
    fasta_nucleotide_count = 0
    fasta_dinucleotide_count = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        G_count = sequence.count('G') + sequence.count('g')
        C_count = sequence.count('C') + sequence.count('c')
        CG_count = sequence.count('CG') + sequence.count('cg') + sequence.count('Cg') + sequence.count('cG')
        fasta_C_count += C_count
        fasta_G_count += G_count
        fasta_CG_count += CG_count
        fasta_nucleotide_count += sequence.count('G') + sequence.count('C') + sequence.count('T') + sequence.count('A') + sequence.count('g') + sequence.count('c') + sequence.count('t') + sequence.count('a')
        fasta_dinucleotide_count += sequence.count('GG') + sequence.count('GC') + sequence.count('GT') + sequence.count('GA') + sequence.count('CG') + sequence.count('CC') + sequence.count('CT') + sequence.count('CA') + sequence.count('TG') + sequence.count('TC') + sequence.count('TT') + sequence.count('TA') + sequence.count('AG') + sequence.count('AC') + sequence.count('AT') + sequence.count('AA') + sequence.count('gg') + sequence.count('gc') + sequence.count('gt') + sequence.count('ga') + sequence.count('cg') + sequence.count('cc') + sequence.count('ct') + sequence.count('ca') + sequence.count('tg') + sequence.count('tc') + sequence.count('tt') + sequence.count('ta') + sequence.count('ag') + sequence.count('ac') + sequence.count('at') + sequence.count('aa') + sequence.count('gG') + sequence.count('gC') + sequence.count('gT') + sequence.count('gA') + sequence.count('cG') + sequence.count('cC') + sequence.count('cT') + sequence.count('cA') + sequence.count('tG') + sequence.count('tC') + sequence.count('tT') + sequence.count('tA') + sequence.count('aG') + sequence.count('aC') + sequence.count('aT') + sequence.count('aA') + sequence.count('Gg') + sequence.count('Gc') + sequence.count('Gt') + sequence.count('Ga') + sequence.count('Cg') + sequence.count('Cc') + sequence.count('Ct') + sequence.count('Ca') + sequence.count('Tg') + sequence.count('Tc') + sequence.count('Tt') + sequence.count('Ta') + sequence.count('Ag') + sequence.count('Ac') + sequence.count('At') + sequence.count('Aa')

    return fasta_C_count, fasta_G_count, fasta_CG_count, fasta_nucleotide_count, fasta_dinucleotide_count


def main(yaml_file):
    total_C_count = 0
    total_G_count = 0
    total_CG_count = 0
    total_nucleotide_count = 0
    total_dinucleotide_count = 0

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
        fasta_CG_count = 0
        fasta_nucleotide_count = 0
        fasta_dinucleotide_count = 0
        #calculate counts for individual fastas
        fasta_C_count, fasta_G_count, fasta_CG_count, fasta_nucleotide_count, fasta_dinucleotide_count = calculate_OE_CpG_content(fasta_file)
        #sum counts across fasta files
        total_C_count += fasta_C_count
        total_G_count += fasta_G_count
        total_CG_count += fasta_CG_count
        total_nucleotide_count += fasta_nucleotide_count
        total_dinucleotide_count += fasta_dinucleotide_count
    #compute Observed/Expected CpG ratio
    Expected_CG_perc = (total_C_count / total_nucleotide_count) * (total_G_count / total_nucleotide_count)
    Observed_CG_perc = total_CG_count / total_dinucleotide_count
    OE_CpG_ratio = Observed_CG_perc / Expected_CG_perc
    print(Expected_CG_perc)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <YAML>")
        sys.exit(1)

    yaml_file = sys.argv[1]

    main(yaml_file)
