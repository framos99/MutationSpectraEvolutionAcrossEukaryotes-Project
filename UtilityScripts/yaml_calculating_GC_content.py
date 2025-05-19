from Bio import SeqIO
import sys
import yaml

def calculate_gc_content(fasta_file):
    fasta_gc_count = 0
    fasta_nucleotide_count = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        gc_count = sequence.count('G') + sequence.count('C') + sequence.count('g') + sequence.count('c')
        fasta_gc_count += gc_count
        fasta_nucleotide_count += sequence.count('G') + sequence.count('C') + sequence.count('T') + sequence.count('A') + sequence.count('g') + sequence.count('c') + sequence.count('t') + sequence.count('a')

    return fasta_gc_count, fasta_nucleotide_count

def main(yaml_file):
    total_gc_count = 0
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
        fasta_gc_count = 0
        fasta_nucleotide_count = 0

        fasta_gc_count, fasta_nucleotide_count = calculate_gc_content(fasta_file)
        total_gc_count += fasta_gc_count
        total_nucleotide_count += fasta_nucleotide_count
    gc_content = total_gc_count / total_nucleotide_count
    print(gc_content)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <YAML>")
        sys.exit(1)

    yaml_file = sys.argv[1]

    main(yaml_file)
