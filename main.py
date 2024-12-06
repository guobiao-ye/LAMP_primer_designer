import argparse
import heuristic_primer_design
import specificity

def main(input_file):
    heuristic_primer_design.find_primers(input_file)
    specificity.specificity_screening()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find primers for a given FASTA file.')
    parser.add_argument('input_file', type=str, help='Path to the input FASTA file.')
    
    args = parser.parse_args()
    
    main(args.input_file)