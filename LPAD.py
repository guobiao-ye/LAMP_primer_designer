import argparse
from heuristic_primer_design import find_primers
from specificity import specificity_screening
from batch_process import process_input_file, save_full_scores
from final_output import merge_results
 
DEFAULT_INPUT_FILE = './data/input/example.fasta'
DEFAULT_REF_FILE = '"./data/resource/hg38.fa"'
DEFAULT_OUTPUT_DIR = './data/output/Final_score'

def main(input_file, ref_file, output_dir):
    find_primers(input_file)
    specificity_screening(input_file, ref_file)
    process_input_file()
    save_full_scores()
    merge_results(output_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find primers for a given FASTA file.')

    parser.add_argument('-i', '--input_file', type=str, default=DEFAULT_INPUT_FILE, help='Path to the input FASTA file (default: %(default)s).')
    parser.add_argument('-r', '--ref_file', type=str, default=DEFAULT_REF_FILE, help='Path to the reference FASTA file (default: %(default)s).')
    parser.add_argument('-o', '--output_dir', type=str, default=DEFAULT_OUTPUT_DIR, help='Directory where output will be saved (default: %(default)s).')

    args = parser.parse_args()
    
    main(args.input_file, args.ref_file, args.output_dir)
