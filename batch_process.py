import csv
from calc_features_single import calc_tm, gc_percent, calc_terminal_delta_g, if_secondary_structure
from calc_features_multi import primers_distance, if_dimer
from primer_score import score_primers  # Ensure the scoring function is correctly imported
from Bio.Seq import Seq


def reverse_complement(sequence):
    """
    Get the reverse complement of a DNA sequence.
    :param sequence: Original DNA sequence
    :return: Reverse complement sequence
    """
    return str(Seq(sequence).reverse_complement())


def parse_primer_info(primer_info):
    """
    Parse primer information from a string into a dictionary.
    :param primer_info: Primer information in string format
    :return: Primer information as a dictionary
    """
    return eval(primer_info)  # Ensure the input format is correct


def process_input_file(input_file='data/output/Intermediate_file/specific_primer.csv',
                       output_file='data/output/Intermediate_file/score_into_model.csv'):
    """
    Process the input file, calculate scores for each primer group, and save the results.
    :param input_file: Input CSV file path
    :param output_file: Output result file path
    """
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)

        # Define the fields to keep
        primer_scores = ['gc_score', 'tm_score', 'delta_g_score', 'hairpin_score']
        primers_names = ['F3', 'F2', 'F1', 'B1', 'B2', 'B3']
        fieldnames = (
                ['lamp_id'] +
                [f"{primer}_{score}" for primer in primers_names for score in primer_scores] +
                ['distance_score', 'dimer_score']
        )

        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for idx, row in enumerate(reader):
            # Extract primer information and adjust the order
            primers = [
                parse_primer_info(row['forward_outer_info'])['sequence'],  # F3
                parse_primer_info(row['forward_middle_info'])['sequence'],  # F2
                parse_primer_info(row['forward_inner_info'])['sequence'],  # F1
                parse_primer_info(row['reverse_inner_info'])['sequence'],  # B1
                parse_primer_info(row['reverse_middle_info'])['sequence'],  # B2
                parse_primer_info(row['reverse_outer_info'])['sequence'],  # B3
            ]

            # Extract starting positions
            start_positions = [
                parse_primer_info(row['forward_outer_info'])['position'],  # F3
                parse_primer_info(row['forward_middle_info'])['position'],  # F2
                parse_primer_info(row['forward_inner_info'])['position'],  # F1
                parse_primer_info(row['reverse_inner_info'])['position'] -
                parse_primer_info(row['reverse_inner_info'])['length'],  # B1
                parse_primer_info(row['reverse_middle_info'])['position'] -
                parse_primer_info(row['reverse_middle_info'])['length'],  # B2
                parse_primer_info(row['reverse_outer_info'])['position'] -
                parse_primer_info(row['reverse_outer_info'])['length'],  # B3
            ]

            # Call the scoring function
            scores = score_primers(primers, start_positions)
            single_primer_scores = scores['single_primer_scores']

            # Construct output dictionary, keeping only the necessary columns
            output_row = {'lamp_id': f'Group_{idx + 1}'}
            for primer, details in zip(primers_names, single_primer_scores):
                for score_key in primer_scores:
                    output_row[f"{primer}_{score_key}"] = details.get(score_key, None)

            # Add distance_score and dimer_score
            output_row['distance_score'] = scores['distance_score']
            output_row['dimer_score'] = scores['dimer_score']

            # Write the result
            writer.writerow(output_row)


def save_full_scores(input_file='data/output/Intermediate_file/specific_primer.csv',
                     full_output_file='data/output/Intermediate_file/full_score_results.csv'):
    """
    Save the complete output of score_primers to a CSV file.
    :param input_file: Input CSV file path
    :param full_output_file: Output file path for complete scores
    """
    with open(input_file, 'r') as infile, open(full_output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)

        # Define the fields for full results
        fieldnames = ['lamp_id', 'full_scores']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for idx, row in enumerate(reader):
            # Extract primer information and adjust the order
            primers = [
                parse_primer_info(row['forward_outer_info'])['sequence'],  # F3
                parse_primer_info(row['forward_middle_info'])['sequence'],  # F2
                parse_primer_info(row['forward_inner_info'])['sequence'],  # F1
                parse_primer_info(row['reverse_inner_info'])['sequence'],  # B1
                parse_primer_info(row['reverse_middle_info'])['sequence'],  # B2
                parse_primer_info(row['reverse_outer_info'])['sequence'],  # B3
            ]

            # Extract starting positions
            start_positions = [
                parse_primer_info(row['forward_outer_info'])['position'],  # F3
                parse_primer_info(row['forward_middle_info'])['position'],  # F2
                parse_primer_info(row['forward_inner_info'])['position'],  # F1
                parse_primer_info(row['reverse_inner_info'])['position'] -
                parse_primer_info(row['reverse_inner_info'])['length'],  # B1
                parse_primer_info(row['reverse_middle_info'])['position'] -
                parse_primer_info(row['reverse_middle_info'])['length'],  # B2
                parse_primer_info(row['reverse_outer_info'])['position'] -
                parse_primer_info(row['reverse_outer_info'])['length'],  # B3
            ]

            # Call the scoring function
            scores = score_primers(primers, start_positions)

            # Write the full result
            writer.writerow({'lamp_id': f'Group_{idx + 1}', 'full_scores': scores})


if __name__ == "__main__":
    input_csv = 'data/output/Intermediate_file/specific_primer.csv'
    full_output_csv = 'data/output/Intermediate_file/full_score_results.csv'  # Replace with actual full score output path
    save_full_scores(input_csv, full_output_csv)

if __name__ == "__main__":
    input_csv = 'data/output/Intermediate_file/specific_primer.csv'
    output_csv = 'data/output/Intermediate_file/score_into_model.csv'  # Replace with actual output file path
    process_input_file(input_csv, output_csv)
