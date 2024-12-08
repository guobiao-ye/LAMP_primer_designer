import csv
from collections import defaultdict
from primer_score import score_primers
import re

def process_csv(input_file, output_file):
    # Initialize a dictionary to store primer data for each lamp_id
    lamp_data = defaultdict(list)

    # Read the input CSV with UTF-8 encoding to handle special characters
    with open(input_file, 'r', encoding='utf-8') as infile:
        reader = csv.DictReader(infile)

        # Iterate through the rows of the CSV
        for row in reader:
            lamp_id = row['lamp_id']
            primer = row['name']
            sequence = row['sequence']
            start_pos = extract_start_position(row['position'])
            
            # Check if the sequence contains only A, T, C, G (case insensitive)
            if not re.match("^[ATCGatcg]*$", sequence):  # Regex to check valid nucleotides
                #print(f"Skipping lamp_id {lamp_id} because the primer sequence contains invalid characters.")
                continue  # Skip this lamp_id if the sequence is invalid

            # Check if any of the necessary columns are empty
            if not lamp_id or not primer or not sequence or start_pos is None:
                #print(f"Skipping lamp_id {lamp_id} because one or more required fields are empty.")
                continue  # Skip this lamp_id if any required field is empty
            
            # Store the primer sequence and its start position
            lamp_data[lamp_id].append((primer, sequence, start_pos))
    
    # Initialize an empty list to store the final output rows
    output_rows = []

    # Define custom order for primers (F3, F2, F1, B1, B2, B3)
    primer_order = {'F3': 0, 'F2': 1, 'F1': 2, 'B1': 3, 'B2': 4, 'B3': 5}

    # Process each lamp_id
    for lamp_id, primers in lamp_data.items():
        # Check if there are exactly 6 primers for this lamp_id
        if len(primers) != 6:
            #print(f"Skipping lamp_id {lamp_id} because it does not have exactly 6 primers.")
            continue  # Skip this lamp_id if it does not have exactly 6 primers

        # Map primers to their type (F3, F2, F1, B1, B2, B3) and sort by custom order
        primers_with_type = []
        for primer, sequence, start_pos in primers:
            if 'F3' in primer:
                primer_type = 'F3'
            elif 'F2' in primer:
                primer_type = 'F2'
            elif 'F1' in primer:
                primer_type = 'F1'
            elif 'B1' in primer:
                primer_type = 'B1'
            elif 'B2' in primer:
                primer_type = 'B2'
            elif 'B3' in primer:
                primer_type = 'B3'
            else:
                continue  # Skip primers that do not match F1/F2/F3/B1/B2/B3 pattern

            primers_with_type.append((primer_type, primer, sequence, start_pos))

        # Sort primers by the custom order (F3 > F2 > F1 > B1 > B2 > B3) and start position
        primers_sorted = sorted(primers_with_type, key=lambda x: (primer_order[x[0]], x[3]))

        # Separate primers into different lists based on the order (F3-F1, B1-B3)
        primer_list = [primer[2] for primer in primers_sorted]  # Get sequences
        start_pos_list = [primer[3] for primer in primers_sorted]  # Get start positions

        # Call the score_primers function
        scores = score_primers(primer_list, start_pos_list)

        # Create the row for this lamp_id, with primer-specific scores
        row = {
            'lamp_id': lamp_id,
            'F3_gc_score': scores['single_primer_scores'][0]['gc_score'],
            'F3_tm_score': scores['single_primer_scores'][0]['tm_score'],
            'F3_delta_g_score': scores['single_primer_scores'][0]['delta_g_score'],
            'F3_hairpin_score': scores['single_primer_scores'][0]['hairpin_score'],
            'F2_gc_score': scores['single_primer_scores'][1]['gc_score'],
            'F2_tm_score': scores['single_primer_scores'][1]['tm_score'],
            'F2_delta_g_score': scores['single_primer_scores'][1]['delta_g_score'],
            'F2_hairpin_score': scores['single_primer_scores'][1]['hairpin_score'],
            'F1_gc_score': scores['single_primer_scores'][2]['gc_score'],
            'F1_tm_score': scores['single_primer_scores'][2]['tm_score'],
            'F1_delta_g_score': scores['single_primer_scores'][2]['delta_g_score'],
            'F1_hairpin_score': scores['single_primer_scores'][2]['hairpin_score'],
            'B1_gc_score': scores['single_primer_scores'][3]['gc_score'],
            'B1_tm_score': scores['single_primer_scores'][3]['tm_score'],
            'B1_delta_g_score': scores['single_primer_scores'][3]['delta_g_score'],
            'B1_hairpin_score': scores['single_primer_scores'][3]['hairpin_score'],
            'B2_gc_score': scores['single_primer_scores'][4]['gc_score'],
            'B2_tm_score': scores['single_primer_scores'][4]['tm_score'],
            'B2_delta_g_score': scores['single_primer_scores'][4]['delta_g_score'],
            'B2_hairpin_score': scores['single_primer_scores'][4]['hairpin_score'],
            'B3_gc_score': scores['single_primer_scores'][5]['gc_score'],
            'B3_tm_score': scores['single_primer_scores'][5]['tm_score'],
            'B3_delta_g_score': scores['single_primer_scores'][5]['delta_g_score'],
            'B3_hairpin_score': scores['single_primer_scores'][5]['hairpin_score'],
            'distance_score': scores['distance_score'],
            'dimer_score': scores['dimer_score']
        }

        # Append this row to the output list
        output_rows.append(row)

    # Write the output CSV
    with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
        fieldnames = ['lamp_id', 'F3_gc_score', 'F3_tm_score', 'F3_delta_g_score', 'F3_hairpin_score',
                      'F2_gc_score', 'F2_tm_score', 'F2_delta_g_score', 'F2_hairpin_score',
                      'F1_gc_score', 'F1_tm_score', 'F1_delta_g_score', 'F1_hairpin_score',
                      'B1_gc_score', 'B1_tm_score', 'B1_delta_g_score', 'B1_hairpin_score',
                      'B2_gc_score', 'B2_tm_score', 'B2_delta_g_score', 'B2_hairpin_score',
                      'B3_gc_score', 'B3_tm_score', 'B3_delta_g_score', 'B3_hairpin_score',
                      'distance_score', 'dimer_score']

        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Processed data has been saved to {output_file}")


# Helper function to extract start position
def extract_start_position(position_str):
    # Replace en dash with a regular hyphen
    position_str = position_str.replace('–', '-')
    
    # Now split the position string and convert the first part to an integer
    return int(position_str.split('-')[0])


if __name__ == "__main__":
    input_file = 'data/restructured_output.csv'  # Input CSV path
    output_file = 'data/primers_scoring_output.csv'  # Output CSV path
    process_csv(input_file, output_file)
#
# if __name__ == "__main__":
#     input_file = 'data/nonsense_primers.csv'  # Input CSV path
#     output_file = 'data/nonsense_primers_scoring_output.csv'  # Output CSV path
#     process_csv(input_file, output_file)

if __name__ == "__main__":
    input_file = 'data/nonsense_primers_10w.csv'  # Input CSV path
    output_file = 'data/nonsense_primers_10w_scoring_output.csv'  # Output CSV path
    process_csv(input_file, output_file)