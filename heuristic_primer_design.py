import primer3
from Bio import SeqIO
import csv


"""
------------------------------------------------------------------------------------------------
Prepare parameters for LAMP primers
------------------------------------------------------------------------------------------------
"""
outer_primer_params = {
    "PRIMER_INTERNAL_OPT_SIZE": 20,  # outer_primer_target_length
    "PRIMER_INTERNAL_MIN_SIZE": 18,  # outer_primer_min_length
    "PRIMER_INTERNAL_MAX_SIZE": 23,  # outer_primer_max_length
    "PRIMER_INTERNAL_OPT_TM": 60.0,  # outer_primer_target_tm
    "PRIMER_INTERNAL_MAX_POLY_X": 5,  # max_poly_bases
    "PRIMER_INTERNAL_DNA_CONC": 400,  # dna_conc
    "PRIMER_INTERNAL_SALT_MONOVALENT": 50,  # salt_monovalent
    "PRIMER_INTERNAL_SALT_DIVALENT": 8,  # salt_divalent
    "PRIMER_INTERNAL_DNTP_CONC": 1.4,  # dntp_conc
    "PRIMER_NUM_RETURN": 8888,  # max_primer_gen
}


loop_primer_params = {
    "PRIMER_INTERNAL_OPT_SIZE": 20,  # loop_primer_target_length
    "PRIMER_INTERNAL_MIN_SIZE": 18,  # loop_primer_min_length
    "PRIMER_INTERNAL_MAX_SIZE": 23,  # loop_primer_max_length
    "PRIMER_INTERNAL_OPT_TM": 60.0,  # loop_primer_target_tm
    "PRIMER_INTERNAL_MAX_POLY_X": 5,  # max_poly_bases
    "PRIMER_INTERNAL_DNA_CONC": 400,  # dna_conc
    "PRIMER_INTERNAL_SALT_MONOVALENT": 50,  # salt_monovalent
    "PRIMER_INTERNAL_SALT_DIVALENT": 8,  # salt_divalent
    "PRIMER_INTERNAL_DNTP_CONC": 1.4,  # dntp_conc
    "PRIMER_NUM_RETURN": 8888,  # max_primer_gen
}

middle_primer_params = {
    "PRIMER_INTERNAL_OPT_SIZE": 20,  # middle_primer_target_length
    "PRIMER_INTERNAL_MIN_SIZE": 18,  # middle_primer_min_length
    "PRIMER_INTERNAL_MAX_SIZE": 23,  # middle_primer_max_length
    "PRIMER_INTERNAL_OPT_TM": 60.0,  # middle_primer_target_tm
    "PRIMER_INTERNAL_MAX_POLY_X": 5,  # max_poly_bases
    "PRIMER_INTERNAL_DNA_CONC": 400,  # dna_conc
    "PRIMER_INTERNAL_SALT_MONOVALENT": 50,  # salt_monovalent
    "PRIMER_INTERNAL_SALT_DIVALENT": 8,  # salt_divalent
    "PRIMER_INTERNAL_DNTP_CONC": 1.4,  # dntp_conc
    "PRIMER_NUM_RETURN": 8888,  # max_primer_gen
}

inner_primer_params = {
    "PRIMER_INTERNAL_OPT_SIZE": 23,  # inner_primer_target_length
    "PRIMER_INTERNAL_MIN_SIZE": 20,  # inner_primer_min_length
    "PRIMER_INTERNAL_MAX_SIZE": 26,  # inner_primer_max_length
    "PRIMER_INTERNAL_OPT_TM": 62.0,  # inner_primer_target_tm
    "PRIMER_INTERNAL_MAX_POLY_X": 5,  # max_poly_bases
    "PRIMER_INTERNAL_DNA_CONC": 400,  # dna_conc
    "PRIMER_INTERNAL_SALT_MONOVALENT": 50,  # salt_monovalent
    "PRIMER_INTERNAL_SALT_DIVALENT": 8,  # salt_divalent
    "PRIMER_INTERNAL_DNTP_CONC": 1.4,  # dntp_conc
    "PRIMER_NUM_RETURN": 8888,  # max_primer_gen
}




"""
------------------------------------------------------------------------------------------------
Step 1: Generate the DNA oligonucleotide primer sets independently

1. Consider the target sequence as the sense strand. Use primer3 to generate four sets of DNA 
oligonucleotide primers (F2, B2, F3, B3) for the sense strand and the remaining two sets of DNA
oligonucleotide primers (F1c, B1c) for the antisense strand. Save the information for each 
potential primer within each set of results.

2. Perform the above operations on quasi-species sequences and compare them with the results 
from the target sequence to ensure that the primer sets for the target sequence match regions 
that are evolutionarily conserved.

3. Sort the six primer sets for the F1c, B1c, F2, B2, F3, and B3 regions internally based on 
their position and penalty score, and remove primers with excessive overlap according to the 
preset threshold.
------------------------------------------------------------------------------------------------
"""
def run_primer3(sequence, primer_params, reverse=False):
    """
    Generate oligos for the coding strand or template strand based on the reverse flag.

    :param sequence: The DNA sequence (coding strand) as a string.
    :param primer_params: Dictionary of Primer3 parameters for designing primers.
    :param reverse: Boolean indicating whether to design primers for the reverse strand.
    :return: List of oligos designed for the specified strand.
    """

    def reverse_complement(seq):
        """Generate the reverse complement of a DNA sequence."""
        complement_map = str.maketrans("ACGTacgt", "TGCAtgca")
        return seq.translate(complement_map)[::-1]

    # Default global Primer3 parameters
    default_global_params = {
        'PRIMER_TASK': 'pick_hyb_probe_only',
        'PRIMER_INTERNAL_OPT_SIZE': 20,
        'PRIMER_INTERNAL_MIN_SIZE': 18,
        'PRIMER_INTERNAL_MAX_SIZE': 27,
        'PRIMER_INTERNAL_OPT_TM': 60.0,
        'PRIMER_INTERNAL_MIN_TM': 50.0,
        'PRIMER_INTERNAL_MAX_TM': 69.0,
        'PRIMER_INTERNAL_MAX_POLY_X': 4,
        'PRIMER_INTERNAL_DNA_CONC': 400,
        'PRIMER_INTERNAL_SALT_MONOVALENT': 50,
        'PRIMER_INTERNAL_SALT_DIVALENT': 8,
        'PRIMER_INTERNAL_DNTP_CONC': 1.4,
        'PRIMER_NUM_RETURN': 8888,
    }

    # Merge user-provided primer_params with default_global_params
    for key, value in default_global_params.items():
        primer_params.setdefault(key, value)

    # Reverse the sequence if designing for the template strand
    if reverse:
        sequence = reverse_complement(sequence)

    # Base inputs for Primer3
    seq_args = {
        'SEQUENCE_ID': 'target_sequence',
        'SEQUENCE_TEMPLATE': sequence,
    }

    # Run Primer3 using design_primers
    results = primer3.design_primers(seq_args, global_args=primer_params)

    return results

def extract_oligos(primer3_results, sequence_length=None, is_reverse=False):
    """
    Extract oligos from Primer3 results with sequence information, 
    adjusting position for reverse primers.

    :param primer3_results: Dictionary returned by Primer3 with oligo details.
    :param sequence_length: Length of the sequence, required if is_reverse=True.
    :param is_reverse: Boolean indicating whether the oligos are reverse primers.
    :return: List of dictionaries with oligo fields including adjusted position.
    """
    if is_reverse and sequence_length is None:
        raise ValueError("sequence_length must be provided for reverse primers.")

    oligos = []
    for result in primer3_results['PRIMER_INTERNAL']:
        position = result['COORDS'][0]
        length = result['COORDS'][1]

        # Adjust position if reverse primer
        if is_reverse:
            position = sequence_length - position - 1

        oligo = {
            'position': position,
            'length': length,
            'penalty': result['PENALTY'],
            'TM': result['TM'],
            'sequence': result['SEQUENCE']
        }
        oligos.append(oligo)

    return oligos

def filter_forward_oligos(oligos, sequence, quasispecies_sequences):
    """
    Filter oligos based on their sequence matching across all sequences in the list. Only used when evolutionary conservation is considered

    :param oligos: List of oligos with position, length, and sequence.
    :param sequence: Target sequence (e.g., coding strand).
    :param quasispecies_sequences: List of quasispecies sequences (e.g., multiple sequence alignment).
    :return: List of oligos that match in all sequences.
    """
    filtered_oligos = oligos.copy()

    for quasispecies_sequence in quasispecies_sequences:
        surviving_oligos = []
        for oligo in filtered_oligos:
            position = oligo['position']
            length = oligo['length']
            expected_sequence = sequence[position:position + length]

            # Extract the corresponding region from the current sequence
            actual_sequence = quasispecies_sequence[position:position + length]

            # Check if the sequence matches
            if actual_sequence == expected_sequence:
                surviving_oligos.append(oligo)

        # Update the filtered_oligos to only include those that survived this sequence
        filtered_oligos = surviving_oligos

    return filtered_oligos

def filter_reverse_oligos(oligos, sequence, quasispecies_sequences):
    """
    Filter oligos based on their sequence matching across all sequences in the list. Only used when evolutionary conservation is considered

    :param oligos: List of oligos with position, length, and sequence.
    :param sequence: Target sequence (e.g., coding strand).
    :param quasispecies_sequences: List of quasispecies sequences (e.g., multiple sequence alignment).
    :return: List of oligos that match in all sequences.
    """
    filtered_oligos = oligos.copy()

    for quasispecies_sequence in quasispecies_sequences:
        surviving_oligos = []
        for oligo in filtered_oligos:
            position = oligo['position']
            length = oligo['length']
            expected_sequence = sequence[position:position - length]

            # Extract the corresponding region from the current sequence
            actual_sequence = quasispecies_sequence[position:position - length]

            # Check if the sequence matches
            if actual_sequence == expected_sequence:
                surviving_oligos.append(oligo)

        # Update the filtered_oligos to only include those that survived this sequence
        filtered_oligos = surviving_oligos

    return filtered_oligos

def sort_oligos(oligos):
    """
    Sort oligos by location and penalty.

    :param oligos: List of oligo dictionaries with 'position', 'length', 'penalty', and 'TM'.
    :return: Tuple of sorted lists:
             - sorted by location
             - sorted by penalty
    """

    # Sort by location (ascending)
    sort_by_location = sorted(oligos, key=lambda o: o['position'])

    # Sort by penalty (ascending)
    sort_by_penalty = sorted(oligos, key=lambda o: o['penalty'])

    return (
        sort_by_location,
        sort_by_penalty
    )

def reduce_primers_by_overlap(oligos, max_overlap_percent, is_reverse=False):
    """
    Reduce primers by overlap threshold.

    :param oligos_by_location: List of oligos sorted by position.
    :param oligos_by_penalty: List of oligos sorted by penalty.
    :param max_overlap_percent: Maximum allowable overlap percentage.
    :param is_last: Boolean indicating if position represents the last base of the oligo.
    :return: List of reduced primers.
    """
    oligos_by_location, oligos_by_penalty = sort_oligos(oligos)

    if len(oligos_by_location) != len(oligos_by_penalty):
        raise ValueError("Lists must have the same length.")

    # Allow 100% overlap: return location-sorted oligos directly
    if max_overlap_percent == 100:
        return sorted(oligos_by_location, key=lambda x: x['position'])

    unavailable = set()
    selected_primers = []

    for primer in oligos_by_penalty:
        # Determine start and end based on is_reverse
        if is_reverse:
            start = primer['position'] - primer['length'] + 1
            end = primer['position']
        else:
            start = primer['position']
            end = primer['position'] + primer['length'] - 1

        if (start, end) in unavailable:
            continue

        # Add primer to the final list
        selected_primers.append(primer)

        # Calculate overlap range
        upstream_start = start - (2 * primer['length'])
        downstream_end = end

        for oligo in oligos_by_location:
            # Determine oligo start and end based on is_reverse
            if is_reverse:
                oligo_start = oligo['position'] - oligo['length'] + 1
                oligo_end = oligo['position']
            else:
                oligo_start = oligo['position']
                oligo_end = oligo['position'] + oligo['length'] - 1

            # Skip unavailable oligos
            if (oligo_start, oligo_end) in unavailable:
                continue

            # Check if oligo is within range
            if oligo_end < upstream_start or oligo_start > downstream_end:
                continue

            # Calculate overlap percentage
            overlap = max(0, min(end, oligo_end) - max(start, oligo_start) + 1)
            overlap_percent = (overlap / min(primer['length'], oligo['length'])) * 100

            # Mark oligo as unavailable if overlap exceeds threshold
            if overlap_percent > max_overlap_percent:
                unavailable.add((oligo_start, oligo_end))

    return sorted(selected_primers, key=lambda x: x['position'])

def move_position(oligos, mode='end_to_start'):
    """
    Convert the 'position' of each oligo between end and start based on the mode,
    and update the sequence to its reverse complement.

    :param oligos: List of oligos, where each oligo is a dictionary with 'position', 'length', and 'sequence'.
    :param mode: Conversion mode, either 'end_to_start' or 'start_to_end'.
    :return: A new list of oligos with 'position' and 'sequence' updated.
    """
    if not oligos:
        print("No oligos to process.")
        return []

    if mode not in ['end_to_start', 'start_to_end']:
        raise ValueError("Invalid mode. Use 'end_to_start' or 'start_to_end'.")

    # DNA base complement mapping
    complement_map = str.maketrans("ACGTacgt", "TGCAtgca")

    min_before = float('inf')
    min_after = float('inf')

    # Create a new list to avoid modifying the original input
    updated_oligos = []

    for primer in oligos:
        if 'position' not in primer or 'length' not in primer or 'sequence' not in primer:
            raise ValueError("Each oligo must have 'position', 'length', and 'sequence' keys.")

        # Track minimum position before modification
        if primer['position'] < min_before:
            min_before = primer['position']

        # Calculate the new position based on the mode
        if mode == 'end_to_start':
            new_position = primer['position'] - primer['length'] + 1
        elif mode == 'start_to_end':
            new_position = primer['position'] + primer['length'] - 1

        # Track minimum position after modification
        if new_position < min_after:
            min_after = new_position

        # Convert sequence to reverse complement
        reverse_complement_sequence = primer['sequence'].translate(complement_map)[::-1]

        # Update the primer dictionary
        updated_oligo = primer.copy()
        updated_oligo['position'] = new_position
        updated_oligo['sequence'] = reverse_complement_sequence
        updated_oligos.append(updated_oligo)

    # Debugging output for min_before and min_after
    print(f"Minimum position before: {min_before}, Minimum position after: {min_after}")

    return updated_oligos




"""
------------------------------------------------------------------------------------------------
Step 2: Combine the primers to generate the LAMP primer sets

Given a series of primers and a set of constraints, use a nested loop structure to sequentially 
iterate through inner primers, middle primers, outer primers, and loop primers to find the 
optimal primer combinations.

1. According to the distance restrictions between LAMP primers, separately obtain possible 
primer combination sets for both ends of the DNA (F3, F2, F1c on one end and B3, B2, B1c on the
other end).

2. For each combination, calculate the sum of the penalty scores for the DNA oligonucleotide 
primers and the distance penalties to determine the information of the optimal primer 
combinations along with their corresponding penalty scores. Remove primer combinations with 
excessive overlap according to a predefined threshold.

3. Finally, combine the two primer sets (F3, F2, F1c and B3, B2, B1c) based on using a LAMP
primer combination scoring system to obtain the final set of LAMP primer combinations.
------------------------------------------------------------------------------------------------
"""
def find_best_forward_combinations(
    inner_forward, loop_forward, middle_forward, outer_forward,
    signature_max_length, min_primer_spacing, loop_min_gap,
    include_loop_primers, distance_penalties,
    penalty_weights,  # [inner, loop, middle, outer]
    to_penalty_weights  # [inner_to_loop, loop_to_middle, middle_to_outer]
):
    """
    Find the best forward combinations based on penalties and spacing constraints.

    :param inner_forward: List of inner primers.
    :param loop_forward: List of loop primers.
    :param middle_forward: List of middle primers.
    :param outer_forward: List of outer primers.
    :param signature_max_length: Maximum allowed length of a signature.
    :param min_primer_spacing: Minimum spacing between primers.
    :param loop_min_gap: Minimum gap for loop primers.
    :param include_loop_primers: Whether to include loop primers.
    :param distance_penalties: List of penalties for inter-primer distances.
    :param penalty_weights: List of weights for [inner, loop, middle, outer] primers.
    :param to_penalty_weights: List of weights for [inner_to_loop, loop_to_middle, middle_to_outer] distances.
    :return: Tuple of best forward infos, penalties, and the count of forward sets.
    """
    # Unpack weights for individual penalties
    inner_penalty_weight, loop_penalty_weight, middle_penalty_weight, outer_penalty_weight = penalty_weights

    # Unpack weights for inter-primer penalties
    inner_to_loop_penalty_weight, loop_to_middle_penalty_weight, middle_to_outer_penalty_weight = to_penalty_weights

    best_forward_infos = []
    best_forward_penalties = []
    forward_set_count = 0

    for inner_index, inner_info in enumerate(inner_forward):
        inner_location = inner_info['position']
        inner_length = inner_info['length']
        inner_penalty = inner_info['penalty']

        best_set_penalty = float('inf')  # Start with a large value

        # Calculate search range for loop primers
        search_start_at = max(0, inner_location - signature_max_length + inner_length + 20)
        loop_start_at = search_start_at
        loop_end_at = max(0, inner_location - 1 - min_primer_spacing)

        # Handle loop primer placeholder if loop primers are not included
        if not include_loop_primers:
            placeholder_primer = {'position': loop_end_at + 1, 'length': 1, 'penalty': 0}
            loop_forward = [placeholder_primer]

        for loop_index, loop_info in enumerate(loop_forward):
            loop_location = loop_info['position']
            loop_length = loop_info['length']
            loop_penalty = loop_info['penalty']

            # Check loop primer position range
            if loop_location < loop_start_at and loop_length != 1:
                continue
            if loop_location > loop_end_at and loop_length != 1:
                break

            # Calculate range for middle primers
            middle_start_at = search_start_at
            middle_end_at = min(
                loop_location - loop_length - min_primer_spacing,
                inner_location - loop_min_gap - 1
            )
            middle_end_at = max(0, middle_end_at)

            inner_to_loop_distance = inner_location - (loop_location + 1)

            for middle_index, middle_info in enumerate(middle_forward):
                middle_location = middle_info['position']
                middle_length = middle_info['length']
                middle_penalty = middle_info['penalty']

                # Check middle primer position range
                if middle_location < middle_start_at:
                    continue
                if middle_location > middle_end_at:
                    break

                # Check spacing constraints for middle primers
                if (middle_location + middle_length + min_primer_spacing > loop_location - loop_length + 1) or \
                   (middle_location + middle_length + loop_min_gap > inner_location):
                    continue

                # Calculate range for outer primers
                outer_start_at = search_start_at
                outer_end_at = middle_location - 1 - min_primer_spacing

                loop_to_middle_distance = (loop_location - loop_length + 1) - (middle_location + middle_length)

                for outer_index, outer_info in enumerate(outer_forward):
                    outer_location = outer_info['position']
                    outer_length = outer_info['length']
                    outer_penalty = outer_info['penalty']

                    # Check outer primer position range
                    if outer_location < outer_start_at:
                        continue
                    if outer_location > outer_end_at:
                        break

                    # Check spacing constraints for outer primers
                    if outer_location + outer_length + min_primer_spacing > middle_location:
                        continue

                    middle_to_outer_distance = middle_location - (outer_location + outer_length)
                    inner_to_middle_distance = inner_location - (middle_location + middle_length)

                    # Calculate penalties
                    if include_loop_primers:
                        spacing_penalty = (
                            distance_penalties[inner_to_loop_distance - 15] * inner_to_loop_penalty_weight +
                            distance_penalties[loop_to_middle_distance - 15] * loop_to_middle_penalty_weight +
                            distance_penalties[middle_to_outer_distance - 30] * middle_to_outer_penalty_weight
                        )
                        primer3_penalty = (
                            inner_penalty * inner_penalty_weight +
                            loop_penalty * loop_penalty_weight +
                            middle_penalty * middle_penalty_weight +
                            outer_penalty * outer_penalty_weight
                        )
                    else:
                        spacing_penalty = (
                            distance_penalties[inner_to_middle_distance - 30] * inner_to_loop_penalty_weight +
                            distance_penalties[middle_to_outer_distance - 30] * middle_to_outer_penalty_weight
                        )
                        primer3_penalty = (
                            inner_penalty * inner_penalty_weight +
                            middle_penalty * middle_penalty_weight +
                            outer_penalty * outer_penalty_weight
                        )

                    forward_set_penalty = spacing_penalty + primer3_penalty

                    # Update best combination if penalty is lower
                    if forward_set_penalty < best_set_penalty:
                        best_set_penalty = forward_set_penalty
                        while len(best_forward_infos) <= inner_index:
                            best_forward_infos.append(None)
                        while len(best_forward_penalties) <= inner_index:
                            best_forward_penalties.append(None)
                        best_forward_infos[inner_index] = [loop_info, middle_info, outer_info]
                        best_forward_penalties[inner_index] = [spacing_penalty, primer3_penalty]
    
    forward_set_count = sum(1 for item in best_forward_infos if item is not None)

    return best_forward_infos, best_forward_penalties, forward_set_count

def find_best_reverse_combinations(
    inner_reverse, loop_reverse, middle_reverse, outer_reverse,
    signature_max_length, min_primer_spacing, loop_min_gap,
    include_loop_primers, distance_penalties,
    penalty_weights,  # [inner, loop, middle, outer]
    to_penalty_weights  # [inner_to_loop, loop_to_middle, middle_to_outer]
):
    """
    Find the best reverse primer combinations based on penalties and spacing constraints.

    :param inner_reverse: List of inner reverse primers (list of dicts).
    :param loop_reverse: List of loop reverse primers (list of dicts).
    :param middle_reverse: List of middle reverse primers (list of dicts).
    :param outer_reverse: List of outer reverse primers (list of dicts).
    :param signature_max_length: Maximum allowed length of a signature.
    :param min_primer_spacing: Minimum spacing between primers.
    :param loop_min_gap: Minimum gap for loop primers.
    :param include_loop_primers: Whether to include loop primers.
    :param distance_penalties: List of penalties for inter-primer distances.
    :param penalty_weights: List of weights for [inner, loop, middle, outer] primers.
    :param to_penalty_weights: List of weights for [inner_to_loop, loop_to_middle, middle_to_outer] distances.
    :return: Tuple of best reverse infos and penalties.
    """
    # Unpack weights for individual penalties
    inner_penalty_weight, loop_penalty_weight, middle_penalty_weight, outer_penalty_weight = penalty_weights

    # Unpack weights for inter-primer penalties
    inner_to_loop_penalty_weight, loop_to_middle_penalty_weight, middle_to_outer_penalty_weight = to_penalty_weights

    best_reverse_infos = []
    best_reverse_penalties = []
    reverse_set_count = 0

    for inner_index, inner_info in enumerate(inner_reverse):
        inner_location = inner_info['position']
        inner_length = inner_info['length']
        inner_penalty = inner_info['penalty']

        best_set_penalty = float('inf')  # Start with a large initial penalty

        # Calculate search range for loop primers
        search_end_at = inner_location + signature_max_length - inner_length - 20
        loop_start_at = inner_location + 1 + min_primer_spacing
        loop_end_at = search_end_at

        # Handle loop primer placeholder if loop primers are not included
        if not include_loop_primers:
            placeholder_primer = {'position': loop_start_at - 1, 'length': 1, 'penalty': 0}
            loop_reverse = [placeholder_primer]

        for loop_info in loop_reverse:
            loop_location = loop_info['position']
            loop_length = loop_info['length']
            loop_penalty = loop_info['penalty']

            # Check loop primer position range
            if loop_location < loop_start_at and loop_length != 1:
                continue
            if loop_location > loop_end_at and loop_length != 1:
                break

            # Calculate range for middle primers
            middle_start_at = max(loop_location + loop_length + min_primer_spacing,
                                  inner_location + loop_min_gap + 1)
            middle_end_at = search_end_at

            inner_to_loop_distance = loop_location - (inner_location + 1)

            for middle_info in middle_reverse:
                middle_location = middle_info['position']
                middle_length = middle_info['length']
                middle_penalty = middle_info['penalty']

                # Check middle primer position range
                if middle_location < middle_start_at:
                    continue
                if middle_location > middle_end_at:
                    break

                # Check spacing constraints for middle primers
                if (middle_location - middle_length - min_primer_spacing < loop_location + loop_length - 1) or \
                   (middle_location - middle_length - loop_min_gap < inner_location):
                    continue

                # Calculate range for outer primers
                outer_start_at = middle_location + min_primer_spacing + 1
                outer_end_at = search_end_at

                loop_to_middle_distance = (middle_location - middle_length + 1) - \
                                          (loop_location + loop_length)

                for outer_info in outer_reverse:
                    outer_location = outer_info['position']
                    outer_length = outer_info['length']
                    outer_penalty = outer_info['penalty']

                    # Check outer primer position range
                    if outer_location < outer_start_at:
                        continue
                    if outer_location > outer_end_at:
                        break

                    # Check spacing constraints for outer primers
                    if outer_location - outer_length - min_primer_spacing < middle_location:
                        continue

                    middle_to_outer_distance = (outer_location - outer_length) - middle_location
                    inner_to_middle_distance = (middle_location - middle_length) - inner_location

                    # Calculate penalties
                    if include_loop_primers:
                        spacing_penalty = (
                            distance_penalties[inner_to_loop_distance - 15] * inner_to_loop_penalty_weight +
                            distance_penalties[loop_to_middle_distance - 15] * loop_to_middle_penalty_weight +
                            distance_penalties[middle_to_outer_distance - 30] * middle_to_outer_penalty_weight
                        )
                        primer3_penalty = (
                            inner_penalty * inner_penalty_weight +
                            loop_penalty * loop_penalty_weight +
                            middle_penalty * middle_penalty_weight +
                            outer_penalty * outer_penalty_weight
                        )
                    else:
                        spacing_penalty = (
                            distance_penalties[inner_to_middle_distance-30] * inner_to_loop_penalty_weight +
                            distance_penalties[middle_to_outer_distance-30] * middle_to_outer_penalty_weight
                        )
                        primer3_penalty = (
                            inner_penalty * inner_penalty_weight +
                            middle_penalty * middle_penalty_weight +
                            outer_penalty * outer_penalty_weight
                        )

                    reverse_set_penalty = spacing_penalty + primer3_penalty

                    # Update best combination if penalty is lower
                    if reverse_set_penalty < best_set_penalty:
                        best_set_penalty = reverse_set_penalty
                        while len(best_reverse_infos) <= inner_index:
                            best_reverse_infos.append(None)
                        while len(best_reverse_penalties) <= inner_index:
                            best_reverse_penalties.append(None)
                        best_reverse_infos[inner_index] = [loop_info, middle_info, outer_info]
                        best_reverse_penalties[inner_index] = [spacing_penalty, primer3_penalty]
    
    reverse_set_count = sum(1 for item in best_reverse_infos if item is not None)

    return best_reverse_infos, best_reverse_penalties, reverse_set_count

def reduce_result_by_overlap(possible_result, max_overlap_percent=99, sort_by_score=True, sort_by_location=False):
    """
    Reduce possible results based on overlap threshold.

    :param possible_result: List of result dictionaries with primer information and penalties.
    :param max_overlap_percent: Maximum allowable overlap percentage.
    :param sort_by_score: Whether to sort by total penalty score.
    :param sort_by_location: Whether to sort by starting location.
    :return: Reduced list of results.
    """
    # Default sorting preference
    if not sort_by_score and not sort_by_location:
        sort_by_score = True
    if sort_by_score and sort_by_location:
        sort_by_location = False

    # Allow 100% overlap: return sorted by penalty
    if max_overlap_percent == 100:
        return sorted(possible_result, key=lambda x: x["penalty"])

    # Prepare sorted lists
    sorted_by_penalty = sorted(
        enumerate(possible_result),
        key=lambda x: x[1]["penalty"]
    )
    sorted_by_location = sorted(
        enumerate(possible_result),
        key=lambda x: x[1]["forward_outer_info"]["position"]
    )

    unavailable_indices = set()
    selected_results = []

    for index, result in sorted_by_penalty:
        if index in unavailable_indices:
            continue

        # Add current result to selected list
        selected_results.append(result)
        unavailable_indices.add(index)

        # Calculate overlap bounds
        pos, end = result["forward_outer_info"]["position"], result["reverse_outer_info"]["position"]
        length = end - pos + 1
        upstream_start = max(0, pos - 5 * length)
        downstream_end = pos + length - 1

        # Mark overlapping results as unavailable
        for other_index, other_result in sorted_by_location:
            if other_index in unavailable_indices:
                continue

            other_pos, other_end = other_result["forward_outer_info"]["position"], other_result["reverse_outer_info"]["position"]
            other_length = other_pos - other_end + 1

            # No overlap, skip
            if other_pos + other_length - 1 < upstream_start or other_pos > downstream_end:
                continue

            # Calculate overlap percentage
            overlap_start = max(pos, other_pos)
            overlap_end = min(pos + length - 1, other_pos + other_length - 1)
            overlap_length = max(0, overlap_end - overlap_start + 1)
            overlap_percent = (overlap_length / min(length, other_length)) * 100

            if overlap_percent > max_overlap_percent:
                unavailable_indices.add(other_index)

    # Sort final results
    if sort_by_location:
        return sorted(selected_results, key=lambda x: x["forward_outer_info"]["position"])
    return sorted(selected_results, key=lambda x: x["penalty"])

def find_possible_result(
    forward_inner_candidates,
    reverse_inner_candidates,
    best_forward_infos,
    best_reverse_infos,
    best_forward_penalties,
    best_reverse_penalties,
    signature_max_length,
    min_inner_pair_spacing,
    opt_inner_pair_spacing,
    Tm_penalties,
    distance_penalties,
    include_loop_primers,
    output_file_path
):
    # Initialize variables
    possible_result = []
    previous_first_compatible_index = 0  # Bound the lower end of the inner loop

    # Iterate over the forward inner candidates
    for i in range(len(forward_inner_candidates)):
        if i >= len(best_forward_infos) or not best_forward_infos[i]:
            continue

        finner_info = forward_inner_candidates[i]
        floop_info, fmiddle_info, fouter_info = best_forward_infos[i]
        forward_spacing_penalty, forward_primer3_penalty = best_forward_penalties[i]

        forward_start = fouter_info["position"]
        forward_end = finner_info["position"] + finner_info["length"] - 1

        # Bound the upper end of the inner loop search
        max_reverse_location = forward_start + signature_max_length - 1
        previous_compatible_index_found = False

        # Iterate over the reverse inner candidates
        for j in range(previous_first_compatible_index, len(reverse_inner_candidates)):
            if j >= len(best_reverse_infos) or not best_reverse_infos[j]:
                continue

            binner_info = reverse_inner_candidates[j]
            bloop_info, bmiddle_info, bouter_info = best_reverse_infos[j]
            reverse_spacing_penalty, reverse_primer3_penalty = best_reverse_penalties[j]

            reverse_end = bouter_info["position"]
            reverse_start = binner_info["position"] - binner_info["length"] + 1

            # Skip primers located too far 5' with respect to the forward primer
            if not previous_compatible_index_found:
                if reverse_start <= forward_end:
                    continue
                else:
                    previous_first_compatible_index = j
                    previous_compatible_index_found = True

            # Stop searching if the inner loop bounds are exceeded
            if reverse_start > max_reverse_location:
                break

            # Enforce minimum inner spacing distance
            inner_spacing = reverse_start - (forward_end + 1)
            if inner_spacing < min_inner_pair_spacing:
                continue

            # Enforce max signature length
            if reverse_end - (forward_start + 1) > signature_max_length:
                continue

            # Calculate penalties
            inner_spacing_penalty = distance_penalties[inner_spacing - opt_inner_pair_spacing] * 1  # Example weight: 1

            if include_loop_primers:
                Total_Tm_diff = (abs(finner_info['TM'] - binner_info['TM']) + abs(fmiddle_info['TM'] - bmiddle_info['TM']) + abs(fouter_info['TM'] - bouter_info['TM']) + abs(floop_info['TM'] - bloop_info['TM'])) * Tm_penalties
            else:
                Total_Tm_diff = (abs(finner_info['TM'] - binner_info['TM']) + abs(fmiddle_info['TM'] - bmiddle_info['TM']) + abs(fouter_info['TM'] - bouter_info['TM'])) * Tm_penalties

            total_penalty = forward_spacing_penalty + inner_spacing_penalty + reverse_spacing_penalty + forward_primer3_penalty + reverse_primer3_penalty + Total_Tm_diff

            # Construct result object (dictionary-based structure for simplicity)
            result = {
                "forward_inner_info": finner_info,
                "reverse_inner_info": binner_info,
                "forward_middle_info": fmiddle_info,
                "reverse_middle_info": bmiddle_info,
                "forward_outer_info": fouter_info,
                "reverse_outer_info": bouter_info,
                "penalty": total_penalty,
            }

            if include_loop_primers:
                result["floop_info"] = floop_info
                result["bloop_info"] = bloop_info
                result["has_loop_primers"] = True

            possible_result.append(result)

    # Filter result by overlap (mock function, you can replace it with your actual logic)
    possible_result = reduce_result_by_overlap(possible_result)

    if len(possible_result) == 0:
        print("Failed to find result.")
        return

    print(f"Found {len(possible_result)} possible LAMP primer combinations")

    # Sort result by penalty
    possible_result.sort(key=lambda x: x["penalty"])

    # Write result to CSV
    fieldnames = [
        "penalty",
        "forward_inner_info", "reverse_inner_info",
        "forward_middle_info", "reverse_middle_info",
        "forward_outer_info", "reverse_outer_info",
        "floop_info", "bloop_info", "has_loop_primers"
    ]
    with open(output_file_path, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for result in possible_result:
            writer.writerow({
                "penalty": result["penalty"],
                "forward_inner_info": result.get("forward_inner_info"),
                "reverse_inner_info": result.get("reverse_inner_info"),
                "forward_middle_info": result.get("forward_middle_info"),
                "reverse_middle_info": result.get("reverse_middle_info"),
                "forward_outer_info": result.get("forward_outer_info"),
                "reverse_outer_info": result.get("reverse_outer_info"),
                "floop_info": result.get("floop_info", None),
                "bloop_info": result.get("bloop_info", None),
                "has_loop_primers": result.get("has_loop_primers", False),
            })

    print(f"Output written to {output_file_path}")

def generate_distance_penalties(a=20, b=1.01, max_penalty=100):
    """
    Generate distance penalties based on the maximum distance using the custom function 
    f(x) = a * (b^x - 1), and limit the penalty values to a maximum value.

    :param a: Scaling factor for the penalty values.
    :param b: Base of the exponent, b > 1 controls the growth rate.
    :param max_penalty: Maximum allowable penalty value (default is 25).
    :return: A list of distance penalties.
    """
    # Generate penalties using f(x) = a * (b^x - 1) and limit to max_penalty
    penalties = [min(a * (b**x - 1), max_penalty) for x in range(200)]

    return penalties




"""
------------------------------------------------------------------------------------------------
Now we can design LAMP primers by calling the above functions!!!
------------------------------------------------------------------------------------------------
"""
def find_primers(fasta_file):    
    sequence = None
    quasispecies_sequences = []

    # Open the FASTA file and parse the records
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if sequence is None:
                # Set the first sequence as the main sequence
                sequence = str(record.seq)
            else:
                # Add subsequent sequences to the quasispecies list
                quasispecies_sequences.append(str(record.seq))
                
        sequence_length = len(sequence)
        forward_inner_c_oligos = extract_oligos(run_primer3(sequence, inner_primer_params, reverse=True), sequence_length, is_reverse=True)
        reverse_inner_c_oligos = extract_oligos(run_primer3(sequence, inner_primer_params, reverse=False), sequence_length, is_reverse=False)
        forward_inner_oligos = move_position(forward_inner_c_oligos, 'end_to_start')
        reverse_inner_oligos = move_position(reverse_inner_c_oligos, 'start_to_end')
        forward_loop_oligos = extract_oligos(run_primer3(sequence, loop_primer_params, reverse=True), sequence_length, is_reverse=True)
        reverse_loop_oligos = extract_oligos(run_primer3(sequence, loop_primer_params, reverse=False), sequence_length, is_reverse=False)
        forward_middle_oligos = extract_oligos(run_primer3(sequence, middle_primer_params, reverse=False), sequence_length, is_reverse=False)
        reverse_middle_oligos = extract_oligos(run_primer3(sequence, middle_primer_params, reverse=True), sequence_length, is_reverse=True)
        forward_outer_oligos = extract_oligos(run_primer3(sequence, outer_primer_params, reverse=False), sequence_length, is_reverse=False)
        reverse_outer_oligos = extract_oligos(run_primer3(sequence, outer_primer_params, reverse=True), sequence_length, is_reverse=True)

        if quasispecies_sequences:
            forward_outer_oligos = filter_forward_oligos(forward_outer_oligos, sequence, quasispecies_sequences)
            reverse_outer_oligos = filter_reverse_oligos(reverse_outer_oligos, sequence, quasispecies_sequences)
            forward_outer_oligos = filter_forward_oligos(forward_outer_oligos, sequence, quasispecies_sequences)
            reverse_outer_oligos = filter_reverse_oligos(reverse_outer_oligos, sequence, quasispecies_sequences)
            forward_outer_oligos = filter_forward_oligos(forward_outer_oligos, sequence, quasispecies_sequences)
            reverse_outer_oligos = filter_reverse_oligos(reverse_outer_oligos, sequence, quasispecies_sequences)
            forward_outer_oligos = filter_forward_oligos(forward_outer_oligos, sequence, quasispecies_sequences)
            reverse_outer_oligos = filter_reverse_oligos(reverse_outer_oligos, sequence, quasispecies_sequences)
            forward_outer_oligos = filter_forward_oligos(forward_outer_oligos, sequence, quasispecies_sequences)
            reverse_outer_oligos = filter_reverse_oligos(reverse_outer_oligos, sequence, quasispecies_sequences)
            
        print("All oligos found!")
        print(f"Number of oligos are: forward_inner - {len(forward_inner_oligos)}, reverse_inner - {len(reverse_inner_oligos)}, forward_loop - {len(forward_loop_oligos)}, reverse_loop - {len(reverse_loop_oligos)}, forward_middle - {len(forward_middle_oligos)}, reverse_middle - {len(reverse_middle_oligos)}, forward_outer - {len(forward_outer_oligos)}, reverse_outer - {len(reverse_outer_oligos)}")

        forward_inner_candidates = reduce_primers_by_overlap(forward_inner_oligos, 60, False)
        reverse_inner_candidates = reduce_primers_by_overlap(reverse_inner_oligos, 60, True)
        forward_loop_candidates = reduce_primers_by_overlap(forward_loop_oligos, 60, True)
        reverse_loop_candidates = reduce_primers_by_overlap(reverse_loop_oligos, 60, False)
        forward_middle_candidates = reduce_primers_by_overlap(forward_middle_oligos, 60, False)
        reverse_middle_candidates = reduce_primers_by_overlap(reverse_middle_oligos, 60, True)
        forward_outer_candidates = reduce_primers_by_overlap(forward_outer_oligos, 60, False)
        reverse_outer_candidates = reduce_primers_by_overlap(reverse_outer_oligos, 60, True)

        print("All candidates prepared by reducing overlaps!")
        print(f"Number of candidates are:  forward_inner - {len(forward_inner_candidates)}, reverse_inner - {len(reverse_inner_candidates)}, forward_loop - {len(forward_loop_candidates)}, reverse_loop - {len(reverse_loop_candidates)}, forward_middle - {len(forward_middle_candidates)}, reverse_middle - {len(reverse_middle_candidates)}, forward_outer - {len(forward_outer_candidates)}, reverse_outer - {len(reverse_outer_candidates)}")

        distance_penalties = generate_distance_penalties()
        best_forward_infos, best_forward_penalties, forward_set_count = find_best_forward_combinations(forward_inner_candidates,
                                                                                                    forward_loop_candidates,
                                                                                                    forward_middle_candidates,
                                                                                                    forward_outer_candidates,
                                                                                                    320, 1, 25, True, distance_penalties, 
                                                                                                    [1.2, 0.7, 1.1, 1.0], [1, 1, 1])
        print(f'{forward_set_count} sets of forward combinations are generated.')

        best_reverse_infos, best_reverse_penalties, reverse_set_count = find_best_reverse_combinations(reverse_inner_candidates,
                                                                                                    reverse_loop_candidates,
                                                                                                    reverse_middle_candidates,
                                                                                                    reverse_outer_candidates,
                                                                                                    320, 1, 25, True, distance_penalties, 
                                                                                                    [1.2, 0.7, 1.1, 1.0], [1, 1, 1])
        print(f'{reverse_set_count} sets of reverse combinations are generated.')

        find_possible_result(forward_inner_candidates, reverse_inner_candidates, best_forward_infos, best_reverse_infos,
                            best_forward_penalties, best_reverse_penalties, 320, 1, 4, 0.6, distance_penalties, True, "./data/output/Intermediate_file/candidate_result.csv")
            