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
    #print(f"Minimum position before: {min_before}, Minimum position after: {min_after}")

    return updated_oligos