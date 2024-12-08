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
