def generate_mutation_list(sequence, quasispecies_sequences):
    """
    Generates a mutation count list with case-insensitive comparison.

    Parameters:
    - sequence (str): The main sequence.
    - quasispecies_sequences (list of str): List of quasispecies sequences.

    Returns:
    - list of int: Mutation count list.
    """
    sequence_length = len(sequence)
    mutation_list = [0] * sequence_length  # Initialize list with the same length as the main sequence
    
    count = 1
    for qs_seq in quasispecies_sequences:
        qs_seq = qs_seq.upper()
        if len(qs_seq) != sequence_length:
            # Skip sequences with different lengths
            continue
        count += 1
        for i in range(sequence_length):
            if qs_seq[i] != sequence[i]:
                mutation_list[i] += 1

    return [x / count for x in mutation_list]

def filter_forward_oligos(oligos, sequence, mutation_list, weight = 0.1):
    """
    Filters oligos based on sequence matching and adjusts their penalties based on mutation frequencies.

    This function performs two main tasks:
    1. Filters out oligos whose sequences do not match the target sequence at their specified positions.
    2. Adjusts the 'penalty' of each surviving oligo by subtracting the sum of mutation frequencies
       at the positions covered by the oligo, multiplied by the given weight.

    Parameters:
    - oligos (list of dict): A list of oligo dictionaries.
    - sequence (str): 
        The target sequence (e.g., coding strand) to which oligos are compared.
    - mutation_list (list of int or float): 
        A list where each index corresponds to a position in the target sequence and 
        the value represents the mutation frequency at that position.
    - weight (float): 
        The weight factor used to adjust the penalty based on mutation frequencies.

    Returns:
    - list of dict: 
        A list of oligo dictionaries that match the target sequence and have adjusted penalties.
    """
    filtered_oligos = []

    position = oligos[0]['position']
    length = oligos[0]['length']
    oligo_seq = oligos[0]['sequence']

    # Extract the expected sequence from the target sequence
    expected_sequence = sequence[position:position + length]

    # Check if the oligo's sequence matches the expected sequence
    assert oligo_seq == expected_sequence

    for oligo in oligos:
        position = oligo['position']
        length = oligo['length']

        
        # Calculate the sum of mutation frequencies within the oligo's range
        mutation_sum = sum(mutation_list[position:position + length])

        # Adjust the oligo's penalty
        adjusted_penalty = oligo.get('penalty') + (mutation_sum * weight)

        # Update the oligo's penalty in the dictionary
        oligo['penalty'] = adjusted_penalty

        # Add the adjusted oligo to the filtered list
        filtered_oligos.append(oligo)

    return filtered_oligos


def filter_reverse_oligos(oligos, sequence, mutation_list, weight = 0.1):
    """
    Filters oligos based on sequence matching and adjusts their penalties based on mutation frequencies.

    This function performs two main tasks:
    1. Filters out oligos whose sequences do not match the target sequence at their specified positions.
    2. Adjusts the 'penalty' of each surviving oligo by subtracting the sum of mutation frequencies
       at the positions covered by the oligo, multiplied by the given weight.

    Parameters:
    - oligos (list of dict): A list of oligo dictionaries.
    - sequence (str): 
        The target sequence (e.g., coding strand) to which oligos are compared.
    - mutation_list (list of int or float): 
        A list where each index corresponds to a position in the target sequence and 
        the value represents the mutation frequency at that position.
    - weight (float): 
        The weight factor used to adjust the penalty based on mutation frequencies.

    Returns:
    - list of dict: 
        A list of oligo dictionaries that match the target sequence and have adjusted penalties.
    """
    filtered_oligos = []

    position = oligos[0]['position']
    length = oligos[0]['length']
    oligo_seq = oligos[0]['sequence']
    
    complement_map = str.maketrans("ACGTacgt", "TGCAtgca")

    # Extract the expected sequence from the target sequence
    expected_sequence = sequence[position - length + 1:position + 1]
    expected_sequence = expected_sequence.translate(complement_map)[::-1]
    
    # Check if the oligo's sequence matches the expected sequence
    assert oligo_seq == expected_sequence

    for oligo in oligos:
        position = oligo['position']
        length = oligo['length']

        # Calculate the sum of mutation frequencies within the oligo's range
        mutation_sum = sum(mutation_list[position - length + 1:position + 1])

        # Adjust the oligo's penalty
        adjusted_penalty = oligo.get('penalty') + (mutation_sum * weight)

        # Update the oligo's penalty in the dictionary
        oligo['penalty'] = adjusted_penalty

        # Add the adjusted oligo to the filtered list
        filtered_oligos.append(oligo)

    return filtered_oligos
