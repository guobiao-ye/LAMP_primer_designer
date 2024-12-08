import primer3

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