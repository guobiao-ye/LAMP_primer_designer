import math
from typing import Tuple, Any

import primer3

# Dictionary for delta H values  
nn_h = {
    'AA': 91, 'AC': 65, 'AG': 78, 'AT': 86, 'AN': 80,
    'CA': 58, 'CC': 110, 'CG': 119, 'CT': 78, 'CN': 91,
    'GA': 56, 'GC': 111, 'GG': 110, 'GT': 65, 'GN': 85,
    'TA': 60, 'TC': 56, 'TG': 58, 'TT': 91, 'TN': 66,
    'NA': 66, 'NC': 85, 'NG': 91, 'NT': 80, 'NN': 80,
    'aa': 91, 'ac': 65, 'ag': 78, 'at': 86, 'an': 80,
    'ca': 58, 'cc': 110, 'cg': 119, 'ct': 78, 'cn': 91,
    'ga': 56, 'gc': 111, 'gg': 110, 'gt': 65, 'gn': 85,
    'ta': 60, 'tc': 56, 'tg': 58, 'tt': 91, 'tn': 66,
    'na': 66, 'nc': 85, 'ng': 91, 'nt': 80, 'nn': 80,
}

# Dictionary for delta S values  
nn_s = {
    'AA': 240, 'AC': 173, 'AG': 208, 'AT': 239, 'AN': 215,
    'CA': 129, 'CC': 266, 'CG': 278, 'CT': 208, 'CN': 220,
    'GA': 135, 'GC': 267, 'GG': 266, 'GT': 173, 'GN': 210,
    'TA': 169, 'TC': 135, 'TG': 129, 'TT': 240, 'TN': 168,
    'NA': 168, 'NC': 210, 'NG': 220, 'NT': 215, 'NN': 203,
    'aa': 240, 'ac': 173, 'ag': 208, 'at': 239, 'an': 215,
    'ca': 129, 'cc': 266, 'cg': 278, 'ct': 208, 'cn': 220,
    'ga': 135, 'gc': 267, 'gg': 266, 'gt': 173, 'gn': 210,
    'ta': 169, 'tc': 135, 'tg': 129, 'tt': 240, 'tn': 168,
    'na': 168, 'nc': 210, 'ng': 220, 'nt': 215, 'nn': 203,
}


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))


def gc_percent(primer):
    # Ensure primer is not empty before attempting to calculate GC content
    if not primer:
        raise ValueError("Primer sequence is empty.")

    g_count = primer.lower().count('g')
    c_count = primer.lower().count('c')
    total_count = len(primer)

    if total_count == 0:
        raise ValueError("Primer sequence has no valid nucleotides.")

    return round(((g_count + c_count) / total_count) * 100, 2)


# def calc_tm_long(dna):
#     """Calculate Tm for long primers (length > 36)."""
#     gc_pct = gc_percent(dna)
#     return round(81.5 + (16.6 * (math.log(50 / 1000.0) / math.log(10))) + (41.0 * (gc_pct / 100)) - (600.0 / len(dna)),2)
#
# def calc_tm_short(primer):
#     """Calculate Tm for short primers (length <= 36)."""
#     dH = 0
#     dS = 108
#     for i in range(len(primer) - 1):
#         pair = primer[i:i+2]
#         dH += nn_h.get(pair, 0)
#         dS += nn_s.get(pair, 0)
#     dH *= -100.0
#     dS *= -0.1
#     tm = (dH / (dS + 1.987 * math.log(100 / 4000000000.0))) - 273.15 + 16.6 * (math.log(50 / 1000.0) / math.log(10))
#     return round(tm,2)


def calc_tm(primer):
    # """Calculate Tm based on primer length."""
    # if len(primer) > 36:
    #     return calc_tm_long(primer)
    # elif 10 <= len(primer) <= 36:
    #     return calc_tm_short(primer)
    # else:
    #     return None
    tm = primer3.calc_tm(
        seq=primer,
        mv_conc=50.0,  # Sodium ion concentration(mM)
        dv_conc=4.0,  # Magnesium ion concentration(mM)
        dntp_conc=0.2,  # dNTP concentration(mM)
        dna_conc=100.0,  # DNA concentration(nM)
    )
    return round(tm, 2)


def calc_terminal_delta_g(primer, end_length=6):
    """
    Calculate the free energy (ΔG) of the primer terminal region, evaluating both the 5' and 3' ends' stability.
    :param primer: Primer sequence
    :param end_length: Length of the terminal region (default is 6 bases)
    :return: (5' end free energy, 3' end free energy)
    """
    # Ensure primer length is greater than or equal to end_length
    if len(primer) < end_length:
        raise ValueError("Primer length must be greater than or equal to the end_length.")

    # Calculate terminal region sequences
    five_prime_seq = primer[:end_length]  # Get the first end_length bases for 5' end
    three_prime_seq = primer[-end_length:]  # Get the last end_length bases for 3' end

    def calculate_delta_g(sequence):
        """
        Calculate the free energy (ΔG) for a given sequence
        :param sequence: Input sequence (5' or 3' end)
        :return: Free energy value
        """
        dH = 0  # Enthalpy change
        dS = 0  # Entropy change

        # Calculate dH and dS for the sequence
        for i in range(len(sequence) - 1):
            pair = sequence[i:i + 2]  # Get adjacent base pairs
            dH += nn_h.get(pair, 0)  # Lookup enthalpy value for the pair
            dS += nn_s.get(pair, 0)  # Lookup entropy value for the pair

        # Convert units: dH from cal/mol to kcal/mol, dS from cal/mol·K to kcal/mol·K
        dH *= -0.1
        dS *= -0.1

        # Calculate free energy (ΔG) in kcal/mol
        T = 333.15  # 60°C
        delta_g = dH - T * dS / 1000.0  # ΔG formula

        return round(delta_g, 2)

    # Calculate ΔG for 5' and 3' ends separately
    five_prime_delta_g = round(calculate_delta_g(five_prime_seq), 2)
    three_prime_delta_g = round(calculate_delta_g(three_prime_seq), 2)

    return five_prime_delta_g, three_prime_delta_g


def if_secondary_structure(primer: str) -> Tuple[bool, Any]:
    # Define minimum and maximum complementary region lengths for hairpin structure
    min_stem_length = 4
    max_stem_length = 12

    # Define the length range for loop regions
    min_loop_length = 4
    max_loop_length = 8

    if_hairpin = False

    # Scan through the primer sequence for self-complementary regions
    for length in range(min_stem_length, max_stem_length + 1):
        for i in range(len(primer) - length):
            subseq = primer[i:i + length]
            rev_complement_subseq = reverse_complement(subseq)

            # Check if the remaining sequence contains the reverse complementary sequence
            remaining_seq = primer[i + length:]
            index = remaining_seq.find(rev_complement_subseq)

            if index != -1:
                # If complementary sequence found, calculate the distance between the start of the second sequence
                # and the end of the first sequence
                distance = (i + length + index) - (i + length)

                # Check if the loop region length is within the specified range
                if min_loop_length <= distance <= max_loop_length:
                    if_hairpin = True

    hairpin = primer3.calc_hairpin(seq=primer)
    hairpin_dG = round(hairpin.dg / 1000, 2)

    return if_hairpin, hairpin_dG
