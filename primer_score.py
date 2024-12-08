import math
from calc_features_single import calc_tm, gc_percent, calc_terminal_delta_g, if_secondary_structure
from calc_features_multi import primers_distance, if_dimer

def score_primers(primer_list, start_pos_list):
    """
    Calculate the weighted score for a combination of 6 primers and return data for visualization.
    :param primer_list: List of 6 primer sequences
    :param start_pos_list: List of starting positions for the 6 primers
    :return: A dictionary containing various primer features for subsequent chart plotting
    """
    primer_set = ['F3', 'F2', 'F1c', 'B1c', 'B2', 'B3']

    # Initialize the return data dictionary
    primer_scores = []

    # Define scoring constants
    GC_CONTENT_MIN = 40  # Minimum GC content
    GC_CONTENT_MAX = 65  # Maximum GC content
    GC_CONTENT_OPT = (50, 60)  # Optimal GC content range
    Tm_F1c_B1c_MIN = 64  # Minimum Tm for F1c and B1c
    Tm_F1c_B1c_MAX = 66  # Maximum Tm for F1c and B1c
    Tm_F2_B2_F3_B3_LOOP_MIN = 59  # Minimum Tm for F2, B2, F3, B3, and Loop primers
    Tm_F2_B2_F3_B3_LOOP_MAX = 61  # Maximum Tm for F2, B2, F3, B3, and Loop primers
    DELTA_G_MAX = -4  # Free energy requirement for the 3' and 5' ends (in kcal/mol)
    DISTANCE_F2_B2_MIN = 120  # Minimum distance between F2 and B2
    DISTANCE_F2_B2_MAX = 160  # Maximum distance between F2 and B2
    DISTANCE_F2_F1_MIN = 40  # Minimum distance between F2 and F1
    DISTANCE_F2_F1_MAX = 60  # Maximum distance between F2 and F1
    DISTANCE_F2_F3_MAX = 60  # Maximum distance between F2 and F3
    DISTANCE_B2_B1_MIN = 40  # Minimum distance between B2 and B1
    DISTANCE_B2_B1_MAX = 60  # Maximum distance between B2 and B1
    DISTANCE_B2_B3_MAX = 60  # Maximum distance between B2 and B3

    # Calculate the features for each primer
    for i, primer in enumerate(primer_list):
        tm = calc_tm(primer)
        gc = gc_percent(primer)
        terminal_delta_g_5p, terminal_delta_g_3p = calc_terminal_delta_g(primer)
        secondary_structure = if_secondary_structure(primer)

        # Calculate the GC score for the primer
        gc_score = 0.7 if GC_CONTENT_MIN <= gc <= GC_CONTENT_MAX else 0
        if GC_CONTENT_OPT[0] <= gc <= GC_CONTENT_OPT[1]:
            gc_score += 0.3

        def calculate_tm_score(tm, tm_min, tm_max):
            if tm_min <= tm <= tm_max:
                return 1  # Full score within the range
            elif tm < tm_min - 5:
                return 0  # 5 degrees below the minimum, score 0
            elif tm > tm_max + 5:
                return 0  # 5 degrees above the maximum, score 0
            else:
                # Calculate the degree of deviation
                if tm < tm_min:
                    # Deviation below the minimum, construct a non-linear function
                    score = 1 - ((tm_min - tm) ** 2 / (5 ** 2))
                else:
                    # Deviation above the maximum
                    score = 1 - ((tm - tm_max) ** 2 / (5 ** 2))

                return round(max(0, min(1, score)), 2)  # Ensure the score is between 0 and 1

        # Apply Tm score calculation
        if i == 2 or i == 3:  # F1c, B1c
            tm_score = calculate_tm_score(tm, Tm_F1c_B1c_MIN, Tm_F1c_B1c_MAX)
        else:  # F2, B2, F3, B3, Loop
            tm_score = calculate_tm_score(tm, Tm_F2_B2_F3_B3_LOOP_MIN, Tm_F2_B2_F3_B3_LOOP_MAX)

        # Calculate the free energy score for the primer
        delta_g_score = 0
        if i == 2 or i == 3:
            if terminal_delta_g_5p <= DELTA_G_MAX:
                delta_g_score += 1
        else:
            if terminal_delta_g_3p <= DELTA_G_MAX:
                delta_g_score += 1

        # Hairpin structure score
        hairpin_score = 0 if secondary_structure[0] else 0.5
        if secondary_structure[1] < 4.5:
            hairpin_score += 0.5

        # Save the score for each primer
        primer_scores.append({
            'primer': primer_set[i],
            'sequence': primer,
            'gc': gc,
            'tm': tm,
            'delta_g_5p': terminal_delta_g_5p,
            'delta_g_3p': terminal_delta_g_3p,
            'hairpin': secondary_structure,
            'gc_score': gc_score,
            'tm_score': tm_score,
            'delta_g_score': delta_g_score,
            'hairpin_score': hairpin_score,
        })

    # Calculate the distances between primers (multi-primer features)
    distances = primers_distance(primer_list, start_pos_list)
    distance_F2_B2_score = round(max(0, min(1, (DISTANCE_F2_B2_MAX - distances['distance_F2_to_B2_end']) / (
            DISTANCE_F2_B2_MAX - DISTANCE_F2_B2_MIN))), 2)
    distance_F2_F1_score = round(max(0, min(1, (DISTANCE_F2_F1_MAX - distances['distance_F2_to_F1']) / (
            DISTANCE_F2_F1_MAX - DISTANCE_F2_F1_MIN))), 2)
    distance_F2_F3_score = round(
        max(0, min(1, (DISTANCE_F2_F3_MAX - distances['distance_F2_to_F3']) / DISTANCE_F2_F3_MAX)), 2)
    distance_B2_B1_score = round(max(0, min(1, (DISTANCE_B2_B1_MAX - distances['distance_B2_to_B1']) / (
            DISTANCE_B2_B1_MAX - DISTANCE_B2_B1_MIN))), 2)
    distance_B2_B3_score = round(
        max(0, min(1, (DISTANCE_B2_B3_MAX - distances['distance_B2_to_B3']) / DISTANCE_B2_B3_MAX)), 2)

    distance_score = round(
        (distance_F2_B2_score + distance_F2_F1_score + distance_F2_F3_score + distance_B2_B1_score + distance_B2_B3_score)/5,
        2)

    # Calculate the dimer score
    dimer_count = if_dimer(primer_list)
    dimer_score = round(max(0, 1 - dimer_count / 3), 2)  # Full score when no dimers, subtract score if dimers exist

    # Calculate dimer status
    dimer_count = if_dimer(primer_list)

    # Return all data for visualization
    return {
        'single_primer_scores': primer_scores,
        'distances': distances,
        'distance_score': distance_score,
        'dimer_count': dimer_count,
        'dimer_score': dimer_score,
    }


if __name__ == "__main__":
    # Sample test: 6 primers and starting positions
    primer_list = ["CTGGTTGTCAAACAACTGG",
                   "TAATAATCTTGGTGGCGTTG",
                   "TACCATAACCAGCTGCTGCG",
                   "TCAAGTGCAAACTCTAGCCGT",
                   "CAGCAGCACCAAGAACTG",
                   "TTCTCTTTCTGGTCAAGGTA"]
    start_pos_list = [50, 69, 142, 171, 200, 226]
    score = score_primers(primer_list, start_pos_list)
    print(f"Total Primer Score: {score}")
