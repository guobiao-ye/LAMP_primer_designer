def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))


def primers_distance(primer_list, start_pos_list):
    # Primer lengths
    lengths = [len(primer) for primer in primer_list]

    # Calculate the required distances
    distance_F2_to_B2_end = (start_pos_list[4] + lengths[4]) - (start_pos_list[1])  # From F2 to the end of B2
    distance_F2_to_F1 = start_pos_list[2] - start_pos_list[1]  # From the 5' end of F2 to the 5' end of F1
    distance_F2_to_F3 = start_pos_list[1] - (start_pos_list[0] + lengths[0])  # From the 5' end of F2 to the end of F3
    distance_B2_to_B1 = (start_pos_list[4] + lengths[4]) - (
                start_pos_list[3] + lengths[3])  # From the end of B2 to the end of B1
    distance_B2_to_B3 = start_pos_list[5] - (start_pos_list[4] + lengths[4])  # From the end of B2 to the 5' end of B3

    return {
        "distance_F2_to_B2_end": distance_F2_to_B2_end,
        "distance_F2_to_F1": distance_F2_to_F1,
        "distance_F2_to_F3": distance_F2_to_F3,
        "distance_B2_to_B1": distance_B2_to_B1,
        "distance_B2_to_B3": distance_B2_to_B3
    }


def if_dimer(primer_list):
    # Define the minimum and maximum threshold lengths for complementary segments
    min_dimer_length = 8
    max_dimer_length = 16

    # Count the number of dimers that meet the criteria
    dimer_count = 0

    # Iterate through each primer pair (including a primer with itself)
    for i in range(len(primer_list)):
        for j in range(i, len(primer_list)):
            primer1 = primer_list[i]
            primer2 = primer_list[j]

            # Get the reverse complement of primer2 (used to check complementarity)
            rev_complement_primer2 = reverse_complement(primer2)

            # Check for complementary segments
            for length in range(min_dimer_length, max_dimer_length + 1):
                # Find complementary segments between primer1 and rev_complement_primer2
                for k in range(len(primer1) - length + 1):
                    subseq1 = primer1[k:k + length]
                    if subseq1 in rev_complement_primer2:
                        dimer_count += 1  # If a complementary segment is found, increment the count
                        break  # Stop searching for complementary segments once one is found

    return dimer_count  # Return the count of dimers found
