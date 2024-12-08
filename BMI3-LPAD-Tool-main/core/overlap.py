from utils.oligos import sort_oligos

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