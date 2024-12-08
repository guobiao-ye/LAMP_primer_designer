import csv
from core.overlap import reduce_result_by_overlap

def find_best_forward_combinations(
        inner_forward, loop_forward, middle_forward, outer_forward,
        signature_max_length, min_primer_spacing, loop_min_gap,
        include_loop_primers, ideal_gap, distance_penalties,
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
                                distance_penalties[abs(inner_to_loop_distance-ideal_gap)] * inner_to_loop_penalty_weight +
                                distance_penalties[abs(loop_to_middle_distance-ideal_gap)] * loop_to_middle_penalty_weight +
                                distance_penalties[abs(middle_to_outer_distance-ideal_gap)] * middle_to_outer_penalty_weight
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
                                distance_penalties[abs(middle_to_outer_distance-ideal_gap)] * middle_to_outer_penalty_weight
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
        include_loop_primers, ideal_gap, distance_penalties,
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
                                distance_penalties[abs(inner_to_loop_distance-ideal_gap)] * inner_to_loop_penalty_weight +
                                distance_penalties[abs(loop_to_middle_distance-ideal_gap)] * loop_to_middle_penalty_weight +
                                distance_penalties[abs(middle_to_outer_distance-ideal_gap)] * middle_to_outer_penalty_weight
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
                                distance_penalties[abs(middle_to_outer_distance-ideal_gap)] * middle_to_outer_penalty_weight
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
                Total_Tm_diff = (abs(finner_info['TM'] - binner_info['TM']) + abs(
                    fmiddle_info['TM'] - bmiddle_info['TM']) + abs(fouter_info['TM'] - bouter_info['TM']) + abs(
                    floop_info['TM'] - bloop_info['TM'])) * Tm_penalties
            else:
                Total_Tm_diff = (abs(finner_info['TM'] - binner_info['TM']) + abs(
                    fmiddle_info['TM'] - bmiddle_info['TM']) + abs(
                    fouter_info['TM'] - bouter_info['TM'])) * Tm_penalties

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