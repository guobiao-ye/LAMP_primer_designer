from Bio import SeqIO
from utils import oligos
from utils import runprimer3
from utils import penalty
from utils.config_loader import load_config
from core import overlap
from core import combination
from core import mutant

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

def find_primers(fasta_file):
    sequence = None
    quasispecies_sequences = []
    outer_primer_params = load_config("configs/outer_primer.json")
    middle_primer_params = load_config("configs/middle_primer.json")
    loop_primer_params = load_config("configs/loop_primer.json")
    inner_primer_params = load_config("configs/inner_primer.json")

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
        forward_inner_c_oligos = oligos.extract_oligos(runprimer3.run_primer3(sequence, inner_primer_params, reverse=True),
                                                sequence_length, is_reverse=True)
        reverse_inner_c_oligos = oligos.extract_oligos(runprimer3.run_primer3(sequence, inner_primer_params, reverse=False),
                                                sequence_length, is_reverse=False)
        forward_inner_oligos = oligos.move_position(forward_inner_c_oligos, 'end_to_start')
        reverse_inner_oligos = oligos.move_position(reverse_inner_c_oligos, 'start_to_end')
        forward_loop_oligos = oligos.extract_oligos(runprimer3.run_primer3(sequence, loop_primer_params, reverse=True), sequence_length,
                                             is_reverse=True)
        reverse_loop_oligos = oligos.extract_oligos(runprimer3.run_primer3(sequence, loop_primer_params, reverse=False), sequence_length,
                                             is_reverse=False)
        forward_middle_oligos = oligos.extract_oligos(runprimer3.run_primer3(sequence, middle_primer_params, reverse=False),
                                               sequence_length, is_reverse=False)
        reverse_middle_oligos = oligos.extract_oligos(runprimer3.run_primer3(sequence, middle_primer_params, reverse=True),
                                               sequence_length, is_reverse=True)
        forward_outer_oligos = oligos.extract_oligos(runprimer3.run_primer3(sequence, outer_primer_params, reverse=False),
                                              sequence_length, is_reverse=False)
        reverse_outer_oligos = oligos.extract_oligos(runprimer3.run_primer3(sequence, outer_primer_params, reverse=True), sequence_length,
                                              is_reverse=True)
        
        # This chunk of code add penalty for mutant
        if quasispecies_sequences:
            mutation_list = mutant.generate_mutation_list(sequence, quasispecies_sequences)
            forward_inner_oligos = mutant.filter_forward_oligos(forward_inner_oligos, sequence, mutation_list)
            reverse_inner_oligos = mutant.filter_reverse_oligos(reverse_inner_oligos, sequence, mutation_list)
            forward_loop_oligos = mutant.filter_reverse_oligos(forward_loop_oligos, sequence, mutation_list)
            reverse_loop_oligos = mutant.filter_forward_oligos(reverse_loop_oligos, sequence, mutation_list)
            forward_middle_oligos = mutant.filter_forward_oligos(forward_middle_oligos, sequence, mutation_list)
            reverse_middle_oligos = mutant.filter_reverse_oligos(reverse_middle_oligos, sequence, mutation_list)
            forward_outer_oligos = mutant.filter_forward_oligos(forward_outer_oligos, sequence, mutation_list)
            reverse_outer_oligos = mutant.filter_reverse_oligos(reverse_outer_oligos, sequence, mutation_list)

        print("All oligos found!")
        print(
            f"Number of oligos are: forward_inner - {len(forward_inner_oligos)}, reverse_inner - {len(reverse_inner_oligos)}, forward_loop - {len(forward_loop_oligos)}, reverse_loop - {len(reverse_loop_oligos)}, forward_middle - {len(forward_middle_oligos)}, reverse_middle - {len(reverse_middle_oligos)}, forward_outer - {len(forward_outer_oligos)}, reverse_outer - {len(reverse_outer_oligos)}")
        
        # This chunk of code reduces overlap
        forward_inner_candidates = overlap.reduce_primers_by_overlap(forward_inner_oligos, 60, False)
        reverse_inner_candidates = overlap.reduce_primers_by_overlap(reverse_inner_oligos, 60, True)
        forward_loop_candidates = overlap.reduce_primers_by_overlap(forward_loop_oligos, 60, True)
        reverse_loop_candidates = overlap.reduce_primers_by_overlap(reverse_loop_oligos, 60, False)
        forward_middle_candidates = overlap.reduce_primers_by_overlap(forward_middle_oligos, 60, False)
        reverse_middle_candidates = overlap.reduce_primers_by_overlap(reverse_middle_oligos, 60, True)
        forward_outer_candidates = overlap.reduce_primers_by_overlap(forward_outer_oligos, 60, False)
        reverse_outer_candidates = overlap.reduce_primers_by_overlap(reverse_outer_oligos, 60, True)

        print("All candidates prepared by reducing overlaps!")
        print(
            f"Number of candidates are:  forward_inner - {len(forward_inner_candidates)}, reverse_inner - {len(reverse_inner_candidates)}, forward_loop - {len(forward_loop_candidates)}, reverse_loop - {len(reverse_loop_candidates)}, forward_middle - {len(forward_middle_candidates)}, reverse_middle - {len(reverse_middle_candidates)}, forward_outer - {len(forward_outer_candidates)}, reverse_outer - {len(reverse_outer_candidates)}")

        distance_penalties = penalty.generate_distance_penalties(max_distance=320)
        best_forward_infos, best_forward_penalties, forward_set_count = combination.find_best_forward_combinations(
            forward_inner_candidates,
            forward_loop_candidates,
            forward_middle_candidates,
            forward_outer_candidates,
            320, 1, 25, True, 0, distance_penalties,
            [1, 0.6, 1.2, 1], [1, 1, 1])
        print(f'{forward_set_count} sets of forward combinations are generated.')

        best_reverse_infos, best_reverse_penalties, reverse_set_count = combination.find_best_reverse_combinations(
            reverse_inner_candidates,
            reverse_loop_candidates,
            reverse_middle_candidates,
            reverse_outer_candidates,
            320, 1, 25, True, 0, distance_penalties,
            [1, 0.6, 1.2, 1], [1, 1, 1])
        print(f'{reverse_set_count} sets of reverse combinations are generated.')

        combination.find_possible_result(forward_inner_candidates, reverse_inner_candidates, best_forward_infos, best_reverse_infos,
                             best_forward_penalties, best_reverse_penalties, 320, 1, 4, 0.6, distance_penalties, True,
                             "./data/output/Intermediate_file/candidate_result.csv")

find_primers('data/input/example.fasta')
