import time
import numpy as np
from scipy.special import logsumexp
from Seq import *
from Bio.SeqIO import parse
import os
import sys

FASTQ_PATH = '../fastq'

NOT_MOTIF_STATES = 5

STATES_AT_START = 3

FIRST_MOTIF_STATE = 3

LINE_LENGTH = 50

START_SIGN = '^'
END_SIGN = '$'

emission_dict = {START_SIGN: 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4, END_SIGN: 5}

INITIAL_STATE_INDEX = 3

ADDITIONAL_LETTERS = 2

NUM_OF_BASES = 4

base_dict = {START_SIGN: 0, "A": 1, "C": 2, "G": 3, "T": 4, END_SIGN: 5}

DNA_PROBABILITY = [0.22, 0.28, 0.28, 0.22]


def initial_emissions_mat_calc(seeds, alpha):
    """
    calc initial emissions mat
    :param seed: motif
    :param alpha: Set the initial emission probabilities for the motif states to have 1-3 * alpha for the motif letter,
    and alpha for others letters.
    :return: numpy array
    """
    n = 0
    for seed in seeds:
        n += len(seed)
    emissions_mat = np.zeros((NUM_OF_BASES + ADDITIONAL_LETTERS, n + NOT_MOTIF_STATES))
    # B start
    emissions_mat[0, 0] = 1
    # B end
    emissions_mat[-1, -1] = 1
    # B_pre_telo
    emissions_mat[1:-1, 1] = DNA_PROBABILITY
    # B_norm
    emissions_mat[1:-1, -2] = DNA_PROBABILITY
    # B_telo
    emissions_mat[1:-1, 2] = DNA_PROBABILITY
    seed = ''.join(seeds)
    for column, letter in enumerate(seed):
        emissions_mat[1:-1, column + INITIAL_STATE_INDEX] = alpha
        emissions_mat[
            base_dict[letter], column + INITIAL_STATE_INDEX] = 1 - 3 * alpha
    return wrap_log(emissions_mat).T, n


def build_transition_matrix(seeds, telo_in_seq, p_enter_telo_from_background, p_exit_telo, p_end_of_seq,
                            prob_for_primary_motif,
                            p_same_motif_block,
                            p_exit_from_motif_to_backgroud, p_telo_background_to_motif):
    """
    build transition matrix according to number of states and given q,p
    :param k_counter: num of states
    :param q: probability of motif in seq
    :param p: probability to enter motif
    :return: transition matrix
    """
    # chance to enter motif from telomere background
    p_telo_no_exit = 1 - p_exit_telo
    p_telo_background = p_telo_no_exit * (1 - p_telo_background_to_motif)
    p_primary_motif = p_telo_no_exit * p_telo_background_to_motif * prob_for_primary_motif
    p_secondary_motif = p_telo_no_exit * p_telo_background_to_motif * (1 - prob_for_primary_motif) / 2
    k_counter = 0
    for seed in seeds:
        k_counter += len(seed)
    dim = k_counter + NOT_MOTIF_STATES
    transition_mat = np.zeros([dim, dim])
    # seq contain contain_telo
    transition_mat[0, 1] = telo_in_seq
    transition_mat[0, -2] = 1 - telo_in_seq
    # chance to enter contain_telo
    transition_mat[1, 1] = 1 - p_enter_telo_from_background
    transition_mat[1, 2] = p_enter_telo_from_background
    # B_telo
    transition_mat[2, 2] = p_telo_background
    transition_mat[2, -2] = p_exit_telo
    # transition_mat[2, -1] = p_end_of_seq
    # B_norm
    transition_mat[-2, -2] = 1
    # transition_mat[-2, -1] = p_end_of_seq
    previous_length = 0
    for i in range(len(seeds)):
        # np_eye
        transition_mat[STATES_AT_START + previous_length:STATES_AT_START + previous_length + len(seeds[i]) - 1,
        STATES_AT_START + previous_length + 1:STATES_AT_START + previous_length + len(seeds[i])] = np.eye(
            len(seeds[i]) - 1)
        # motif_last - back to B_telo
        transition_mat[STATES_AT_START + previous_length + len(seeds[i]) - 1, 2] = p_exit_from_motif_to_backgroud
        # one more copy of same motif
        transition_mat[STATES_AT_START + previous_length + len(
            seeds[i]) - 1, STATES_AT_START + previous_length] = p_same_motif_block
        prev_prev_len = 0
        for j in range(len(seeds)):
            if j != i:
                # prob to switch motifs immediately
                transition_mat[
                    STATES_AT_START + previous_length + len(seeds[i]) - 1, STATES_AT_START + prev_prev_len] = \
                    (1 - p_same_motif_block - p_exit_from_motif_to_backgroud) / (len(seeds) - 1)
            prev_prev_len += len(seeds[j])
        # B_telo to motif
        if i == 0:
            transition_mat[2, STATES_AT_START + previous_length] = p_primary_motif
        else:
            transition_mat[2, STATES_AT_START + previous_length] = p_secondary_motif
        previous_length += len(seeds[i])
    transition_mat = transition_mat * (1 - p_end_of_seq)
    transition_mat[:, -1] = p_end_of_seq
    transition_mat[-1, -1] = 1
    for row in transition_mat:
        if np.around(np.sum(row), decimals=5) != 1:
            raise Exception()
    return wrap_log(transition_mat)


def forward(seq, emission_mat, transition_mat, k_counter):
    """
    calculate the forward table for a given seq
    :param seq: sequence
    :param emission_mat
    :param transition_mat
    :param k_counter: number of states
    :return: Forward table
    """
    k_dim = k_counter + NOT_MOTIF_STATES
    N = len(seq)
    forward_table = wrap_log(np.zeros([k_dim, N]))
    forward_table[0, 0] = wrap_log(1)
    for j in range(1, N):
        curr_letter = forward_table[:, j - 1].reshape(-1, 1)
        forward_table[:, j] = logsumexp(curr_letter + transition_mat, axis=0) + emission_mat[:, emission_dict[seq[j]]]
    return forward_table


def backward(seq, emission_mat, transition_mat, k_counter):
    """
    calculate the backward table for a given seq
    :param seq: sequence
    :param emission_mat
    :param transition_mat
    :param k_counter: number of states
    :return: Backward table
    """
    k_dim = k_counter + NOT_MOTIF_STATES
    N = len(seq)
    backward_table = wrap_log(np.zeros([k_dim, N]))
    backward_table[-1, -1] = wrap_log(1)
    for j in range(N - 2, -1, -1):
        curr_letter = backward_table[:, j + 1].reshape(-1, 1)
        backward_table[:, j] = logsumexp(
            curr_letter + transition_mat.T + emission_mat[:, emission_dict[seq[j + 1]]].reshape((-1, 1)), axis=0)
    return backward_table


def posterior(seq, emission_mat, transition_mat, k_counter, seeds, rec_num, counter):
    """
    calculates the most probable state for every base in seq
    :param seq: sequence
    :param emission_mat
    :param transition_mat
    :param k_counter: num of states
    :return: seq of states, aligned to original seq
    """
    N = len(seq)
    forward_table = forward(seq, emission_mat, transition_mat, k_counter)
    backward_table = backward(seq, emission_mat, transition_mat, k_counter)
    posterior_table = forward_table + backward_table
    # motif_order = EMPTY_STRING

    seq_obj = Seq(N - 2, rec_num)

    last_motif_0 = FIRST_MOTIF_STATE + len(seeds[0]) - 1
    first_motif_1 = last_motif_0 + 1
    last_motif_1 = first_motif_1 + len(seeds[1]) - 1
    first_motif_2 = last_motif_1 + 1
    last_motif_2 = first_motif_2 + len(seeds[2]) - 1

    # decide states
    for j in range(1, N - 1):
        curr_k = int(np.argmax(posterior_table[:, j]))

        if FIRST_MOTIF_STATE <= curr_k <= last_motif_0:
            # motif_order += MOTIF_0
            seq_obj.add_motif_base(0, (seq[j], curr_k - FIRST_MOTIF_STATE), j - 1)

        elif first_motif_1 <= curr_k <= last_motif_1:
            # motif_order += MOTIF_1
            seq_obj.add_motif_base(1, (seq[j], curr_k - first_motif_1), j - 1)

        elif first_motif_2 <= curr_k <= last_motif_2:
            # motif_order += MOTIF_2
            seq_obj.add_motif_base(2, (seq[j], curr_k - first_motif_2), j - 1)

        elif curr_k == 2:
            # motif_order += TELO_BACKGROUND
            seq_obj.add_telo_background(seq[j], j - 1)
        elif curr_k == 1:
            # motif_order += 'P'
            seq_obj.add_pre_telo((seq[j], curr_k))
        else:
            # motif_order += BACKGROUND
            seq_obj.add_normal_dna_base((seq[j], curr_k))

    # print_results(seq[1:-1], motif_order)
    seq_obj.print_statistics(doc=None, counter=counter)

    seq_obj.save_to_file()

    return


# def viterbi(seq, emission_mat, transition_mat, k_counter, seeds):
#     """
#     calculates the most probable motif location
#     :param seq: sequence
#     :param emission_mat
#     :param transition_mat
#     :param k_counter: num of states
#     :return: seq of states, aligned to original seq
#     """
#     k_dim = k_counter + NOT_MOTIF_STATES
#     N = len(seq)
#     prob_mat = wrap_log(np.zeros([k_dim, N]))
#     trace_mat = np.zeros([k_dim, N])
#     prob_mat[0, 0] = wrap_log(1)
#     for j in range(1, N):
#         curr_letter = prob_mat[:, j - 1].reshape((-1, 1))
#         potential_trans = curr_letter + transition_mat
#         max_values = np.max(potential_trans, axis=0).T
#         trace_mat[:, j] = np.argmax(potential_trans, axis=0).T
#         prob_mat[:, j] = max_values + emission_mat[:, emission_dict[seq[j]]]
#     # begin trace
#     motif_order = EMPTY_STRING
#     curr_k = int(np.argmax(prob_mat[:, -1]))
#     for j in range(N - 1, -1, -1):
#         last_motif_0 = FIRST_MOTIF_STATE + len(seeds[0])
#         first_motif_1 = last_motif_0 + 1
#         last_motif_1 = first_motif_1 + len(seeds[1])
#         first_motif_2 = last_motif_1 + 1
#         last_motif_2 = first_motif_2 + len(seeds[2])
#         if FIRST_MOTIF_STATE <= curr_k <= last_motif_0:
#             motif_order = MOTIF_0 + motif_order
#         elif first_motif_1 <= curr_k <= last_motif_1:
#             motif_order = MOTIF_1 + motif_order
#         elif first_motif_2 <= curr_k <= last_motif_2:
#             motif_order = MOTIF_2 + motif_order
#         elif curr_k == 2:
#             motif_order = TELO_BACKGROUND + motif_order
#         else:
#             motif_order = BACKGROUND + motif_order
#         curr_k = int(trace_mat[curr_k, j])
#     return motif_order[1:-1]


def wrap_log(to_wrap):
    """
    helper func to avoid log warnings for np.log
    :param to_wrap: element to log
    :return: np.log(to_wrap)
    """
    with np.errstate(divide='ignore'):
        result = np.log(to_wrap)
    return result


def parse_seqs(fasta_path):
    sequences = []
    file = open(fasta_path, 'r')
    try:
        parsed_seqs = parse(file, format='fastq')
        for seq in parsed_seqs:
            sequences.append((seq.id, START_SIGN + seq.seq._data + END_SIGN))
    except ValueError as err:
        msg = 'File '+fasta_path + ' - ' + str(err)
        print("\033[1;31m" + msg + '\033[0m')
    finally:
        file.close()
        return sequences


def main(fastq_files_path):
    # M W D
    seeds = ['TTAGGG', 'TTAAAA', 'TTGGGG']

    emission_mat, k_counter = initial_emissions_mat_calc(seeds, alpha=(0.08 / 3))
    transition_mat = build_transition_matrix(seeds, telo_in_seq=0.0005, p_enter_telo_from_background=0.45,
                                             p_exit_telo=0.0002,
                                             p_end_of_seq=0.00005, prob_for_primary_motif=0.86,
                                             p_same_motif_block=0.75,
                                             p_exit_from_motif_to_backgroud=0.15, p_telo_background_to_motif=0.7)

    files_visited = []

    for filename in os.listdir(fastq_files_path):
        file_num = int(filename.split("_")[1].split(".")[0])
        if not (int(sys.argv[1]) <= file_num < int(sys.argv[2])):
            continue
        path = os.path.join(fastq_files_path, filename)
        print(path)
        seqs = parse_seqs(path)
        if not seqs:
            continue
        files_visited.append(filename)
        with open(sys.argv[1] + "_" + sys.argv[2] + '.txt', 'w') as file_handler:
            for item in files_visited:
                file_handler.write("%s\n" % item)
        print('len' + str(len(seqs)))
        for value, seq in enumerate(seqs):
            posterior(seq[1], emission_mat, transition_mat, k_counter, seeds, seq[0], value)


if __name__ == '__main__':
    start = time.time()
    main(FASTQ_PATH)
    end = time.time()
    print("time:" + str(end - start))
