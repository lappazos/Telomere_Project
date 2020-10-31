import Telo_Detecter as td
import numpy as np
from scipy.special import logsumexp
import os
import pickle
from Seq import Telomere

TELO_PATH = '.'

INITIAL_LETTER_INDEX = 1

INITIAL_STATE_INDEX = 2

ADDITIONAL_LETTERS = 2

NON_MOTIF_STATES = 4

NUM_OF_BASES = 4

base_dict = {td.START_SIGN: 0, "A": 1, "C": 2, "G": 3, "T": 4, td.END_SIGN: 5}

alphabet = [td.START_SIGN, 'A', 'C', 'G', 'T', td.END_SIGN]


# from stackoverflow
def log_space_product(A, B):
    Astack = np.stack([A] * A.shape[0]).transpose(2, 1, 0)
    Bstack = np.stack([B] * B.shape[1]).transpose(1, 0, 2)
    return logsumexp(Astack + Bstack, axis=0)


def write_ll_history(ll_iterations):
    """
    writes log likelihood to file
    :param ll_iterations: array of log likelihood per iteration
    """
    f = open("ll_history.txt", 'w')
    for i, ll in enumerate(ll_iterations):
        if i != len(ll_iterations) - 1:
            f.write(str(ll) + "\n")
        else:
            f.write(str(ll))
    f.close()


def write_motif_profile(emissions_mat, tuple):
    """
    prints to file according to agreed format the final Emission and transition results
    :param emissions_mat: final emission matrix
    :param q: final q
    :param p: final p
    """
    np.savetxt('motif_profile.txt', emissions_mat[1:-1, 1:-1], fmt='%.2f', delimiter='\t')
    f = open('motif_profile.txt', 'a')
    p_enter_telo_from_background, p_exit_telo, prob_for_primary_motif, p_same_motif_block, p_exit_from_motif_to_backgroud, p_telo_background_to_motif = np.exp(
        tuple)
    f.write("p_enter_telo_from_background" + str("%.2f" % p_enter_telo_from_background) + "\n")
    f.write("p_exit_telo" + str("%.2f" % p_exit_telo) + "\n")
    f.write("prob_for_primary_motif" + str("%.2f" % prob_for_primary_motif) + "\n")
    f.write("p_same_motif_block" + str("%.2f" % p_same_motif_block) + "\n")
    f.write("p_exit_from_motif_to_backgroud" + str("%.2f" % p_exit_from_motif_to_backgroud) + "\n")
    f.write("p_telo_background_to_motif" + str("%.2f" % p_telo_background_to_motif) + "\n")
    f.close()


def EM(emissions_mat, transitions_mat, sequences, threshold, k_counter, seeds):
    """
    performs EM algorithm to learn emission and transition matrices
    :param emissions_mat: initial emission mat
    :param transitions_mat: initial transition mat
    :param sequences: sequences to learn from
    :param threshold: convergences threshold
    :param k_counter: number of unique motifs
    :return: final emission mat, final transition mat, Log likelihood history
    """
    dim = k_counter + td.NOT_MOTIF_STATES
    diff = np.NINF
    ll_history = []
    first_iter = True
    while True:
        curr_diff = 0

        # holds matirx of each letter seperately for N_kx
        positions_dict = {"A": np.array([]).reshape(dim, -1),
                          "C": np.array([]).reshape(dim, -1),
                          "G": np.array([]).reshape(dim, -1),
                          "T": np.array([]).reshape(dim, -1),
                          td.START_SIGN: np.array([]).reshape(dim, -1),
                          td.END_SIGN: np.array([]).reshape(dim, -1)}

        N_kl = td.wrap_log(np.zeros((dim, dim)))
        counter = 0
        length = len(sequences)
        for seq in sequences:
            # print(counter/length)
            counter += 1
            F = td.forward(seq, emissions_mat, transitions_mat, k_counter)
            seq_likelihood = F[-1][-1]
            curr_diff += seq_likelihood
            B = td.backward(seq, emissions_mat, transitions_mat,
                            k_counter)
            B -= seq_likelihood  # -seq is instaed of /p(x_j), relevant for both n_xx n_kl
            P = F + B

            # for N_kx

            seq_to_list = np.array(list(seq))
            for letter in positions_dict:
                letter_indices = np.argwhere(seq_to_list == letter).T

                # for N_kx
                positions_dict[letter] = np.concatenate(
                    (positions_dict[letter], np.take(P, letter_indices, axis=1).reshape(dim, -1)), axis=1)

                # for N_kl - mult by relevant emission
                B[:, letter_indices.reshape(-1, )] += emissions_mat[:, base_dict[letter]].reshape(
                    -1, 1)

            # for N_kl
            padding_line = td.wrap_log(np.zeros((dim, 1)))
            F = np.concatenate((padding_line, F), axis=1)
            B = np.concatenate((B, padding_line), axis=1).T
            N_kl = np.logaddexp(N_kl, log_space_product(F, B))

        ll_history.append(curr_diff)
        print(curr_diff)

        if curr_diff - diff <= threshold:
            return emissions_mat, data, ll_history
        elif not first_iter:
            write_ll_history(ll_history)
            write_motif_profile(np.exp(emissions_mat.T), data)
            diff = curr_diff
        else:
            pass

        N_kx = np.empty((dim, 0))
        for letter in alphabet:
            N_kx = np.concatenate((N_kx, logsumexp(positions_dict[letter], axis=1).reshape(-1, 1)), axis=1)

        N_kl += transitions_mat

        N_kx_relevant_cols = N_kx[:, 1:-1]
        emissions_mat[2:21, 1:-1] = (N_kx_relevant_cols - logsumexp(N_kx_relevant_cols, axis=1).reshape(-1, 1))[2:21, :]
        p_enter_telo_from_background = N_kl[1, 2] - logsumexp(N_kl[1, :])
        logsum2 = logsumexp(N_kl[2, :])
        p_exit_telo = N_kl[2, -2] - logsum2
        p_telo_no_exit = np.log(-np.expm1(p_exit_telo))
        p_telo_background = N_kl[2, 2] - logsum2
        p_telo_background_to_motif = np.log(-np.expm1(p_telo_background - p_telo_no_exit))
        p_primary_motif = N_kl[2, 3] - logsum2
        prob_for_primary_motif = p_primary_motif - (p_telo_no_exit + p_telo_background_to_motif)
        logsum8 = logsumexp(N_kl[8, :])
        p_same_motif_block = N_kl[8, 3] - logsum8
        p_exit_from_motif_to_backgroud = N_kl[8, 2] - logsum8
        transitions_mat = td.build_transition_matrix(seeds, telo_in_seq=0.0005,
                                                     p_enter_telo_from_background=np.exp(p_enter_telo_from_background),
                                                     p_exit_telo=np.exp(p_exit_telo),
                                                     p_end_of_seq=0.00005,
                                                     prob_for_primary_motif=np.exp(prob_for_primary_motif),
                                                     p_same_motif_block=np.exp(p_same_motif_block),
                                                     p_exit_from_motif_to_backgroud=np.exp(
                                                         p_exit_from_motif_to_backgroud),
                                                     p_telo_background_to_motif=np.exp(p_telo_background_to_motif))
        data = (
            p_enter_telo_from_background, p_exit_telo, prob_for_primary_motif, p_same_motif_block,
            p_exit_from_motif_to_backgroud, p_telo_background_to_motif)
        first_iter = False


def main(path):
    """
    main function
    """
    np.seterr(all='warn', invalid='ignore')

    seeds = ['TTAGGG', 'TTAAAA', 'TTGGGG']

    emission_mat, k_counter = td.initial_emissions_mat_calc(seeds, alpha=(0.08 / 3))
    transition_mat = td.build_transition_matrix(seeds, telo_in_seq=0.0005, p_enter_telo_from_background=0.1,
                                                p_exit_telo=0.0002,
                                                p_end_of_seq=0.00005, prob_for_primary_motif=0.8,
                                                p_same_motif_block=0.65,
                                                p_exit_from_motif_to_backgroud=0.25, p_telo_background_to_motif=0.6)

    seqs = []
    parse_dir('HH', seqs)
    parse_dir('Normal', seqs)
    emissions_mat, tuple, ll_history = EM(emission_mat, transition_mat, seqs, 3,
                                          k_counter, seeds)
    write_ll_history(ll_history)
    write_motif_profile(np.exp(emissions_mat.T), tuple)


def parse_dir(path, seqs):
    for folder in os.listdir(path):
        if os.path.isdir(os.path.join(path, folder, 'cluster')):
            for filename in os.listdir(os.path.join(path, folder, 'cluster')):
                add = True
                _, file_extension = os.path.splitext(filename)
                if file_extension != '.obj':
                    continue
                file = open(os.path.join(path, folder, 'cluster', filename), 'rb')
                telo = pickle.load(file)
                file.close()
                telo.generate_seq()
                if telo.more_then_one_telo:
                    continue
                for elem in telo.seq:
                    if isinstance(elem, Telomere):
                        if elem.num_of_motifs <= 25:
                            add = False
                        if (25 < elem.num_of_motifs < 275) and (elem.motif_types_num[2] > 0.79):
                            add = False
                if add:
                    seqs.append('^' + telo.sequence + '$')


if __name__ == "__main__":
    main(TELO_PATH)
