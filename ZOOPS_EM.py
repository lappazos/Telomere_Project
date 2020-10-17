import argparse
import HMM as mf
import numpy as np
from Bio import SeqIO
from scipy.special import logsumexp
import os
import pickle

INITIAL_LETTER_INDEX = 1

INITIAL_STATE_INDEX = 2

ADDITIONAL_LETTERS = 2

NON_MOTIF_STATES = 4

NUM_OF_BASES = 4

base_dict = {mf.START_SIGN: 0, "A": 1, "C": 2, "G": 3, "T": 4, mf.END_SIGN: 5}

alphabet = [mf.START_SIGN, 'A', 'C', 'G', 'T', mf.END_SIGN]


# from stackoverflow
def log_space_product(A, B):
    Astack = np.stack([A] * A.shape[0]).transpose(2, 1, 0)
    Bstack = np.stack([B] * B.shape[1]).transpose(1, 0, 2)
    return logsumexp(Astack + Bstack, axis=0)


def initial_emissions_mat_calc(seed, alpha):
    """
    calc initial emissions mat
    :param seed: motif
    :param alpha: Set the initial emission probabilities for the motif states to have 1-3 * alpha for the motif letter,
    and alpha for others letters.
    :return: numpy array
    """
    n = len(seed)
    emissions_mat = np.zeros((NUM_OF_BASES + ADDITIONAL_LETTERS, n + NON_MOTIF_STATES))
    # B start
    emissions_mat[0, 0] = 1
    # B end
    emissions_mat[-1, -1] = 1
    # B_1
    emissions_mat[1:-1, 1] = mf.UNIFORM_PROB
    # B_2
    emissions_mat[1:-1, -2] = mf.UNIFORM_PROB
    for column, letter in enumerate(seed):
        emissions_mat[1:-1, column + INITIAL_STATE_INDEX] = alpha
        emissions_mat[
            base_dict[letter], column + INITIAL_STATE_INDEX] = 1 - 3 * alpha
    return mf.wrap_log(emissions_mat).T


def parse_seqs(fasta_path):
    """
    :param fasta_path: path of fasta file
    :return: array of seqs
    """
    sequences = []
    for seq in SeqIO.parse(fasta_path, format='fasta'):
        sequences.append(mf.START_SIGN + seq.seq._data + mf.END_SIGN)
    return sequences


def calc_q(sequences, seed):
    """
    calculate initial q
    :param sequences: the data provided - dna sequences
    :param seed: motif to search
    :return: q probability
    """
    counter = 0
    for s in sequences:
        if seed in s:
            counter += 1
    return counter / len(sequences)


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
    np.savetxt('motif_profile.txt', emissions_mat[1:-1, 2:-2], fmt='%.2f', delimiter='\t')
    f = open('motif_profile.txt', 'a')
    p_enter_telo_from_background, p_exit_telo, prob_for_primary_motif, p_same_motif_block, p_exit_from_motif_to_backgroud, p_telo_background_to_motif = np.exp(tuple)
    f.write("p_enter_telo_from_background" + str("%.2f" % p_enter_telo_from_background) + "\n")
    f.write("p_exit_telo" + str("%.2f" % p_exit_telo) + "\n")
    f.write("prob_for_primary_motif" + str("%.2f" % prob_for_primary_motif) + "\n")
    f.write("p_same_motif_block" + str("%.2f" % p_same_motif_block) + "\n")
    f.write("p_exit_from_motif_to_backgroud" + str("%.2f" % p_exit_from_motif_to_backgroud) + "\n")
    f.write("p_telo_background_to_motif" + str("%.2f" % p_telo_background_to_motif) + "\n")
    f.close()


def write_motif_positions(viterbi_output_list):
    """
    prints to file the index of every motif per sequence. -1 if there isn't any
    :param viterbi_output_list: list of viterbi motif sequences
    """
    f = open('motif_positions.txt', 'w')
    for i, seq in enumerate(viterbi_output_list):
        if i != len(viterbi_output_list) - 1:
            f.write(str(seq.find(mf.MOTIF)) + "\n")
        else:
            f.write(str(seq.find(mf.MOTIF)))
    f.close()


def EM(emissions_mat, transitions_mat, sequences, threshold, k_counter,seeds):
    """
    performs EM algorithm to learn emission and transition matrices
    :param emissions_mat: initial emission mat
    :param transitions_mat: initial transition mat
    :param sequences: sequences to learn from
    :param threshold: convergences threshold
    :param k_counter: number of unique motifs
    :return: final emission mat, final transition mat, Log likelihood history
    """
    dim = k_counter + mf.NOT_MOTIF_STATES
    diff = np.NINF
    ll_history = []
    while True:
        curr_diff = 0

        # holds matirx of each letter seperately for N_kx
        positions_dict = {"A": np.array([]).reshape(dim, -1),
                          "C": np.array([]).reshape(dim, -1),
                          "G": np.array([]).reshape(dim, -1),
                          "T": np.array([]).reshape(dim, -1),
                          mf.START_SIGN: np.array([]).reshape(dim, -1),
                          mf.END_SIGN: np.array([]).reshape(dim, -1)}

        N_kl = mf.wrap_log(np.zeros((dim, dim)))
        counter = 0
        length = len(sequences)
        for seq in sequences:
            # print(counter/length)
            counter+=1
            F = mf.forward(seq, emissions_mat, transitions_mat, k_counter)
            seq_likelihood = F[-1][-1]
            curr_diff += seq_likelihood
            B = mf.backward(seq, emissions_mat, transitions_mat,
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
            padding_line = mf.wrap_log(np.zeros((dim, 1)))
            F = np.concatenate((padding_line, F), axis=1)
            B = np.concatenate((B, padding_line), axis=1).T
            N_kl = np.logaddexp(N_kl, log_space_product(F, B))

        ll_history.append(curr_diff)
        print(curr_diff)

        if curr_diff - diff <= threshold:
            return emissions_mat, (p_enter_telo_from_background,p_exit_telo,prob_for_primary_motif,p_same_motif_block,
                                                 p_exit_from_motif_to_backgroud, p_telo_background_to_motif), ll_history
        else:
            diff = curr_diff

        N_kx = np.empty((dim, 0))
        for letter in alphabet:
            N_kx = np.concatenate((N_kx, logsumexp(positions_dict[letter], axis=1).reshape(-1, 1)), axis=1)

        N_kl += transitions_mat

        N_kx_relevant_cols = N_kx[:, 1:-1]
        emissions_mat[2:21, 1:-1] = (N_kx_relevant_cols - logsumexp(N_kx_relevant_cols, axis=1).reshape(-1, 1))[2:21, :]
        p_enter_telo_from_background = N_kl[1,2]-logsumexp(N_kl[1, :])
        logsum2=logsumexp(N_kl[2, :])
        p_exit_telo = N_kl[2,-2]-logsum2
        p_telo_no_exit = np.log(-np.expm1(p_exit_telo))
        p_telo_background = N_kl[2,2]-logsum2
        p_telo_background_to_motif= np.log(-np.expm1(p_telo_background - p_telo_no_exit))
        p_primary_motif = N_kl[2,3]-logsum2
        prob_for_primary_motif = p_primary_motif - (p_telo_no_exit + p_telo_background_to_motif)
        logsum8 = logsumexp(N_kl[8, :])
        p_same_motif_block = N_kl[8,3]-logsum8
        p_exit_from_motif_to_backgroud = N_kl[8,2]-logsum8
        transitions_mat = mf.build_transition_matrix(seeds, telo_in_seq=0.0005, p_enter_telo_from_background=np.exp(p_enter_telo_from_background),
                                                 p_exit_telo=np.exp(p_exit_telo),
                                                 p_end_of_seq=0.00005, prob_for_primary_motif=np.exp(prob_for_primary_motif),
                                                 p_same_motif_block=np.exp(p_same_motif_block),
                                                 p_exit_from_motif_to_backgroud=np.exp(p_exit_from_motif_to_backgroud), p_telo_background_to_motif=np.exp(p_telo_background_to_motif))


def parse_args():
    """
    parse arguments
    :return: parsed args
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='File path with list of sequences (e.g. yeastGenes.fasta)')
    parser.add_argument('seed', help='Guess for the motif (e.g. ATTA)')
    parser.add_argument('p', type=float,
                        help='Initial guess for the p transition probability (e.g. 0.01)')
    parser.add_argument('alpha', type=float,
                        help='Softening parameter for the initial profile (e.g. 0.1)')
    parser.add_argument('convergenceThr', type=float,
                        help='ll improvement threshold for the stopping condition'
                             ' (e.g. 0.1)')
    return parser.parse_args()


def main():
    """
    main function
    """
    np.seterr(all='warn', invalid='ignore')

    seeds = ['TTAGGG', 'TTAAAA', 'TTGGGG']

    emission_mat, k_counter = mf.initial_emissions_mat_calc(seeds, alpha=(0.08 / 3))
    transition_mat = mf.build_transition_matrix(seeds, telo_in_seq=0.0005, p_enter_telo_from_background=0.45,
                                             p_exit_telo=0.0002,
                                             p_end_of_seq=0.00005, prob_for_primary_motif=0.86,
                                             p_same_motif_block=0.75,
                                             p_exit_from_motif_to_backgroud=0.15, p_telo_background_to_motif=0.7)

    seqs = []
    for filename in os.listdir('join'):
        file = open(os.path.join('join', filename), 'rb')
        telo = pickle.load(file)
        seq, _ = telo.generate_seq()
        seqs.append('^' + seq + '$')
    emissions_mat, tuple, ll_history = EM(emission_mat, transition_mat, seqs, 2,
                                                    k_counter,seeds)
    write_ll_history(ll_history)
    write_motif_profile(np.exp(emissions_mat.T),tuple)


if __name__ == "__main__":
    main()
