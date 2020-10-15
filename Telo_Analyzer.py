import pickle
import Seq
import os
import matplotlib.pyplot as plt

MIN_SEQ_LEN_TO_PRINT = 1500

MIN_SUBTELO_LENGTH = 50

MIN_TELO_LENGTH = 500


def analyze_telos(path):
    telos_lengths = []
    subtelos_length = []
    motif_dic = {}
    for filename in os.listdir(path):
        print(filename)
        file = open(os.path.join(path, filename), 'rb')
        telo = pickle.load(file)
        file.close()
        if not isinstance(telo, Seq.Seq):
            continue
        if telo.longest_telomere_len > MIN_TELO_LENGTH:
            telos_lengths.append(telo.longest_telomere_len)
            if telo.non_telomeric_parts:
                subtelo_len = max(telo.non_telomeric_parts)
                if subtelo_len > MIN_SUBTELO_LENGTH:
                    subtelos_length.append(max(telo.non_telomeric_parts))
        if len(telo.len) > MIN_SEQ_LEN_TO_PRINT:
            telo.print_doc()
    plt.hist(telos_lengths, bins=50, color='orange')
    # plt.xscale('log')
    plt.title('Telomere length')
    plt.xlabel('# Base-Pairs')
    plt.ylabel('# Reads')
    plt.show()
    print(subtelos_length)
    plt.hist(subtelos_length, bins=50, color='orange')
    plt.title('Sub-Telomere length')
    plt.xlabel('# Base-Pairs')
    # plt.xscale('log')
    plt.ylabel('# Reads')
    plt.show()
    d_view = [(v, k) for k, v in motif_dic.items()]
    d_view.sort(reverse=True)  # natively sort tuples by first element
    val_tot = 0
    for v, k in d_view:
        val_tot += v
    for v, k in d_view:
        print("%s: %d %d" % (k, v, (v / val_tot) * 100))
