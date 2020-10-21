import pickle
import Seq
import os
import matplotlib.pyplot as plt
from Telomere_DNA_Aligner import align_to_chromosomes
from Inter_Telo_Aligner import alignment_between_telos

MIN_SEQ_LEN_TO_PRINT = 150

MIN_SUBTELO_LENGTH = 50

MIN_TELO_LENGTH = 500


def analyze_telos(teloes):
    telos_lengths = []
    subtelos_length = []
    motif_dic = {}
    for telo in teloes:
        # print(telo.rec_num)
        if telo.longest_telomere_len > MIN_TELO_LENGTH:
            telos_lengths.append(telo.longest_telomere_len)
            if telo.non_telomeric_parts:
                subtelo_len = 0
                for s in telo.non_telomeric_parts:
                    subtelo_len = max(subtelo_len,len(s))
                if subtelo_len > MIN_SUBTELO_LENGTH:
                    subtelos_length.append(subtelo_len)
        if telo.len > MIN_SEQ_LEN_TO_PRINT:
            telo.print_doc()
        for motif in telo.motif_dict:
            if motif in motif_dic:
                motif_dic[motif] += telo.motif_dict[motif]
            else:
                motif_dic[motif] = telo.motif_dict[motif]
    plt.hist(telos_lengths, bins=50, color='orange')
    # plt.xscale('log')
    plt.title('Telomere length')
    plt.xlabel('# Base-Pairs')
    plt.ylabel('# Reads')
    # plt.show()
    plt.savefig('Telomere_length.jpg')
    plt.hist(subtelos_length, bins=50, color='orange')
    plt.title('Sub-Telomere length')
    plt.xlabel('# Base-Pairs')
    # plt.xscale('log')
    plt.ylabel('# Reads')
    # plt.show()
    plt.savefig('Sub_Telomere_length.jpg')
    d_view = [(v, k) for k, v in motif_dic.items()]
    d_view.sort(reverse=True)  # natively sort tuples by first element
    val_tot = 0
    for v, k in d_view:
        val_tot += v
    for v, k in d_view:
        print("%s: %d %d" % (k, v, (v / val_tot) * 100))

if __name__ == '__main__':
    path_to_teloes = '.'
    teloes = []
    for filename in os.listdir(path_to_teloes):
        _, file_extension = os.path.splitext(filename)
        if file_extension != '.obj':
            continue
        file = open(os.path.join(path_to_teloes, filename), 'rb')
        telo = pickle.load(file)
        file.close()
        telo.generate_seq()
        teloes.append(telo)
    print("Telo_Analyzer")
    print("--------------------------------------------------------------------")
    analyze_telos(teloes)
    print("--------------------------------------------------------------------")
    print('\n')
    print("Telomere_DNA_Aligner")
    print("--------------------------------------------------------------------")
    align_to_chromosomes(teloes)
    print("--------------------------------------------------------------------")
    print('\n')
    print("Inter_Telo_Aligner")
    print("--------------------------------------------------------------------")
    alignment_between_telos(teloes)
