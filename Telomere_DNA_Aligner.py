# need to be in the same folder with telo objects, Seq.py file, and Genome index
# generates dictionary with alignments
# get the alignment dictionary, and match the different alignments to chromosome location.

import mappy as mp
import Seq
import matplotlib.pyplot as plt

from Bio import Entrez

NUM_CHROMOSOMES = 22

UNKNOWN_INDEX = 0

END_INDEX = 2

BEGIN_INDEX = 1

CHR_END = 0.97

CHR_BEGIN = 0.03

Entrez.email = "Your.Name.Here@example.org"

MIN_ALIGNMENT_LENGTH = 50


def chromosome_matcher(hit, chromo_dict,chromosome_graph):
    handle = Entrez.efetch(db="nucleotide", id=hit.ctg, rettype="gb", retmode="text")
    lines = handle.readlines()
    try:
        chromosome = lines[1].split("chromosome ")[1].split(", ")[0].split(" ")[0]
    except IndexError:
        chromosome = "Unknown"
    print(hit.ctg + " : " + str(hit.r_st) + " - " + str(hit.r_en) + " : " + str(chromosome))
    print('\n')
    chr_len = int(lines[0].split(" bp")[0].split(" ")[-1])
    if hit.r_st >0 and chr_len > 0:
        chromosome_graph[chromosome].append(hit.r_st/chr_len)
    if ('scaffold' in lines[1]) or ('patch' in lines[1]):
        chromo_dict[chromosome][UNKNOWN_INDEX] += 1
        return
    elif hit.r_st < chr_len * CHR_BEGIN:
        chromo_dict[chromosome][BEGIN_INDEX] += 1
        return
    elif hit.r_en > chr_len * CHR_END:
        chromo_dict[chromosome][END_INDEX] += 1
        return
    chromo_dict[chromosome][UNKNOWN_INDEX] += 1
    handle.close()


def align_to_chromosomes(teloes):
    aligner = mp.Aligner('../../GRCh38_latest_genomic.fna.gz', preset='map-ont')
    if not aligner: raise Exception("ERROR: failed to load/build index")

    chromosome_dict = {}
    chromosome_graph = {}
    for i in range(1, NUM_CHROMOSOMES + 1):
        chromosome_dict[str(i)] = [0, 0, 0]
        chromosome_graph[str(i)] = []
    chromosome_dict["X"] = [0, 0, 0]
    chromosome_graph["X"] = []
    chromosome_dict["Y"] = [0, 0, 0]
    chromosome_graph["Y"] = []
    chromosome_dict["Unknown"] = [0, 0, 0]
    chromosome_graph["Unknown"] = []

    for telo in teloes:
        first_print=True
        aligns_dict = {}
        for to_align in telo.non_telomeric_parts:
            if len(to_align) < MIN_ALIGNMENT_LENGTH:
                continue
            for hit in aligner.map(to_align):
                if hit.is_primary:
                    aligns_dict[hit.mlen] = hit
                    if(first_print):
                        print(telo.rec_num)
                        first_print = False
                    print(hit.ctg + " : " + str(hit.r_st) + " - " + str(hit.r_en) + " starnd: " + str(
                        hit.strand) + " blen: " + str(hit.blen) + " mlen: " + str(hit.mlen) + " NM: " + str(hit.NM))
        if aligns_dict:
            best = aligns_dict[max(aligns_dict)]
            chromosome_matcher(best, chromosome_dict,chromosome_graph)

    for j in chromosome_dict:
        if (chromosome_dict[j][0] != 0) or (chromosome_dict[j][1] != 0) or (chromosome_dict[j][2] != 0):
            print(j + " " + str(chromosome_dict[j]))

    for chromosome in chromosome_graph:
        plt.axhline(y=0.5, color='b', linestyle='-')
        for dot in chromosome_graph[chromosome]:
            plt.plot(dot, 0.5, 'ro')
        plt.savefig('Chromosome_'+chromosome+'.jpg')
        plt.close()

