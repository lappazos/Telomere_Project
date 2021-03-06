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


def chromosome_matcher(hit, chromo_dict,chromosome_graph,telo_length,chromo_count):
    handle = Entrez.efetch(db="nucleotide", id=hit.ctg, rettype="gb", retmode="text")
    lines = handle.readlines()
    try:
        chromosome = lines[1].split("chromosome ")[1].split(", ")[0].split(" ")[0]
    except IndexError:
        chromosome = "Unknown"
    print(hit.ctg + " : " + str(hit.r_st) + " - " + str(hit.r_en) + " : " + str(chromosome))
    print('\n')
    chr_len = int(lines[0].split(" bp")[0].split(" ")[-1])
    if hit.r_st <0 or chr_len < 0 or hit.r_st > chr_len:
        return
    strand = 1 if hit.strand ==1 else 0
    chromosome_graph[chromosome][strand].append(hit.r_st/chr_len)
    print('chromosome '+str(chromosome)+' strand '+str(strand) +' relative location  '+str(hit.r_st/chr_len))
    if ('scaffold' in lines[1]) or ('patch' in lines[1]):
        chromo_dict[chromosome][strand][UNKNOWN_INDEX] += 1
        chromo_count[1][0] +=1
        chromo_count[1][1] +=telo_length
        handle.close()
        return
    elif hit.r_st < chr_len * CHR_BEGIN:
        chromo_dict[chromosome][strand][BEGIN_INDEX] += 1
        chromo_count[0][0] += 1
        chromo_count[0][1] += telo_length
        handle.close()
        return
    elif hit.r_en > chr_len * CHR_END:
        chromo_dict[chromosome][strand][END_INDEX] += 1
        chromo_count[0][0] += 1
        chromo_count[0][1] += telo_length
        handle.close()
        return
    chromo_dict[chromosome][strand][UNKNOWN_INDEX] += 1
    chromo_count[1][0] += 1
    chromo_count[1][1] += telo_length
    handle.close()


def align_to_chromosomes(teloes):
    aligner = mp.Aligner('../../GRCh38_latest_genomic.fna.gz', preset='map-ont')
    if not aligner: raise Exception("ERROR: failed to load/build index")

    chromosome_dict = {}
    chromosome_graph = {}
    for i in range(1, NUM_CHROMOSOMES + 1):
        chromosome_dict[str(i)] = [[0, 0, 0],[0, 0, 0]]
        chromosome_graph[str(i)] = [[],[]]
    chromosome_dict["X"] = [[0, 0, 0],[0, 0, 0]]
    chromosome_graph["X"] = [[],[]]
    chromosome_dict["Y"] = [[0, 0, 0],[0, 0, 0]]
    chromosome_graph["Y"] = [[],[]]
    chromosome_dict["Unknown"] = [[0, 0, 0],[0, 0, 0]]
    chromosome_graph["Unknown"] = [[],[]]
    chromo_count = [[0,0],[0,0]]

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
            chromosome_matcher(best, chromosome_dict,chromosome_graph,telo.longest_telomere_len,chromo_count)

    for j in chromosome_dict:
        if (chromosome_dict[j][0] != 0) or (chromosome_dict[j][1] != 0) or (chromosome_dict[j][2] != 0):
            print(j + " " + str(chromosome_dict[j]))

    for chromosome in chromosome_graph:
        plt.axhline(y=1, color='b', linestyle='-')
        plt.axhline(y=0, color='b', linestyle='-')
        plt.axhline(y=8, color='b', linestyle='-')
        plt.axhline(y=-7, color='b', linestyle='-')
        for dot in chromosome_graph[chromosome][1]:
            plt.plot(dot, 0, 'ro')
        for dot in chromosome_graph[chromosome][0]:
            plt.plot(dot, 1, 'ro')
        plt.savefig('Chromosome_'+chromosome+'.jpg')
        plt.close()

    print("Average telo legth on edges - "+str(chromo_count[0][1]/chromo_count[0][0])+"\n")
    print("Average telo legth in center - "+str(chromo_count[1][1]/chromo_count[1][0])+"\n")

