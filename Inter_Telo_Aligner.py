# get telo objects, align them, and check their location compared to the reference

import mappy as mp
import Seq

SEQ_EDGE = 20

MAX_ALIGN_ERR = 0.1

MIN_ALIGN_LEN = 500


def alignment_between_telos(telos):

    for index, telo in enumerate(telos):
        first_print = True
        aligner = mp.Aligner(seq=telo.sequence, preset='map-ont')
        if not aligner:
            raise Exception("ERROR: failed to load/build index")
        for sub_index, align in enumerate(telos):
            if index == sub_index:
                continue
            for hit in aligner.map(align.sequence):
                if ((hit.blen > MIN_ALIGN_LEN) and ((hit.NM / hit.blen) <= MAX_ALIGN_ERR)) and hit.is_primary:
                    ref_start = hit.r_st < SEQ_EDGE
                    ref_end = hit.r_en > (telo.len - SEQ_EDGE)
                    quer_start = hit.q_st < SEQ_EDGE
                    quer_end = hit.q_en > (align.len - SEQ_EDGE)
                    if (ref_start and ref_end) or (quer_start and quer_end) or (ref_start and quer_end) or (
                            ref_end and quer_start):
                        if first_print:
                            print(telo.rec_num)
                            first_print = False
                        print(
                            align.rec_num + " refrence : " + str(hit.r_st) + " - " + str(hit.r_en) + " query : " + str(
                                hit.q_st) + " - " + str(hit.q_en) + " blen: " + str(hit.blen) + " NM: " + str(hit.NM))
        if not first_print:
            print('\n')
