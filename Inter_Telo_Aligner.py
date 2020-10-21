# get telo objects, align them, and check their location compared to the reference

import mappy as mp
import Seq

SEQ_EDGE = 20

MAX_ALIGN_ERR = 0.1

MIN_ALIGN_LEN = 1000


def alignment_between_telos(telos):

    for index, telo in enumerate(telos):
        print(telo.rec_num)
        aligner = mp.Aligner(seq=telo.sequence, preset='map-ont')
        if not aligner:
            raise Exception("ERROR: failed to load/build index")
        for sub_index, align in enumerate(telos):
            if index == sub_index:
                continue
            is_aligned = False
            for hit in aligner.map(align.sequence):
                if ((hit.blen > MIN_ALIGN_LEN) and ((hit.NM / hit.blen) <= MAX_ALIGN_ERR)) and hit.is_primary:
                    ref_start = hit.r_st < SEQ_EDGE
                    ref_end = hit.r_en > (telo.len - SEQ_EDGE)
                    quer_start = hit.q_st < SEQ_EDGE
                    quer_end = hit.q_en > (align.len - SEQ_EDGE)
                    if (ref_start and ref_end) or (quer_start and quer_end) or (ref_start and quer_end) or (
                            ref_end and quer_start):
                        is_aligned = True
                        print(
                            align.rec_num + " refrence : " + str(hit.r_st) + " - " + str(hit.r_en) + " query : " + str(
                                hit.q_st) + " - " + str(hit.q_en) + " blen: " + str(hit.blen) + " NM: " + str(hit.NM))
            if is_aligned:
                print('\n')
        print('\n')
