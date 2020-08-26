__author__ = 'lachesis'
import os, sys
import pdb
import numpy as np


def make_DP_matrix(cds0, cds1, cds2, start_dict, stop_dict, frame_shift_penalty=-20):
    N = len(cds0)  # note: cds1 and cds2 could be 1 shorter
    E = np.zeros((N, 3), dtype=np.int)
    O = {}

    CDS = {0: cds0, 1: cds1, 2: cds2}
    # initialize E[1, j] = CDS[1, j]
    E[0, 0] = cds0[0] + (10 * 0 in start_dict[0]) * cds0[0]
    E[0, 1] = cds1[0] + (10 * 0 in start_dict[1]) * cds1[0]
    E[0, 2] = cds2[0] + (10 * 0 in start_dict[2]) * cds2[0]

    for i in range(1, N):
        for j in range(3):
            if i >= len(CDS[j]): continue
            #pdb.set_trace()
            if i-1 not in stop_dict[j]:
                E[i, j] = E[i-1, j]
                O[(i, j)] = (i-1, j)
            else:
                E[i, j] = 0
                O[(i, j)] = '*'
            for k in range(3):
                if k!=j and i-1 not in stop_dict[k]:
                    tmp = E[i-1, k] + frame_shift_penalty
                    if tmp > E[i, j]:
                        O[(i, j)] = (i-1, k)
                        E[i, j] = tmp
            for k in range(j):
                if i not in stop_dict[k]:
                    tmp = E[i, k] + frame_shift_penalty
                    if tmp > E[i, j]:
                        O[(i, j)] = (i, k)
                        E[i, j] = tmp
            if CDS[j][i] == 1:
                E[i, j] += 1 #+ (10 * i in start_dict[j])
            else:
                E[i, j] -= 1

    return E, O
