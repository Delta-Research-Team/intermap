# Created by gonzalezroy at 3/11/25
import itertools as it

import numpy as np
import numpy_indexed as npi

from intermap.interactions import geometry as aot


def wb1(ijf_chunk, inters_chunk, waters, idxs, resconv, atomic=True):

    if not atomic:
        waters = resconv[waters]

    # Detect waters in i and j rows
    i_is_wat = aot.isin(ijf_chunk[:, 0], waters)
    j_is_wat = aot.isin(ijf_chunk[:, 1], waters)

    # Select the rows that have at least one water, that is connected by a hb

    a1_xor_a2 = i_is_wat ^ j_is_wat
    hb_inters = inters_chunk[:, idxs[0]] | inters_chunk[:, idxs[1]]
    selected = a1_xor_a2 & hb_inters
    i_is_wat = i_is_wat[selected]
    j_is_wat = j_is_wat[selected]
    ijf_wat = ijf_chunk[selected]

    # Swap the rows to have the water in the first column
    t1_i = ijf_wat[i_is_wat][:, 0]
    t1_j = ijf_wat[i_is_wat][:, 1]
    t2_i = ijf_wat[j_is_wat][:, 1]
    t2_j = ijf_wat[j_is_wat][:, 0]
    i_len = t1_i.size

    ijf_wat[:i_len, 0] = t1_i
    ijf_wat[i_len:, 0] = t2_i
    ijf_wat[:i_len, 1] = t1_j
    ijf_wat[i_len:, 1] = t2_j

    # Group the interactions by the water atoms
    if ijf_wat.size == 0:
        return np.zeros((0, 4), dtype=np.int32)

    groups = npi.group_by(ijf_wat[:, [0, 2]])
    sorter = groups.index.sorter
    slices = groups.index.slices
    indices = np.split(sorter, slices[1:-1])

    # Get the interactions that have more than one water
    doubles = [x for x in indices if len(x) > 1]
    sizes = [int(len(x) * (len(x) - 1) / 2) for x in doubles]
    ijkf = np.empty((sum(sizes), 4), dtype=np.int32)
    c = 0
    for i, double in enumerate(doubles):
        section = ijf_wat[double]
        k = section[0][0]
        f = section[0][2]
        combinations = list(it.combinations(section[:, 1], 2))
        for (v, j) in combinations:
            ijkf[c] = (v, j, k, f)
            c += 1
    return ijkf
