# Created by rglez at 1/28/25
from argparse import Namespace

import numpy as np
from numba import njit
from numba_kdtree import KDTree as nckd

import intermap.commons as cmn
import managers.config as conf
import managers.cutoffs as cf
from intermap.interactions.aro import aro


def test_kdtree(iman, conf_path):
    # generate a set of 10 xyz coordinates
    @njit(parallel=False, cache=True)
    def insider(xyz, indices, dist_cut):
        s2_tree = nckd(xyz[indices])
        ball_1 = s2_tree.query_radius(xyz[indices], dist_cut)
        return ball_1

    dist_cut = 4.5
    xyz = np.random.rand(10, 3)
    indices = np.random.randint(0, 10, 5)
    ball_1 = insider(xyz, indices, dist_cut)
    assert sum([len(x) for x in ball_1]) != 0


def test_aro(iman, conf_path):
    xyz = iman.universe.atoms.positions
    k = 0
    s1_indices, s2_indices = iman.s1_idx, iman.s2_idx
    cations = iman.cations
    rings = iman.rings

    config = conf.ConfigManager(conf_path, conf.allowed_parameters)
    args = Namespace(**config.config_args)
    all_inters, all_cutoffs = cf.get_inters_cutoffs()
    to_compute = iman.interactions
    selected_aro, selected_others, cutoffs_aro, cutoffs_others = \
        cmn.get_cutoffs_and_inters(to_compute, all_inters, all_cutoffs)

    iff_aro, inters_aro = aro(xyz, k, s1_indices, s2_indices, cations, rings,
                              cutoffs_aro, selected_aro,,, , ,
    assert iff_aro.any()
    assert inters_aro.any()
