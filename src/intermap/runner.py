# Created by rglez at 12/29/24
"""
Runner for InterMap
"""
import logging
import time
from argparse import Namespace
from os.path import basename, join

import numpy as np
from numba import set_num_threads
from tqdm import tqdm

from intermap import commons as cmn
from intermap.interactions import aro
from intermap.interactions.runners import estimate, runpar
from intermap.interactions.waters import wb1
from intermap.managers.config import ConfigManager
from intermap.managers.container import ContainerManager
from intermap.managers.cutoffs import CutoffsManager
from intermap.managers.indices import IndexManager


# High Priority
# todo: update filling dict when water
# todo: implement granularity as a way to condense information

# Medium Priority
# todo: put n_samples / n_factor in config

# Low Priority
# todo: assert identity against  prolif, again
# todo: Reorganize the code
# todo: join into the same function: aro.get_trees, cmn.get_trees
# todo: check logging
# todo: check docstrings
# todo: start writing tests

# done: check hard-coded cutoffs
# done: implement selecting interactions from config
# done: assert changing cfg does not interfere with cache=True
# done: investigate recompilation issues
# done: check the interactions naming / parsing
# done: do not gather balls and trees outside runpar / estimate functions
# done: rename  to CloseContact in all files

# %%

def run():
    """
    Main function to run InterMap
    """
    # %%=======================================================================
    # 1. Parse the configuration file
    # =========================================================================
    start_time = time.time()
    logger = logging.getLogger('InterMapLogger')
    config = ConfigManager(mode='debug')
    # config = ConfigManager()
    args = Namespace(**config.config_args)
    set_num_threads(args.n_procs)

    # =========================================================================
    # 2. Load the indices & interactions to compute
    # =========================================================================
    iman = IndexManager(args)
    (sel_idx, s1_idx, s2_idx, s1_cat, s2_cat, s1_cat_idx, s2_cat_idx, s1_rings,
     s2_rings, s1_rings_idx, s2_rings_idx, s1_aro_idx, s2_aro_idx, xyz_aro_idx,
     vdw_radii, max_vdw, hydroph, met_don, met_acc, hb_hydr, hb_don, hb_acc,
     xb_hal, xb_don, xb_acc, waters, anions, cations, rings, overlap, universe,
     resid_names, atom_names, resindex, n_frames, traj_frames,
     inters_requested) = (

        iman.sel_idx, iman.s1_idx, iman.s2_idx, iman.s1_cat, iman.s2_cat,
        iman.s1_cat_idx, iman.s2_cat_idx, iman.s1_rings, iman.s2_rings,
        iman.s1_rings_idx, iman.s2_rings_idx, iman.s1_aro_idx, iman.s2_aro_idx,
        iman.xyz_aro_idx, iman.vdw_radii, iman.get_max_vdw_dist(),
        iman.hydroph, iman.met_don, iman.met_acc, iman.hb_hydro, iman.hb_don,
        iman.hb_acc, iman.xb_hal, iman.xb_don, iman.xb_acc, iman.waters,
        iman.anions, iman.cations, iman.rings, iman.overlap, iman.universe,
        iman.resid_names, iman.atom_names, iman.resindex, iman.n_frames,
        iman.traj_frames, iman.inters_requested)

    # =========================================================================
    # 3. Parse the interactions & cutoffs
    # =========================================================================
    cuts = CutoffsManager(args, iman)
    (cuts_aro, cuts_others, selected_aro, selected_others, len_aro, len_others,
     max_dist_aro, max_dist_others) = (

        cuts.cuts_aro, cuts.cuts_others, cuts.selected_aro,
        cuts.selected_others, cuts.len_aro, cuts.len_others, cuts.max_dist_aro,
        cuts.max_dist_others)

    # =========================================================================
    # 4. Estimating memory allocation
    # =========================================================================
    ijf_shape, inters_shape = estimate(universe, xyz_aro_idx, args.chunk_size,
                                       s1_idx, s2_idx, cations, s1_cat_idx,
                                       s2_cat_idx, s1_cat, s2_cat, s1_rings,
                                       s2_rings, s1_rings_idx, s2_rings_idx,
                                       s1_aro_idx, s2_aro_idx, cuts_aro,
                                       selected_aro, len_aro, anions, hydroph,
                                       met_don, met_acc, vdw_radii, hb_hydr,
                                       hb_don, hb_acc, xb_hal, xb_don, xb_acc,
                                       cuts_others, selected_others,
                                       len_others, max_dist_aro,
                                       max_dist_others, overlap)

    # =========================================================================
    # 5. Trim the trajectory
    # =========================================================================
    chunk_frames = list(cmn.split_in_chunks(traj_frames, args.chunk_size))
    n_chunks = traj_frames.size // args.chunk_size
    trajiter = cmn.trajiter(universe, chunk_frames, sel_idx)
    contiguous = list(cmn.split_in_chunks(np.arange(traj_frames.size),
                                          args.chunk_size))

    # %%=======================================================================
    # 6. Detect the interactions
    # =========================================================================
    container = ContainerManager(args, iman, cuts)
    hb_idx = container.hb_idx
    detect_wb = (waters.size > 0) and (hb_idx.size > 0)
    total_pairs, total_inters = 0, 0
    for i, xyz_chunk in tqdm(enumerate(trajiter),
                             desc='Detecting Interactions',
                             unit='chunk', total=n_chunks, ):
        s1_centrs, s2_centrs, xyzs_aro = aro.get_aro_xyzs(
            xyz_chunk, s1_rings, s2_rings, s1_cat, s2_cat)

        trees_aro = aro.get_trees(xyzs_aro, s2_aro_idx)
        trees_others = cmn.get_trees(xyz_chunk, s2_idx)

        ijf_chunk, inters_chunk = runpar(
            xyz_chunk, xyzs_aro, xyz_aro_idx, trees_others, trees_aro,
            ijf_shape, inters_shape, s1_idx, s2_idx, anions, cations,
            s1_cat_idx, s2_cat_idx, hydroph, met_don, met_acc, vdw_radii,
            max_vdw, hb_hydr, hb_don, hb_acc, xb_hal, xb_don, xb_acc, s1_rings,
            s2_rings, s1_rings_idx, s2_rings_idx, s1_aro_idx, s2_aro_idx,
            cuts_others, selected_others, cuts_aro, selected_aro, overlap)

        total_pairs += ijf_chunk.shape[0]
        total_inters += inters_chunk.sum()

        # Fill interactions in ijf and wb interactions in ijkf
        frames = contiguous[i]
        if ijf_chunk.shape[0] > 0:
            ijf_chunk[:, 2] = frames[ijf_chunk[:, 2]]
            container.fill(ijf_chunk, inters_chunk)
            if detect_wb:
                ijwf = wb1(ijf_chunk, inters_chunk, waters, hb_idx)
                container.fill(ijwf, inters='wb')

    # %%=======================================================================
    # 7. Save the interactions
    # =========================================================================
    out_name = f"{basename(args.job_name)}_InterMap.tsv"
    pickle_path = join(args.output_dir, out_name)
    container.save(pickle_path)

    # %%=======================================================================
    # 8. Normal termination
    # =========================================================================
    tot = round(time.time() - start_time, 2)
    ldict = len(container.dict)
    print('\n\n')
    logger.info(
        f"Normal termination of InterMap job '{basename(args.job_name)}'\n\n"
        f" Interactions saved in {pickle_path}\n"
        f" Total number of unique atom pairs detected: {ldict}\n"
        f" Total number of interactions detected: {total_inters}\n"
        f" Elapsed time: {tot} s")

# run()
