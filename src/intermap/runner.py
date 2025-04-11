# Created by rglez at 12/29/24
"""
Runner for InterMap
"""
import logging
import time
from argparse import Namespace
from os.path import basename, join
from pathlib import Path

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
# todo: Reorganize the code
# todo: start writing tests
# todo: assert identity against  prolif, again
# todo: check docstrings


def workflow(args):
    """
    Entry point to run the InterMap workflow.
    """
    # %%=======================================================================
    # 1. Start (logging, print logo, set number of threads)
    # =========================================================================
    start_time = time.time()
    logger = logging.getLogger('InterMapLogger')
    set_num_threads(args.n_procs)

    if isinstance(args, dict):
        args = Namespace(**args)
        cmn.print_colored_ascii()

    # %%=======================================================================
    # 2. Load the indices & interactions to compute
    # =========================================================================
    iman = IndexManager(args)
    (sel_idx, s1_idx, s2_idx, s1_cat, s2_cat, s1_cat_idx, s2_cat_idx, s1_rings,
     s2_rings, s1_rings_idx, s2_rings_idx, s1_aro_idx, s2_aro_idx, xyz_aro_idx,
     vdw_radii, max_vdw, hydroph, met_don, met_acc, hb_hydr, hb_don, hb_acc,
     xb_hal, xb_don, xb_acc, waters, anions, cations, rings, overlap, universe,
     resid_names, atom_names, resconv, n_frames, traj_frames,
     inters_requested) = (

        iman.sel_idx, iman.s1_idx, iman.s2_idx, iman.s1_cat, iman.s2_cat,
        iman.s1_cat_idx, iman.s2_cat_idx, iman.s1_rings, iman.s2_rings,
        iman.s1_rings_idx, iman.s2_rings_idx, iman.s1_aro_idx, iman.s2_aro_idx,
        iman.xyz_aro_idx, iman.vdw_radii, iman.get_max_vdw_dist(),
        iman.hydroph, iman.met_don, iman.met_acc, iman.hb_hydro, iman.hb_don,
        iman.hb_acc, iman.xb_hal, iman.xb_don, iman.xb_acc, iman.waters,
        iman.anions, iman.cations, iman.rings, iman.overlap, iman.universe,
        iman.resid_names, iman.atom_names, iman.resconv, iman.n_frames,
        iman.traj_frames, iman.inters_requested)

    # %%=======================================================================
    # 3. Parse the interactions & cutoffs
    # =========================================================================
    cuts = CutoffsManager(args, iman)
    (cuts_aro, cuts_others, selected_aro, selected_others, len_aro, len_others,
     max_dist_aro, max_dist_others) = (

        cuts.cuts_aro, cuts.cuts_others, cuts.selected_aro,
        cuts.selected_others, cuts.len_aro, cuts.len_others, cuts.max_dist_aro,
        cuts.max_dist_others)

    # %%=======================================================================
    # 4. Estimating memory allocation
    # =========================================================================
    atomic = True if args.resolution == 'atom' else False
    ijf_shape, inters_shape = estimate(
        universe, xyz_aro_idx, args.chunk_size, s1_idx, s2_idx, cations,
        s1_cat_idx, s2_cat_idx, s1_cat, s2_cat, s1_rings, s2_rings,
        s1_rings_idx, s2_rings_idx, s1_aro_idx, s2_aro_idx, cuts_aro,
        selected_aro, len_aro, anions, hydroph, met_don, met_acc, vdw_radii,
        hb_hydr, hb_don, hb_acc, xb_hal, xb_don, xb_acc, cuts_others,
        selected_others, len_others, max_dist_aro, max_dist_others, overlap,
        atomic, resconv)

    # %%=======================================================================
    # 5. Trim the trajectory
    # =========================================================================
    chunk_frames = list(cmn.split_in_chunks(traj_frames, args.chunk_size))
    n_chunks = traj_frames.size // args.chunk_size
    trajiter = cmn.trajiter(universe, chunk_frames, sel_idx)
    contiguous = list(cmn.split_in_chunks(np.arange(traj_frames.size),
                                          args.chunk_size))

    # %%=======================================================================
    # 6. Detect the interactions & Fill the container
    # =========================================================================
    total_pairs, total_inters = 0, 0
    container = ContainerManager(args, iman, cuts)
    for i, xyz_chunk in tqdm(enumerate(trajiter),
                             desc='Detecting Interactions',
                             unit='chunk', total=n_chunks, ):
        # 6.1 Get centroids & coordinates of aromatic rings
        s1_ctrs, s2_ctrs, xyzs_aro = aro.get_aro_xyzs(
            xyz_chunk, s1_rings, s2_rings, s1_cat, s2_cat)

        # 6.2 Get the kdtrees of aromatic & non-aromatic coordinates
        trees_aro = cmn.get_trees(xyzs_aro, s2_aro_idx)
        trees_others = cmn.get_trees(xyz_chunk, s2_idx)

        # 6.3 Detect interactions in parallel
        ijf_chunk, inters_chunk = runpar(
            xyz_chunk, xyzs_aro, xyz_aro_idx, trees_others, trees_aro,
            ijf_shape, inters_shape, s1_idx, s2_idx, anions, cations,
            s1_cat_idx, s2_cat_idx, hydroph, met_don, met_acc, vdw_radii,
            max_vdw, hb_hydr, hb_don, hb_acc, xb_hal, xb_don, xb_acc, s1_rings,
            s2_rings, s1_rings_idx, s2_rings_idx, s1_aro_idx, s2_aro_idx,
            cuts_others, selected_others, cuts_aro, selected_aro, overlap,
            atomic, resconv)

        # 6.4 Update counters
        total_pairs += ijf_chunk.shape[0]
        total_inters += inters_chunk.sum()

        # 6.5 Renumber from atom to residue indices if resolution is 'residue'
        if not atomic:
            ijf_chunk[:, :2] = resconv[ijf_chunk[:, :2]]

        # 6.6 Fill the container with the interactions
        frames = contiguous[i]
        if ijf_chunk.shape[0] > 0:
            ijf_chunk[:, 2] = frames[ijf_chunk[:, 2]]
            container.fill(ijf_chunk, inters_chunk)

            # 6.7 Fill the container with the water bridges
            if container.detect_wb:
                ijwf = wb1(ijf_chunk, inters_chunk, waters, container.hb_idx,
                           resconv, atomic=atomic)
                container.fill(ijwf, inters='wb')

    # %%=======================================================================
    # 7. Save the interactions
    # =========================================================================
    out_name = f"{basename(args.job_name)}_InterMap"
    csv_path = join(args.output_dir, f'{out_name}.csv')
    container.save(csv_path)
    packed = container.pack()

    # %%=======================================================================
    # 8. Normal termination
    # =========================================================================
    tot = round(time.time() - start_time, 2)
    ldict = len(container.dict)
    pair_type = 'atom' if atomic else 'residue'
    print('\n\n')
    logger.info(
        f"Normal termination of InterMap job '{basename(args.job_name)}'\n\n"
        f" Interactions saved in {csv_path}\n"
        f" Total number of unique {pair_type} pairs detected: {ldict}\n"
        f" Total number of interactions detected: {total_inters}\n"
        f" Elapsed time: {tot} s")

    return packed


def run(mode='debug'):
    """
    Run the InterMap workflow.

    Args:
        mode (str): Mode of operation. Can be 'debug' or 'production'.
    """
    # >>>> Detect args
    config = ConfigManager(mode=mode)
    arguments = Namespace(**config.config_args)

    # >>>> Run the InterMap workflow
    return workflow(arguments)


# =============================================================================
# Running from .cfg file
# =============================================================================

if __name__ == '__main__':
    run()
