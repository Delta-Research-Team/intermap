# Created by rglez at 12/29/24
"""
Runner for InterMap
"""
import sys
import time
from argparse import Namespace
from os.path import basename, join
from pprint import pformat

import numpy as np
import rgpack.generals as gnl
from numba import set_num_threads
from tqdm import tqdm

import intermap.config as conf
import intermap.cutoffs as cf
import intermap.interdict as idt
import intermap.topo_trajs as tt
from intermap import commons as cmn, macro
from intermap.indices import IndexManager


# todo: check docstrings


def run(mode='production'):
    """
    Main function to run InterMap
    """
    # =========================================================================
    # Parsing the configuration file
    # =========================================================================

    if mode == 'production':
        if len(sys.argv) != 2:
            raise ValueError(
                '\nInterMap syntax is: intermap path-to-config-file')
        config_path = sys.argv[1]
    elif mode == 'debug':
        config_path = 'tests/imaps/imap1.cfg'
    else:
        raise ValueError('Only modes allowed are production and running')
    # %%
    start_time = time.time()
    config = conf.ConfigManager(config_path, conf.allowed_parameters)
    args = Namespace(**config.config_args)
    log_path = join(args.output_dir, f"{basename(args.job_name)}_InterMap.log")
    logger = cmn.start_logger(log_path)
    logger.info(f"Starting InterMap with the following parameters:"
                f"\n Job name: {args.job_name}"
                f"\n Number of processors: {args.n_procs}"
                f"\n Output directory: {args.output_dir}"
                f"\n Topology: {args.topology}"
                f"\n Trajectory: {args.trajectory}"
                f"\n Start frame: {args.start}"
                f"\n Last frame: {args.last}"
                f"\n Stride: {args.stride}"
                f"\n Chunk size: {args.chunk_size}"
                f"\n Selection 1: {args.selection_1}"
                f"\n Selection 2: {args.selection_2}"
                f"\n Min prevalence: {args.min_prevalence}"
                f"\n Report's format: {args.format}"
                )
    # =========================================================================
    # Load the indices & interactions to compute
    # =========================================================================
    iman = IndexManager(args.topology, args.trajectory, args.selection_1,
                        args.selection_2, args.interactions)

    # Indices & selections for computing interactions
    sel_idx, s1_indices, s2_indices = iman.sel_idx, iman.sel1_idx, iman.sel2_idx
    vdw_radii, max_vdw = iman.radii, iman.get_max_vdw_dist()
    hydrophobes = iman.hydroph
    anions, cations = iman.anions, iman.cations
    metal_donors, metal_acceptors = iman.metal_don, iman.metal_acc
    hb_hydros, hb_donors, hb_acc = iman.hb_H, iman.hb_D, iman.hb_A
    xb_halogens, xb_donors, xb_acc = iman.xb_H, iman.xb_D, iman.xb_A
    rings = iman.rings

    # Names of the selected atoms
    universe = iman.universe
    atnames = universe.atoms.names[sel_idx]
    resnames = universe.atoms.resnames[sel_idx]
    resids = universe.atoms.resids[sel_idx]
    names = [f"{resnames[i]}_{resids[i]}_{atnames[i]}" for i, x in
             enumerate(sel_idx)]

    # Chunks of frames to analyze
    n_frames = iman.n_frames
    last = tt.parse_last_param(args.last, n_frames)
    traj_frames = np.arange(args.start, last, args.stride)
    logger.info(f"Number of frames to consider (start:last:stride): "
                f"{traj_frames.size} ({args.start}:{last}:{args.stride})")
    logger.info(f"Number of atoms to consider: {iman.sel_idx.size}")

    # =========================================================================
    # Parsing the interactions & cutoffs
    # =========================================================================
    all_inters, all_cutoffs = cf.get_inters_cutoffs(args.cutoffs)
    to_compute = iman.interactions
    selected_aro, selected_others, cutoffs_aro, cutoffs_others = \
        cmn.get_cutoffs_and_inters(to_compute, all_inters, all_cutoffs)
    len_others = len(selected_others) if not 'None' in selected_others else 0
    len_aro = len(selected_aro) if not 'None' in selected_aro else 0
    cutoffs_str = {x: args.cutoffs[x] for x in args.cutoffs if x in to_compute}
    logger.info(f"Interactions to compute:\n {pformat(to_compute)}")
    logger.debug(f"Cutoffs parsed:\n {pformat(cutoffs_str)}")

    # =========================================================================
    # Estimating memory allocation
    # =========================================================================
    logger.info(f"Estimating memory allocation")
    set_num_threads(args.n_procs)
    n_frames, n_samples = len(universe.trajectory), 100
    sub = universe.trajectory[::n_frames // n_samples]
    f4 = np.float32
    positions = np.asarray([ts.positions.copy() for ts in sub], dtype=f4)

    ijf_shape, inters_shape, mb1, mb2, v_size, h_size = macro.estimate(
        positions, args.chunk_size, s1_indices, s2_indices, cations, rings,
        cutoffs_aro, selected_aro, len_aro, anions, hydrophobes, metal_donors,
        metal_acceptors, vdw_radii, max_vdw, hb_hydros, hb_donors, hb_acc,
        xb_halogens, xb_donors, xb_acc, cutoffs_others, selected_others,
        len_others)

    logger.debug(f"Allocated space for interactions:"
                 f" ~{mb1} MB ({args.chunk_size}, {v_size}, {h_size})")
    logger.debug(f"Allocated space for coordinates:"
                 f" ~{mb2} MB ({args.chunk_size}, {sel_idx.size}, 3) ")

    # =========================================================================
    # Compiling the parallel function
    # =========================================================================
    logger.info("Compiling the parallel function")
    xyz_test = positions[:1]
    tree_test = cmn.get_trees(xyz_test, s2_indices)
    _, _ = macro.runpar(xyz_test, tree_test, ijf_shape, inters_shape,
                        len_others, len_aro, s1_indices, s2_indices, anions,
                        cations, hydrophobes, metal_donors, metal_acceptors,
                        vdw_radii, max_vdw, hb_hydros, hb_donors, hb_acc,
                        xb_halogens, xb_donors, xb_acc, rings, cutoffs_others,
                        selected_others, cutoffs_aro, selected_aro)

    # =========================================================================
    # Fill the interaction dictionary
    # =========================================================================
    logger.info(f"Starting to compute InterMap interactions")
    fmt, min_prev = args.format, args.min_prevalence
    inters = np.asarray([x for x in selected_others if x != 'None'] +
                        [x for x in selected_aro if x != 'None'])
    inter_dict = idt.InterDict(fmt, min_prev, traj_frames, names, inters)
    total_pairs, total_inters = 0, 0
    N = traj_frames.size // args.chunk_size
    chunks = tt.split_in_chunks(traj_frames, args.chunk_size)
    for i, frames in tqdm(enumerate(chunks), total=N,
                          desc='Detecting Interactions', unit='chunk'):

        xyz_chunk = tt.get_coordinates(universe, frames, sel_idx)
        trees_chunk = cmn.get_trees(xyz_chunk, s2_indices)

        ijf_chunk, inters_chunk = macro.runpar(
            xyz_chunk, trees_chunk, ijf_shape, inters_shape, len_others,
            len_aro, s1_indices, s2_indices, anions, cations, hydrophobes,
            metal_donors, metal_acceptors, vdw_radii, max_vdw, hb_hydros,
            hb_donors, hb_acc, xb_halogens, xb_donors, xb_acc, rings,
            cutoffs_others, selected_others, cutoffs_aro, selected_aro)

        total_pairs += ijf_chunk.shape[0]
        total_inters += inters_chunk.sum()

        # Filling the interaction dictionary
        if ijf_chunk.shape[0] > 0:
            col_frames = frames[ijf_chunk[:, 2]]
            ijf_chunk[:, 2] = col_frames
            inter_dict.fill(ijf_chunk, inters_chunk)

    # =========================================================================
    # Save the interactions
    # =========================================================================
    inter_dict.pack()
    job_name = basename(args.job_name)
    pickle_name = f"{job_name}_InterMap.pickle"
    pickle_path = join(args.output_dir, pickle_name)
    logger.info(f"Saving the interactions in {pickle_path}")
    gnl.pickle_to_file(inter_dict.dict, pickle_path)

    # =========================================================================
    # Timing
    # =========================================================================
    tot = round(time.time() - start_time, 2)
    ldict = len(inter_dict.dict)
    logger.info(f"Total number of unique atom pairs detected:, {ldict}")
    logger.info(f"Total number of interactions detected: {total_inters}")
    logger.info(f"Total elapsed time: {tot} s")
    logger.info(f"Normal termination of InterMap job '{job_name}'")


# run(mode='debug')
