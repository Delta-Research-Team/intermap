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
    # Starting logger
    # =========================================================================
    start_time = time.time()

    if mode == 'production':
        if len(sys.argv) != 2:
            raise ValueError(
                '\nInterMap syntax is: intermap path-to-config-file')
        config_path = sys.argv[1]
    elif mode == 'debug':
        config_path = '/home/rglez/RoyHub/intermap/tests/imaps/imap2.cfg'
    else:
        raise ValueError('Only modes allowed are production and running')
    # %%
    f4 = np.float32
    config = conf.InterMapConfig(config_path, conf.allowed_parameters)
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
    # Load the necessary indices for detecting each interaction
    # =========================================================================
    iman = IndexManager(args.topology, args.trajectory, args.selection_1,
                        args.selection_2, args.interactions)

    vdw_radii = iman.radii
    max_vdw = iman.get_max_vdw_dist()
    hydrophobes = iman.hydroph
    anions, cations = iman.anions, iman.cations
    metal_donors, metal_acceptors = iman.metal_don, iman.metal_acc
    hb_hydros, hb_donors, hb_acc = iman.hb_H, iman.hb_D, iman.hb_A
    xb_halogens, xb_donors, xb_acc = iman.xb_H, iman.xb_D, iman.xb_A
    rings = iman.rings

    # =========================================================================
    # Parsing the interactions & cutoffs
    # =========================================================================
    all_inters, all_cutoffs = cf.get_inters_cutoffs(args.cutoffs)
    to_compute = iman.interactions
    selected_aro, selected_others, cutoffs_aro, cutoffs_others = \
        cmn.get_cutoffs_and_inters(to_compute, all_inters, all_cutoffs)

    len_others, len_aro = len(selected_others), len(selected_aro)
    cutoffs_str = {x: args.cutoffs[x] for x in args.cutoffs if x in to_compute}
    logger.info(f"Interactions to compute:\n {pformat(to_compute)}")
    logger.debug(f"Cutoffs parsed:\n {pformat(cutoffs_str)}")

    # =========================================================================
    # Load & trim the trajectory
    # =========================================================================
    n_frames = iman.n_frames
    last = tt.parse_last_param(args.last, n_frames)
    traj_frames = np.arange(args.start, last, args.stride)
    logger.info(f"Number of frames to consider (start:last:stride): "
                f"{traj_frames.size} ({args.start}:{last}:{args.stride})")
    logger.info(f"Number of atoms to consider: {iman.sel_idx.size}")

    # =========================================================================
    # Naming all atoms and interactions
    # =========================================================================
    universe = iman.universe
    sel_idx = iman.sel_idx
    s1_indices = iman.sel1_idx
    s2_indices = iman.sel2_idx
    natoms = sel_idx.size
    atnames = universe.atoms.names[sel_idx]
    resnames = universe.atoms.resnames[sel_idx]
    resids = universe.atoms.resids[sel_idx]
    names = [f"{resnames[i]}_{resids[i]}_{atnames[i]}" for i, x in
             enumerate(sel_idx)]
    inters = np.asarray(selected_others.tolist() + selected_aro.tolist())

    # =========================================================================
    # Estimating memory allocation
    # =========================================================================
    logger.info(f"Estimating memory allocation")
    set_num_threads(args.n_procs)
    n_frames, n_samples = len(universe.trajectory), 10
    sub = universe.trajectory[::n_frames // n_samples]
    positions = np.asarray([ts.positions.copy() for ts in sub], dtype=f4)

    ijf_shape, inters_shape, mb1, mb2, v_size, h_size = macro.estimate(
        positions, args.chunk_size, s1_indices, s2_indices, cations, rings,
        cutoffs_aro, selected_aro, anions, hydrophobes, metal_donors,
        metal_acceptors, vdw_radii, max_vdw, hb_hydros, hb_donors, hb_acc,
        xb_halogens, xb_donors, xb_acc, cutoffs_others, selected_others)

    logger.debug(f"Allocated space for interactions:"
                 f" ~{mb1} MB ({args.chunk_size}, {v_size}, {h_size})")
    logger.debug(f"Allocated space for coordinates:"
                 f" ~{mb2} MB ({args.chunk_size}, {natoms}, 3) ")

    # =========================================================================
    # Compiling the parallel function
    # =========================================================================
    logger.info("Compiling the parallel function")
    xyz_test = positions[:1]
    tree_test = cmn.get_trees(xyz_test, s2_indices)
    _, _ = macro.run_parallel2(
        xyz_test, tree_test, ijf_shape, inters_shape, len_others, len_aro,
        s1_indices, s2_indices, anions, cations, hydrophobes, metal_donors,
        metal_acceptors, vdw_radii, max_vdw, hb_hydros, hb_donors, hb_acc,
        xb_halogens, xb_donors, xb_acc, rings, cutoffs_others,
        selected_others, cutoffs_aro, selected_aro)

    # =========================================================================
    # Fill the interaction dictionary
    # =========================================================================
    logger.info(f"Starting to compute InterMap interactions")
    fmt, min_prev = args.format, args.min_prevalence
    inter_dict = idt.InterDict(fmt, min_prev, traj_frames, names, inters)

    chunks = tt.split_in_chunks(traj_frames, args.chunk_size)
    total_pairs, total_inters = 0, 0
    for i, frames_chunk in enumerate(chunks):
        logger.info(f"Retrieving coordinates for chunk {i}")
        xyz_chunk = tt.get_coordinates(universe, frames_chunk, sel_idx)
        trees_chunk = cmn.get_trees(xyz_chunk, s2_indices)

        # Compute the interactions
        logger.info(f"Computing chunk {i}")
        ijf_chunk, inters_chunk = macro.run_parallel2(
            xyz_chunk, trees_chunk, ijf_shape, inters_shape, len_others,
            len_aro,
            s1_indices, s2_indices, anions, cations, hydrophobes, metal_donors,
            metal_acceptors, vdw_radii, max_vdw, hb_hydros, hb_donors, hb_acc,
            xb_halogens, xb_donors, xb_acc, rings, cutoffs_others,
            selected_others, cutoffs_aro, selected_aro)

        total_pairs += ijf_chunk.shape[0]
        total_inters += inters_chunk.sum()

        # Filling the interaction dictionary
        logger.info(f"Filling the interaction dictionary with chunk {i}")
        if ijf_chunk.shape[0] > 0:
            col_frames = frames_chunk[ijf_chunk[:, 2]]
            ijf_chunk[:, 2] = col_frames
            inter_dict.fill(ijf_chunk, inters_chunk)

    # =========================================================================
    # Saving
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
    logger.info(f"Normal termination of InterMap job '{job_name}'")
    logger.info(f"Total time: {tot} s")

run(mode='debug')
