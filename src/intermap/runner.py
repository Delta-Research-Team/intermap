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
import intermap.starters as start
import intermap.topo_trajs as tt
from intermap import commons as cmn
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
        config_path = 'tests/imaps/imap1.cfg'
    else:
        raise ValueError('Only modes allowed are production and running')
    # %%
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
    # Fill the interaction dictionary
    # =========================================================================
    stamp = time.time()
    logger.info(f"Starting to compute InterMap interactions")
    set_num_threads(args.n_procs)
    inter_dict = idt.InterDict(
        args.format, args.min_prevalence, traj_frames, names, inters)

    chunks = tt.split_in_chunks(traj_frames, args.chunk_size)
    to_allocate = 1
    total = 0
    total_inters = 0

    for i, frames_chunk in enumerate(chunks):
        xyz_chunk = tt.get_coordinates(universe, frames_chunk, sel_idx,
                                       natoms)
        if i == 0:
            # Estimating memory allocation
            logger.info(f"Estimating memory allocation")
            ijf_template, inters_template = start.get_estimation(
                100, iman.universe, s1_indices, s2_indices, cations, rings,
                cutoffs_aro, selected_aro, anions, hydrophobes, metal_donors,
                metal_acceptors, vdw_radii, max_vdw, hb_hydros, hb_donors,
                hb_acc, xb_halogens, xb_donors, xb_acc, cutoffs_others,
                selected_others)
            # %%

            to_allocate += ijf_template.shape[0] * ijf_template.shape[1]
            logger.debug(f"Number of allocated cells: {to_allocate}")

            # Compiling the parallel function
            logger.info("Compiling the parallel function")
            ijf, inters = start.run_parallel(
                xyz_chunk[:1], ijf_template, inters_template, len_others,
                len_aro, s1_indices, s2_indices, anions, cations, hydrophobes,
                metal_donors, metal_acceptors, vdw_radii, max_vdw, hb_hydros,
                hb_donors, hb_acc, xb_halogens, xb_donors, xb_acc, rings,
                cutoffs_others, selected_others, cutoffs_aro, selected_aro)

        # Compute the interactions
        logger.info(f"Computing chunk {i}")
        ijf_chunk, inters_chunk = start.run_parallel(
            xyz_chunk, ijf_template, inters_template, len_others, len_aro,
            s1_indices, s2_indices, anions, cations, hydrophobes, metal_donors,
            metal_acceptors, vdw_radii, max_vdw, hb_hydros, hb_donors, hb_acc,
            xb_halogens, xb_donors, xb_acc, rings, cutoffs_others,
            selected_others, cutoffs_aro, selected_aro)

        total += ijf_chunk.shape[0]
        total_inters += np.sum(inters_chunk)

        # Raise if not enough space has been allocated
        if (occupancy := ijf_chunk.shape[0] / to_allocate) >= 0.98:
            raise ValueError(f"Chunk {i} occupancy: {round(occupancy, 2)}")
        elif occupancy >= 0.90:
            logger.warning(f"Chunk {i} occupancy: {round(occupancy, 2)}")

        # %% Filling the interaction dictionary
        logger.info(f"Filling the interaction dictionary with chunk {i}")
        if ijf_chunk.shape[0] > 0:
            inter_dict.fill(ijf_chunk, inters_chunk)
            del ijf_chunk
            del inters_chunk
    #
    computing = round(time.time() - stamp, 2)

    # =========================================================================
    # Saving
    # =========================================================================
    inter_dict.pack()
    job_name = basename(args.job_name)
    pickle_name = f"{job_name}_InterMap.pickle"
    pickle_path = join(args.output_dir, pickle_name)
    gnl.pickle_to_file(inter_dict.dict, pickle_path)
    tot = round(time.time() - start_time, 2)
    logger.info(f"Normal termination of InterMap job '{job_name}' in {tot} s")
    n_ints = len(inter_dict.dict.keys())
    logger.info(f"Inter_dict size: {n_ints}")
    logger.info(
        f"Total interactions detected in {computing} s: {total}, {total_inters}")

# run(mode='debug')
