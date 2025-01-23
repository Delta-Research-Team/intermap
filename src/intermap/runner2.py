# Created by rglez at 12/29/24
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
import intermap.interactions as its
import intermap.interdict as idt
import intermap.topo_trajs as tt
from intermap import commons as cmn
from intermap.indices import IndexManager


# todo: remove prolif dependency by copying the SMARTS patterns
# todo: check docstrings
# todo: come up with an idea for loading just the necessary indices into memory
# todo: put a logger when potentially time-consuming tasks starts

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
        config_path = '/home/rglez/RoyHub/intermap/tests/imap_lig-prot.cfg'
        config_path = '/home/rglez/RoyHub/intermap/tests/imap_prot-prot.cfg'
        # config_path = '/home/gonzalezroy/RoyHub/intermap/example/imap.cfg'
    else:
        raise ValueError('Only modes allowed are production and running')

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
    # Parsing the interactions & cutoffs
    # =========================================================================
    all_inters, all_cutoffs = cf.get_inters_cutoffs(args.cutoffs)
    if isinstance(args.interactions, str) and args.interactions == 'all':
        to_compute = all_inters
    else:
        to_compute = args.interactions
    selected_aro, selected_others, cutoffs_aro, cutoffs_others = \
        cmn.get_cutoffs_and_inters(to_compute, all_inters, all_cutoffs)

    len_others, len_aro = len(selected_others), len(selected_aro)
    cutoffs_str = {x: args.cutoffs[x] for x in args.cutoffs if x in to_compute}
    logger.info(f"Interactions to compute:\n {pformat(to_compute)}")
    logger.debug(f"Cutoffs parsed:\n {pformat(cutoffs_str)}")

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
    hb_hydrogens, hb_donors, hb_acceptors = iman.hb_H, iman.hb_D, iman.hb_A
    xb_halogens, xb_donors, xb_acceptors = iman.xb_H, iman.xb_D, iman.xb_A
    rings = iman.rings

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
    atnames = universe.atoms.names[sel_idx]
    resnames = universe.atoms.resnames[sel_idx]
    resids = universe.atoms.resids[sel_idx]
    n_sel_atoms = sel_idx.size
    names = [f"{resnames[i]}_{resids[i]}_{atnames[i]}" for i in sel_idx]
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
    to_allocate = 0
    for i, frames_chunk in enumerate(chunks):
        xyz_chunk = tt.get_coordinates(universe, frames_chunk, sel_idx,
                                       n_sel_atoms)

        if i == 0:
            # Compiling main functions
            logger.info("Compiling main functions")

            # Estimating memory allocation
            logger.info(f"Estimating memory allocation")
            ijf_template, inters_template = its.get_estimation(
                xyz_chunk, 5, s1_indices, s2_indices, cations, rings,
                cutoffs_aro, selected_aro, anions, hydrophobes, metal_donors,
                metal_acceptors, vdw_radii, max_vdw, hb_hydrogens, hb_donors,
                hb_acceptors, xb_halogens, xb_donors, xb_acceptors,
                cutoffs_others, selected_others)

            # Compiling the parallel function
            _, _ = its.run_parallel(
                xyz_chunk[:1], ijf_template, inters_template, len_others,
                len_aro, s1_indices, s2_indices, anions, cations, hydrophobes,
                metal_donors, metal_acceptors, vdw_radii, max_vdw,
                hb_hydrogens,
                hb_donors, hb_acceptors, xb_halogens, xb_donors, xb_acceptors,
                rings, cutoffs_others, selected_others, cutoffs_aro,
                selected_aro)

            to_allocate = ijf_template.shape[0] * ijf_template.shape[1]
            logger.debug(f"Number of allocated cells: {to_allocate}")

        logger.info(f"Computing chunk {i}")
        ijf_chunk, inters_chunk = its.run_parallel(
            xyz_chunk, ijf_template, inters_template, len_others, len_aro,
            s1_indices, s2_indices, anions, cations, hydrophobes,
            metal_donors, metal_acceptors, vdw_radii, max_vdw, hb_hydrogens,
            hb_donors,
            hb_acceptors, xb_halogens, xb_donors, xb_acceptors, rings,
            cutoffs_others, selected_others, cutoffs_aro, selected_aro)

        # Raise if not enough space has been allocated
        if (occupancy := ijf_chunk.shape[0] / to_allocate) >= 0.98:
            raise ValueError(f"Chunk {i} occupancy: {round(occupancy, 2)}")
        elif occupancy >= 0.90:
            logger.warning(f"Chunk {i} occupancy: {round(occupancy, 2)}")

        # %% Filling the interaction dictionary
        logger.info(f"Filling the interaction dictionary")
        inter_dict.fill(ijf_chunk, inters_chunk)
    #
    computing = round(time.time() - stamp, 2)
    n_ints = len(inter_dict.dict)
    logger.info(f"Total interactions detected in {computing} s: {n_ints}")

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


run(mode='debug')
