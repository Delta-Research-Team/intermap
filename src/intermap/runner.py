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
from intermap import aro, commons as cmn, geometry as aot, macro
from intermap.indices import IndexManager


# %% todo: check docstrings


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
    metal_donors, metal_acceptors = iman.metal_don, iman.metal_acc
    hb_hydros, hb_donors, hb_acc = iman.hb_H, iman.hb_D, iman.hb_A
    xb_halogens, xb_donors, xb_acc = iman.xb_H, iman.xb_D, iman.xb_A

    anions, cations = iman.anions, iman.cations
    rings = iman.rings
    s1_cat = s1_indices[aot.isin(s1_indices, cations)]
    s2_cat = s2_indices[aot.isin(s2_indices, cations)]
    s1_rings = rings[aot.isin(rings[:, 0], s1_indices)]
    s2_rings = rings[aot.isin(rings[:, 0], s2_indices)]
    # Internal indexing for xyz2 coordinates
    n0 = s1_cat.size + s2_cat.size
    n1 = n0 + s1_rings.shape[0]
    n2 = n1 + s2_rings.shape[0]
    s1_cat_idx = np.arange(0, s1_cat.size, dtype=np.int32)
    s2_cat_idx = np.arange(s1_cat.size, n0, dtype=np.int32)
    s1_rings_idx = np.arange(n0, n1, dtype=np.int32)
    s2_rings_idx = np.arange(n1, n2, dtype=np.int32)
    s1_aro_indices = np.concatenate((s1_cat_idx, s1_rings_idx)).astype(
        np.int32)
    s2_aro_indices = np.concatenate((s2_cat_idx, s2_rings_idx)).astype(
        np.int32)
    xyz_aro_real_idx = np.concatenate(
        (s1_cat, s2_cat, s1_rings[:, 0], s2_rings[:, 0])).astype(np.int32)

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

    # todo: put this in args
    dist_cut_aro = cutoffs_aro[:2].max()
    dist_cut_others = max(cutoffs_others[:2].max(), max_vdw)

    cutoffs_str = {x: args.cutoffs[x] for x in args.cutoffs if x in to_compute}
    logger.info(f"Interactions to compute:\n {pformat(to_compute)}")
    logger.debug(f"Cutoffs parsed:\n {pformat(cutoffs_str)}")

    # =========================================================================
    # Estimating memory allocation
    # =========================================================================
    logger.info(f"Estimating memory allocation")
    set_num_threads(args.n_procs)
    n_frames, n_samples = len(universe.trajectory), 10
    sub = universe.trajectory[::n_frames // n_samples]
    f4 = np.float32
    positions = np.asarray([ts.positions.copy() for ts in sub], dtype=f4)

    ijf_shape, inters_shape, mb1, mb2, v_size, h_size = macro.estimate(
        positions, xyz_aro_real_idx, args.chunk_size, s1_indices, s2_indices,
        cations, s1_cat_idx, s2_cat_idx, s1_cat, s2_cat, s1_rings, s2_rings,
        s1_rings_idx, s2_rings_idx, s1_aro_indices, s2_aro_indices,
        cutoffs_aro, selected_aro, len_aro, anions, hydrophobes, metal_donors,
        metal_acceptors, vdw_radii, max_vdw, hb_hydros, hb_donors, hb_acc,
        xb_halogens, xb_donors, xb_acc, cutoffs_others, selected_others,
        len_others, dist_cut_aro)

    logger.debug(f"Allocated space for interactions:"
                 f" ~{mb1} MB ({args.chunk_size}, {v_size}, {h_size})")
    logger.debug(f"Allocated space for coordinates:"
                 f" ~{mb2} MB ({args.chunk_size}, {sel_idx.size}, 3) ")

    # %%=======================================================================
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
    trajectory = universe.trajectory
    # =========================================================================
    xyz_chunk = None
    trees_chunk = None
    s1_centrs, s2_centrs, xyzs_aro = None, None, None
    aro_balls = None
    ijf_chunk, inters_chunk = None, None
    # =========================================================================
    for i, frames in tqdm(enumerate(chunks), total=N,
                          desc='Detecting Interactions', unit='chunk'):

        xyz_chunk = tt.get_coordinates(trajectory, frames, sel_idx)

        trees_chunk = cmn.get_trees(xyz_chunk, s2_indices)
        s1_centrs, s2_centrs, xyzs_aro = aro.get_aro_xyzs(
            xyz_chunk, s1_rings, s2_rings, s1_cat, s2_cat)
        aro_balls = aro.get_balls(
            xyzs_aro, s1_aro_indices, s2_aro_indices, dist_cut_aro)

        ijf_chunk, inters_chunk = macro.runpar(
            xyz_chunk, xyzs_aro, xyz_aro_real_idx, trees_chunk, aro_balls,
            ijf_shape, inters_shape, len_others, len_aro, s1_indices,
            s2_indices, anions, cations, s1_cat_idx, s2_cat_idx, hydrophobes,
            metal_donors, metal_acceptors, vdw_radii, max_vdw, hb_hydros,
            hb_donors, hb_acc, xb_halogens, xb_donors, xb_acc, s1_rings,
            s2_rings, s1_rings_idx, s2_rings_idx, s1_aro_indices,
            s2_aro_indices, cutoffs_others,
            selected_others, cutoffs_aro, selected_aro)

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
