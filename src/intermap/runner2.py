# Created by rglez at 12/29/24
import time
from argparse import Namespace
from os.path import dirname, join
from pprint import pformat

import numpy as np
import rgpack.generals as gnl
from numba import set_num_threads

import intermap.cutoffs as cf
import intermap.interactions as its
import intermap.interdict as idt
import intermap.topo_trajs as tt
from intermap import commons as cmn
from intermap.indices import IndexManager

# todo: parser, detect interactions
# todo: parser, detect cutoffs
# todo: remove prolif dependency by copying the SMARTS patterns

start_time = time.time()
# %% ==========================================================================
# Parsed arguments
# =============================================================================
args = Namespace()
args.topo = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1.prmtop'
args.traj = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1_sk100.nc'
args.sel1 = "nucleic"
args.sel2 = "protein"
args.interactions = 'all'
args.n_procs = 8
args.job_name = '8oxoGA2_1'
args.start = 0
args.last = 150
args.stride = 1
args.chunk_size = 50

args.format = 'simple'
args.min_prevalence = 1.0

args.dist_cut_CloseContacts = 3.0
args.dist_cut_Ionic = 4.5
args.dist_cut_Hydrophobic = 4.5
args.dist_cut_Metalic = 2.8
args.dist_cut_HA = 2.5
args.dist_cut_DA = 3.9
args.min_ang_DHA = 130
args.max_ang_DHA = 180
args.dist_cut_XA = 2.5
args.dist_cut_XD = 3.9
args.ang_cut_DXA = 90
args.dist_cut_PiCation = 4.5
args.min_ang_PiCation = 0
args.max_ang_PiCation = 30
args.dist_cut_PiStacking = 6.0
args.min_dist_PiStacking = 3.8
args.min_ang_PiStacking = 0
args.max_ang_PiStacking = 90
args.dist_cut_EdgeToFace = 6.0
args.min_dist_EdgeToFace = 3.8
args.min_ang_EdgeToFace = 50
args.max_ang_EdgeToFace = 90
args.dist_cut_FaceToFace = 4.5
args.min_dist_FaceToFace = 3.8
args.min_ang_FaceToFace = 0
args.max_ang_FaceToFace = 40

# =============================================================================
# Starting logger
# =============================================================================
log_path = join(dirname(args.traj), f"{args.job_name}_InterMap.log")
logger = cmn.start_logger(log_path)
logger.info(f"Starting InterMap with the following parameters:"
            f"\n Topology: {args.topo}"
            f"\n Trajectory: {args.traj}"
            f"\n Start frame: {args.start}"
            f"\n Last frame: {args.last}"
            f"\n Stride: {args.stride}"
            f"\n Selection 1: {args.sel1}"
            f"\n Selection 2: {args.sel2}"
            f"\n Requested Interactions: {args.interactions}"
            f"\n Number of processors: {args.n_procs}"
            f"\n Job name: {args.job_name}")

# =============================================================================
# Parsing the interactions & cutoffs
# =============================================================================

parsed_cutoffs = cf.parse_cutoffs(args=args)
all_inters, all_cutoffs = cf.get_inters_cutoffs(parsed_cutoffs)
if args.interactions == 'all':
    to_compute = all_inters
else:
    to_compute = args.interactions

selected_aro, selected_others, cutoffs_aro, cutoffs_others = \
    cmn.get_cutoffs_and_inters(to_compute, all_inters, all_cutoffs)

cutoffs_str = {x: parsed_cutoffs[x] for x in parsed_cutoffs if x in to_compute}
logger.info(f"Interactions to compute:\n {pformat(to_compute)}")
logger.debug(f"Cutoffs parsed:\n {pformat(cutoffs_str)}")

set_num_threads(args.n_procs)

# =============================================================================
# Load the necessary indices for detecting each interaction
# =============================================================================
iman = IndexManager(args.topo, args.traj, args.sel1, args.sel2,
                    args.interactions)
u = iman.universe
n_atoms = iman.n_atoms
s1_indices_raw, s2_indices_raw = iman.sel1_idx, iman.sel2_idx
vdw_radii = iman.radii
hydrophobes = iman.hydroph
anions, cations = iman.anions, iman.cations
metal_donors, metal_acceptors = iman.metal_don, iman.metal_acc
hb_hydrogens, hb_donors, hb_acceptors = iman.hb_H, iman.hb_D, iman.hb_A
xb_halogens, xb_donors, xb_acceptors = iman.xb_H, iman.xb_D, iman.xb_A
rings = iman.rings

# =============================================================================
# Load & trim the trajectory
# =============================================================================
n_frames = iman.n_frames
last = tt.parse_last_param(args.last, n_frames)
traj_frames = np.arange(args.start, last, args.stride)
logger.info(f"Number of frames to consider (start:last:stride): "
            f"{traj_frames.size} ({args.start}:{last}:{args.stride})")

# =============================================================================
# Naming all atoms and interactions
# =============================================================================
atnames = u.atoms.names
resnames = u.atoms.resnames
resids = u.atoms.resids
names = [f"{resnames[i]}_{resids[i]}_{atnames[i]}" for i in range(n_atoms)]
inters = np.asarray(selected_others.tolist() + selected_aro.tolist())

# %% ==========================================================================
# Fill the interaction dictionary
# =============================================================================
logger.info(f"Starting to compute InterMap interactions")
stamp = time.time()

inter_dict = idt.InterDict(args.format, args.min_prevalence, traj_frames,
                           names, inters)

chunks = tt.split_in_chunks(traj_frames, args.chunk_size)
len_others, len_aro = len(selected_others), len(selected_aro)
for i, frames_chunk in enumerate(chunks):

    # Load the coordinates for the chunk
    xyz_chunk = tt.get_coordinates(u, frames_chunk, n_atoms)

    # Compiling the functions & estimating the interactions
    if i == 0:
        stamp = time.time()

        ijf_template, inters_template = its.get_estimation(
            xyz_chunk, 5, s1_indices_raw, s2_indices_raw, cations, rings,
            cutoffs_aro, selected_aro, anions, hydrophobes, metal_donors,
            metal_acceptors, vdw_radii, hb_hydrogens, hb_donors, hb_acceptors,
            xb_halogens, xb_donors, xb_acceptors, cutoffs_others,
            selected_others)

        its.run_parallel(
            xyz_chunk[:1], ijf_template[:1], inters_template[:1], len_others,
            len_aro, s1_indices_raw[:1], s2_indices_raw[:1], anions,
            cations, hydrophobes, metal_donors,
            metal_acceptors, vdw_radii, hb_hydrogens,
            hb_donors, hb_acceptors, xb_halogens, xb_donors,
            xb_acceptors, rings, cutoffs_others, selected_others,
            cutoffs_aro, selected_aro)

        comp = round(time.time() - stamp, 2)
        logger.info(f"Elapsed time compiling: {comp:.2f} s")
        max_allocated = ijf_template.shape[0] * ijf_template.shape[1]
        t1 = round(time.time() - start_time, 2)
        logger.info(f'Elapsed time before start computing interactions: {t1}')

    # Parallel computing of the interactions
    ijf_chunk, inters_chunk = its.run_parallel(
        xyz_chunk, ijf_template, inters_template, len_others, len_aro,
        s1_indices_raw, s2_indices_raw, anions, cations, hydrophobes,
        metal_donors, metal_acceptors, vdw_radii, hb_hydrogens, hb_donors,
        hb_acceptors, xb_halogens, xb_donors, xb_acceptors, rings,
        cutoffs_others, selected_others, cutoffs_aro, selected_aro)

    # Raise if not enogh space has been allocated
    if (occupancy := ijf_chunk.shape[0] / max_allocated) >= 0.98:
        raise ValueError(f"Chunk {i} occupancy: {round(occupancy, 2)}")
    elif occupancy >= 0.90:
        logger.error(f"Chunk {i} occupancy: {round(occupancy, 2)}")

    # Filling the interaction dictionary
    ijf_chunk[:, 2] = frames_chunk[ijf_chunk[:, 2]]
    inter_dict.fill(ijf_chunk, inters_chunk)
computing = round(time.time() - stamp, 2)
logger.info(f"Interactions computed in {computing} s")

# =============================================================================
# Saving
# =============================================================================
logger.info(f"Saving InterMap dictionary")
inter_dict.pack()
pickle_path = f"{args.job_name}_InterMap_Dict.pkl"
gnl.pickle_to_file(inter_dict.dict, pickle_path)

logger.info(f"Normal termination of InterMap job '{args.job_name}' "
            f"in {round(time.time() - start_time, 2)} s")


unpick = gnl.unpickle_from_file(pickle_path)
