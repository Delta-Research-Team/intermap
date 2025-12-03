# Created by rglez at 11/23/25
"""
In this script, we concatenate multiple trajectory files into a single
 trajectory filefor each molecular dynamics simulation system.
"""
from os.path import basename

import MDAnalysis as mda
import numpy as np
from rgpack import generals as gnl


def concatextract(trajs, topo, out_dcd, n_frames=1000):
    """Concatenate multiple trajectory files into a single trajectory file,
    extracting a specified number of frames evenly distributed across all
    input trajectories.

    Parameters
    ----------
    trajs : list of str
        List of paths to the input trajectory files.
    topo : str
        Path to the topology file.
    out_dcd : str
        Path to the output concatenated trajectory file.
    n_frames : int, optional
        Number of frames to extract and concatenate, by default 1000.
    """
    sizes = {}
    universes = {}
    for trj in trajs:
        u = mda.Universe(topo, trj)
        size = len(u.trajectory)
        sizes[basename(trj)] = size
        universes[basename(trj)] = u

    trj_names = np.asarray([z for x, y in sizes.items() for z in [x] * y])
    trj_indices = np.asarray([i for x, y in sizes.items() for i in range(y)])

    N = sum([x for x in sizes.values()])
    selected = np.asarray(range(0, N, N // 1000))[:1000]
    sub_names = trj_names[selected]
    sub_idx = trj_indices[selected]

    to_extract = {}
    for name, idx in zip(sub_names, sub_idx):
        if name not in to_extract:
            to_extract[name] = []
        to_extract[name].append(idx)

    example_u = universes[sub_names[0]]
    with mda.Writer(out_dcd, n_atoms=example_u.atoms.n_atoms) as W:

        for fname in sizes.keys():  # preserve DCD order
            if fname not in to_extract:
                continue
            u = universes[fname]
            frame_list = to_extract[fname]

            print(f"Extracting {len(frame_list)} frames from {fname}")

            for f in frame_list:
                u.trajectory[f]
                W.write(u.atoms)


# =============================================================================
# ACE2 Trajectories Concatenation
# =============================================================================
# ace_dir = '/media/rglez/Roy5T/RoyData/intermap/trajs/ACE2'
# out_dcd = '/media/rglez/Roy5T/RoyData/intermap/trajs/ACE2/ace2_concat.dcd'
# topo = next(gnl.recursive_finder('*.psf', ace_dir))
# trajs = sorted(gnl.recursive_finder('*.dcd', ace_dir))
# concatextract(trajs, topo, out_dcd)

# =============================================================================
# Ec_T4P Trajectories Concatenation
# =============================================================================
# pao_dir = '/media/rglez/Roy5T/RoyData/intermap/trajs/Ec_T4P'
# out_dcd = '/media/rglez/Roy5T/RoyData/intermap/trajs/Ec_T4P/ec_t4p_concat.dcd'
# topo = next(gnl.recursive_finder('*.pdb', pao_dir))
# trajs = sorted(gnl.recursive_finder('*.xtc', pao_dir))
# concatextract(trajs, topo, out_dcd)

# =============================================================================
# IgG3 Trajectories Concatenation
# =============================================================================
# _dir = '/media/rglez/Roy5T/RoyData/intermap/trajs/IgG3-M1'
# out_dcd = '/media/rglez/Roy5T/RoyData/intermap/trajs/IgG3-M1/igg3_m1_concat.dcd'
# topo = next(gnl.recursive_finder('*.pdb', _dir))
# trajs = sorted(gnl.recursive_finder('*.xtc', _dir))
# concatextract(trajs, topo, out_dcd)

# =============================================================================
# NSP13 Trajectories Concatenation
# =============================================================================
# _dir = '/media/rglez/Roy5T/RoyData/intermap/trajs/DESRES-Trajectory_sarscov2-13795965-no-water/sarscov2-13795965-no-water'
# out_dcd = '/media/rglez/Roy5T/RoyData/intermap/trajs/DESRES-Trajectory_sarscov2-13795965-no-water/sarscov2-13795965-no-water/sarscov2_concat.dcd'
# topo = next(gnl.recursive_finder('*.psf', _dir))
# trajs = sorted(gnl.recursive_finder('*.dcd', _dir))
# concatextract(trajs, topo, out_dcd)

# =============================================================================
# mpro Trajectories Concatenation
# =============================================================================
# _dir = '/media/rglez/Roy5T/RoyData/intermap/1000Frames/mpro'
# out_dcd = '/media/rglez/Roy5T/RoyData/intermap/1000Frames/mpro/mpro_concat.dcd'
# topo = next(gnl.recursive_finder('*.pdb', _dir))
# trajs = sorted(gnl.recursive_finder('*.dcd', _dir))
# concatextract(trajs, topo, out_dcd)
