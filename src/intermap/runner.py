# Created by gonzalezroy at 11/14/24
from multiprocessing import Pool

import mdtraj as md

import intermap.topo_trajs as tt
from intermap.indexman import IndexManager as iman
from intermap.interactions import find_hbonds


def run(args):
    """
    Entry point for running the intermap workflow
    """
    traj = args.traj
    topo = args.topo
    sel1 = args.sel1
    sel2 = args.sel2
    chunk_size = args.chunk_size
    nprocs = args.nprocs

    # Load the topology and the trajectory
    master_traj = md.load_frame(traj, top=topo, index=0)

    seles = iman(sel1, sel2, master_traj)
    s1_donors, s1_hydros, s1_acc = seles.s1_donors, seles.s1_hydros, seles.s1_acc
    s2_donors, s2_hydros, s2_acc = seles.s2_donors, seles.s2_hydros, seles.s2_acc

    # Define the function to be parallelized
    hbonds = []
    chunks = tt.get_traj_chunks(topo, [traj], chunk_size=chunk_size)
    for chunk in chunks:
        xyz = chunk.xyz
        with Pool(nprocs) as p:
            arguments = []
            for frame in xyz:
                arguments.append((frame, s1_donors, s1_hydros, s1_acc,
                                  s2_donors, s2_hydros, s2_acc))
            del xyz
            results = p.starmap(find_hbonds, arguments)
            hbonds.extend(results)

    return hbonds

# =============================================================================
# Debugging area
# =============================================================================
# from argparse import Namespace
#
# topo = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A1/8oxoGA1_1_hmr.prmtop'
# traj = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A1/8oxoGA1_1_sk100.nc'
# sel1 = "(resname =~ '(5|3)?D([ATGC])|(8OG){1}(3|5)?$')"
# sel2 = "water"
# nprocs = 8
# chunk_size = 150
# args = Namespace(topo=topo, traj=traj, sel1=sel1, sel2=sel2, nprocs=nprocs, chunk_size=chunk_size)
