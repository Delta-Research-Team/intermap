# Created by gonzalezroy at 11/14/24
from collections import defaultdict
from multiprocessing import Pool

import mdtraj as md
from numpy.core.defchararray import strip

import intermap.topo_trajs as tt
from intermap.indexman import IndexManager as iman
from intermap.interactions import find_hbonds


def run(args):
    """
    Entry point for running the intermap workflow
    """
    # Parse the arguments
    traj = args.traj
    topo = args.topo
    sel1 = args.sel1
    sel2 = args.sel2
    start = args.start
    last = args.last
    stride = args.stride
    chunk_size = args.chunk_size
    nprocs = args.nprocs

    # Load the topology and the trajectory
    master_traj = md.load_frame(traj, top=topo, index=0)

    # Get the selections
    seles = iman(sel1, sel2, master_traj)
    s1_donors, s1_hydros, s1_acc = seles.s1_donors, seles.s1_hydros, seles.s1_acc
    s2_donors, s2_hydros, s2_acc = seles.s2_donors, seles.s2_hydros, seles.s2_acc

    # Yield chunks of the trajectory
    chunks = tt.get_traj_chunks(topo, [traj],
                                start=start, last=last, stride=stride,
                                chunk_size=chunk_size)

    # Find the hydrogen bonds in parallel for each chunk
    hbonds = defaultdict(list)
    frame = 0
    for chunk in chunks:

        # Get the chunk coords and delete the chunk object to free memory
        xyz = chunk.xyz
        del chunk

        # Prepare the arguments; using a generator for delayed evaluation
        arguments = (
            (coord, s1_donors, s1_hydros, s1_acc, s2_donors, s2_hydros, s2_acc)
            for coord in xyz)

        # Run the function in parallel
        with Pool(nprocs) as p:
            hbonds_chunk = p.starmap(find_hbonds, arguments)
            del xyz

            # Convert the results to a dictionary of counts per frame
            for i, hbonds_frame in enumerate(hbonds_chunk):
                for hbond in hbonds_frame:
                    hbonds[tuple(hbond)].append(i + frame)

        # Update the frame counter
        frame += len(hbonds_chunk)
    return seles, hbonds


# =============================================================================
# Debugging area
# =============================================================================
# from argparse import Namespace
#
# topo = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A1/8oxoGA1_1_hmr.prmtop'
# traj = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A1/8oxoGA1_1_sk100.nc'
# sel1 = "(resname =~ '(5|3)?D([ATGC])|(8OG){1}(3|5)?$')"
# sel2 = "water"
# nprocs = 16
# chunk_size = 100
# args = Namespace(topo=topo, traj=traj, sel1=sel1, sel2=sel2, nprocs=nprocs,
#                  chunk_size=chunk_size)
