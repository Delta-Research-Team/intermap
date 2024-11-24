# Created by gonzalezroy at 11/14/24
from collections import defaultdict


import mdtraj as md
from numba import set_num_threads

import intermap.finders as find
import intermap.topo_trajs as tt
from intermap.indexman import IndexManager as iman



def run(args):
    """
    Entry point for running the intermap workflow
    """
    # Parse the arguments
    traj = args.traj
    topo = args.topo
    sel1 = args.sel1
    sel2 = args.sel2
    interactions = args.interactions
    chunk_size = args.chunk_size
    nprocs = args.nprocs
    start = args.start
    last = args.last
    stride = args.stride

    # Get the selections
    master_traj = md.load_frame(traj, top=topo, index=0)
    index_manager = iman(sel1, sel2, master_traj, interactions)
    labels = index_manager.labels

    # ... for close contacts
    s1_indices = index_manager.s1_idx
    s2_indices = index_manager.s2_idx

    # ... for hydrogen bonds
    selections = index_manager.indices
    s1_hydros = selections['hbonds']['s1_hydros']
    s1_donors = selections['hbonds']['s1_donors']
    s1_acc = selections['hbonds']['s1_acc']
    s2_hydros = selections['hbonds']['s2_hydros']
    s2_donors = selections['hbonds']['s2_donors']
    s2_acc = selections['hbonds']['s2_acc']

    # Yield chunks of the trajectory and find the interactions in parallel
    chunks = tt.get_traj_chunks(topo, [traj],
                                start=start, last=last, stride=stride,
                                chunk_size=chunk_size)

    # Compile the interactions computing functions & set the number of threads
    set_num_threads(nprocs)
    _ = find.inters(master_traj.xyz[:1],
                    s1_donors, s1_hydros, s1_acc,
                    s2_donors, s2_hydros, s2_acc,
                    s1_indices, s2_indices)

    frame = 0
    intermap_dict = defaultdict(list)
    for chunk in chunks:
        # Get the chunk coords and mark the chunk object to free memory
        xyz = chunk.xyz
        N = len(xyz)
        # del chunk

        # Run the interactions in parallel
        hb_list, cc_list = find.inters(xyz,
                                       s1_donors, s1_hydros, s1_acc,
                                       s2_donors, s2_hydros, s2_acc,
                                       s1_indices, s2_indices)
        # del xyz

        # Update the intermap dictionary
        intermap_dict = find.update_intermap(intermap_dict, hb_list, 'hb',
                                             labels, frame)
        # del hb_list
        intermap_dict = find.update_intermap(intermap_dict, cc_list, 'cc',
                                             labels, frame)
        # del cc_list

        # Update the frame counter
        frame += N
    intermap = find.intermap2df(intermap_dict, frame)
    return intermap


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
# chunk_size = 100
# interactions = ['hbonds', 'close_contacts']
#
# args = Namespace(topo=topo, traj=traj, sel1=sel1, sel2=sel2, nprocs=nprocs,
#                  chunk_size=chunk_size, interactions=interactions)
# intermap = run(args)
