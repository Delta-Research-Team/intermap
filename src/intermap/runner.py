# Created by gonzalezroy at 11/14/24

import mdtraj as md
import numpy as np
from numba import set_num_threads

import intermap.finders as find
import intermap.topo_trajs as tt
from intermap.indexman import IndexManager as iman

# todo: Parallelization of numba functions could be made by preallocating the
#  array size for the interactions in the frame. Then create a List of arrays
#  and append the arrays to the list.



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
    prev_cutoff = args.prev_cutoff

    # Get the selections
    master_traj = md.load_frame(traj, top=topo, index=0)
    index_manager = iman(sel1, sel2, master_traj, interactions)
    inter_types = np.asarray(['hb', 'cc'])
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
    _ = find.inters(master_traj.xyz[:1], 0,
                    s1_donors, s1_hydros, s1_acc,
                    s2_donors, s2_hydros, s2_acc,
                    s1_indices, s2_indices)

    n_frame = 0
    intermap_lists = []
    for chunk in chunks:
        # Get the chunk coords and mark the chunk object to free memory
        xyz = chunk.xyz

        # # Run the interactions in parallel
        inter_list = find.inters(xyz, n_frame,
                                 s1_donors, s1_hydros, s1_acc,
                                 s2_donors, s2_hydros, s2_acc,
                                 s1_indices, s2_indices)
        intermap_lists.extend(inter_list)
        n_frame += len(xyz)

    intermap_lists = np.concatenate(intermap_lists)

    if intermap_lists.size == 0:
        return None

    intermap = find.get_intermap(intermap_lists, labels, inter_types, n_frame,
                                 prev_cutoff=prev_cutoff)
    return intermap


# =============================================================================
# Debugging area
# =============================================================================
# from argparse import Namespace
#
# topo = '/home/gonzalezroy/RoyHub/oxo-8/data/raw/oligo/A1/A1_21bp_box_dry.prmtop'
# traj = '/home/gonzalezroy/RoyHub/oxo-8/data/raw/oligo/A1/8oxoGA1_21bp_1_dry.nc'
# sel1 = "(resname =~ '(5|3)?D([ATGC])|(8OG){1}(3|5)?$')"
# sel2 = sel1
# nprocs = 8
# chunk_size = 500
# interactions = ['hbonds', 'close_contacts']
# prev_cutoff = 2
# args = Namespace(topo=topo, traj=traj, sel1=sel1, sel2=sel2, nprocs=nprocs,
#                  chunk_size=chunk_size, interactions=interactions, start=0,
#                  last=-1, stride=1, prev_cutoff=prev_cutoff)
# inter_map = run(args)
