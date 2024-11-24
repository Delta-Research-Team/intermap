# Created by rglez at 11/16/24
"""
Functions to find the interactions between two selections in a frame.
"""

import numpy as np
import pandas as pd
from numba import njit, prange
from numba.typed import List
from numba_kdtree import KDTree as nckd

from intermap.geometric import calc_angle, calc_dist


@njit(parallel=False)
def hbonds(xyz, k, s1_donors, s1_hydros, s1_acc, s2_donors, s2_hydros, s2_acc):
    """
    Args:
        xyz: Coordinates of the atoms in the frame.
        k (int): Frame identifier.
        s1_donors: Indices of the donor atoms in s1.
        s1_hydros: Indices of the hydrogen atoms in s1.
        s1_acc: Indices of the acceptor atoms in s1.
        s2_donors: Indices of the donor atoms in s2.
        s2_hydros: Indices of the hydrogen atoms in s2.
        s2_acc: Indices of the acceptor atoms in s2.

    Returns:
        DHA_labels (list): List of DHA labels for each hydrogen bond.
    """
    # Hardcoded cutoffs
    HA_cut = 0.25  # distance between Hydrogen and Acceptor cutoff
    DA_cut = 0.39  # distance between Donor and Acceptor cutoff
    DHA_cut = 90  # angle between Donor, Hydrogen and Acceptor cutoff

    # Create & query the trees
    s2_A_tree = nckd(xyz[s2_acc])
    s1_A_tree = nckd(xyz[s1_acc])

    ball_1 = s2_A_tree.query_radius(xyz[s1_hydros], HA_cut)
    ball_2 = s1_A_tree.query_radius(xyz[s2_hydros], HA_cut)

    # Find the DHA candidates
    n_triples = 0
    for x in ball_1:
        n_triples += x.size
    for x in ball_2:
        n_triples += x.size
    DHA_labels = np.empty((n_triples, 3), dtype=np.int32)
    DHA_coords = np.empty((n_triples, 3, 3), dtype=np.float32)

    counter = -1
    for i in range(len(ball_1)):
        D = s1_donors[i]
        H = s1_hydros[i]
        for j in ball_1[i]:
            A = s2_acc[j]
            counter += 1
            DHA_labels[counter] = [D, A, k]
            DHA_coords[counter][0] = xyz[D]
            DHA_coords[counter][1] = xyz[H]
            DHA_coords[counter][2] = xyz[A]

    for i in range(len(ball_2)):
        D = s2_donors[i]
        H = s2_hydros[i]
        for j in ball_2[i]:
            A = s1_acc[j]
            counter += 1
            DHA_labels[counter] = [D, A, k]
            DHA_coords[counter][0] = xyz[D]
            DHA_coords[counter][1] = xyz[H]
            DHA_coords[counter][2] = xyz[A]

    # Filter the DHA candidates
    DA_dist_correct = calc_dist(DHA_coords[:, 0],
                                DHA_coords[:, 2]) <= DA_cut
    DHA_angle_correct = calc_angle(DHA_coords[:, 0, :],
                                   DHA_coords[:, 1, :],
                                   DHA_coords[:, 2, :]) > DHA_cut

    correct_DHA = DHA_labels[DA_dist_correct & DHA_angle_correct]

    return correct_DHA


@njit(parallel=False)
def close_contacts(xyz, k, s1_indices, s2_indices):
    """
    Args:
        xyz (ndarray): Coordinates of the atoms in the frame.
        k (int): Frame identifier.
        s1_indices (ndarray): Indices of the atoms in s1.
        s2_indices (ndarray): Indices of the atoms in s2.

    Returns:
        CC_labels (list): List of labels for each close contact.
    """
    # Hardcoded cutoffs
    CC_cut = 0.3

    # Create & query the trees
    s2_tree = nckd(xyz[s2_indices])
    ball_1 = s2_tree.query_radius(xyz[s1_indices], CC_cut)

    # Find the CC
    n_doubles = sum([len(x) for x in ball_1])
    CC_labels = np.zeros((n_doubles, 3), dtype=np.int32)

    counter = -1
    for i, x in enumerate(ball_1):
        X1 = s1_indices[i]
        for j in x:
            X2 = s2_indices[j]
            counter += 1
            CC_labels[counter][0] = X1
            CC_labels[counter][1] = X2
            CC_labels[counter][2] = k
    return CC_labels


@njit(parallel=True)
def inters(xyz, s1_donors, s1_hydros, s1_acc, s2_donors, s2_hydros, s2_acc,
           s1_indices, s2_indices):
    """
    Find the interactions between two selections in a frame.

    Args:
        xyz (ndarray): Coordinates of the atoms in the frame.
        s1_donors (ndarray): Indices of the donor atoms in s1.
        s1_hydros (ndarray): Indices of the hydrogen atoms in s1.
        s1_acc (ndarray): Indices of the acceptor atoms in s1.
        s2_donors (ndarray): Indices of the donor atoms in s2.
        s2_hydros (ndarray): Indices of the hydrogen atoms in s2.
        s2_acc (ndarray): Indices of the acceptor atoms in s2.
        s1_indices (ndarray): Indices of the atoms in s1.
        s2_indices (ndarray): Indices of the atoms in s2.

    Returns:
        hb_list (List): List of hydrogen bonds for each frame.
        cc_list (List): List of close contacts for each frame.
    """
    num_frames = xyz.shape[0]
    hb_list = List()
    cc_list = List()
    for i in prange(num_frames):
        hb_list.append(
            hbonds(
                xyz[i], i,
                s1_donors, s1_hydros, s1_acc, s2_donors,
                s2_hydros, s2_acc))
        cc_list.append(
            close_contacts(
                xyz[i], i, s1_indices, s2_indices))
    return hb_list, cc_list


def update_intermap(intermap_dict, inter_list, inter_type, labels, frame):
    """
    Update the intermap dictionary with the new interactions found in the frame.

    Args:
        intermap_dict (dict): Dictionary with the interactions.
        inter_list (List): List of interactions for each frame.
        inter_type (str): Type of interaction.
        labels (ndarray): Labels of the atoms in the selection.
        frame (int): Frame identifier.

    Returns:
        intermap_dict (dict): Updated dictionary with the interactions.
    """
    for i, inter_frame in enumerate(inter_list):
        k = inter_frame[0, 2]
        for r1, r2 in inter_frame[:, :2]:
            intermap_dict[(labels[r1], labels[r2], inter_type)].append(
                frame + k)
    return intermap_dict


def intermap2df(intermap_dict, nframes, prev_cutoff=1):
    """
    Convert the intermap dictionary to a DataFrame.

    Args:
        nframes (int): Number of frames in the trajectory.
        prev_cutoff (float): Prevalence cutoff for the interactions.
        intermap_dict (dict): Dictionary with the interactions.

    Returns:
        df (DataFrame): DataFrame with the interactions.
    """

    def data_generator():
        """
        Generator to yield the data for the DataFrame creation.
        """
        for key, value in intermap_dict.items():
            prevalence = round(len(value) / nframes * 100, 2)
            if prevalence > prev_cutoff:
                yield [y for x in key for y in x.split('-')] + [prevalence] + [
                    sorted(value)]

    df = pd.DataFrame(data_generator(), columns=['resname1', 'resid1', 'name1',
                                                 'resname2', 'resid2', 'name2',
                                                 'inter', 'prevalence',
                                                 'frames'])
    return df

# %% ==========================================================================
# Debugging area
# =============================================================================
# import mdtraj as md
# from intermap.indexman import IndexManager as iman
# from intermap import topo_trajs as tt
# import numpy as np
# import time
# from collections import defaultdict
# import pandas as pd
#
# topo = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1.prmtop'
# traj = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1_sk100.nc'
#
# sel1 = "(resname =~ '(5|3)?D([ATGC])|(8OG){1}(3|5)?$')"
# sel2 = "water"
# chunk_size = 100
# n_cores = 8
#
# # %% ==========================================================================
#
# # Yield chunks of the trajectory
# chunks = tt.get_traj_chunks(topo, [traj], chunk_size=chunk_size)
# chunk = next(chunks)
# xyz = chunk.xyz
# N = len(xyz)
#
# # Get the selections
# master_traj = md.load_frame(traj, top=topo, index=0)
# seles = iman(sel1, sel2, master_traj)
# labels = seles.labels
# idxs = seles.indices
# s1_indices = seles.s1_idx
# s2_indices = seles.s2_idx
# s1_hydros = idxs['hbonds']['s1_hydros']
# s1_donors = idxs['hbonds']['s1_donors']
# s1_acc = idxs['hbonds']['s1_acc']
# s2_hydros = idxs['hbonds']['s2_hydros']
# s2_donors = idxs['hbonds']['s2_donors']
# s2_acc = idxs['hbonds']['s2_acc']
#
# # %% Run the functions in parallel ============================================
# first_timer = time.time()
# set_num_threads(n_cores)
# _ = inters(xyz[:1],
#            s1_donors, s1_hydros, s1_acc,
#            s2_donors, s2_hydros, s2_acc,
#            s1_indices, s2_indices)
#
# hb_list, cc_list = inters(xyz,
#                           s1_donors, s1_hydros, s1_acc,
#                           s2_donors, s2_hydros, s2_acc,
#                           s1_indices, s2_indices)
# print(f"Elapsed time: {time.time() - first_timer:.2f} s")
#
# # %% ==========================================================================
# # Fill intermap
# frame = 0
# intermap_dict = defaultdict(list)
# intermap_dict = update_intermap(intermap_dict, hb_list, 'hb', labels, frame)
# intermap_dict = update_intermap(intermap_dict, cc_list, 'cc', labels, frame)
# frame += N
# intermap = intermap2df(intermap_dict, N, prev_cutoff=2)
