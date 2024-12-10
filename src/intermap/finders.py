# Created by rglez at 11/16/24
"""
Functions to find the interactions between two selections in a frame.
"""

import numpy as np
import pandas as pd
from numba import njit
from numba_kdtree import KDTree as nckd

from intermap.geometric import calc_angle, calc_dist


# todo: Change param names & function names to be more descriptive

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
    DHA_labels = np.empty((n_triples, 4), dtype=np.int32)
    DHA_coords = np.empty((n_triples, 3, 3), dtype=np.float32)

    counter = -1
    for i in range(len(ball_1)):
        D = s1_donors[i]
        H = s1_hydros[i]
        for j in ball_1[i]:
            A = s2_acc[j]
            counter += 1
            DHA_labels[counter] = [D, A, k, 0]
            DHA_coords[counter][0] = xyz[D]
            DHA_coords[counter][1] = xyz[H]
            DHA_coords[counter][2] = xyz[A]

    for i in range(len(ball_2)):
        D = s2_donors[i]
        H = s2_hydros[i]
        for j in ball_2[i]:
            A = s1_acc[j]
            counter += 1
            DHA_labels[counter] = [D, A, k, 0]
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

    return correct_DHA.astype(np.int32)


@njit(parallel=False)
def vdw_contacts(xyz, k, contact_id, s1_indices, s2_indices, cutoff,
                 vdw_radii):
    """
    Args:
        xyz (ndarray): Coordinates of the atoms in the frame.
        k (int): Frame identifier.
        contact_id: Identifier for the contact type.
        s1_indices (ndarray): Indices of the atoms in s1.
        s2_indices (ndarray): Indices of the atoms in s2.
        cutoff: Cutoff distance for the tree query.
        vdw_radii: Van der Waals radii of the atoms.

    Returns:
        CC_labels (list): List of labels for each vdw contact.
    """
    # Create & query the trees
    s2_tree = nckd(xyz[s2_indices])
    ball_1 = s2_tree.query_radius(xyz[s1_indices], cutoff)

    # Find the CC
    n_doubles = sum([len(x) for x in ball_1])
    CC_labels = np.zeros((n_doubles, 4), dtype=np.int32)

    counter = -1
    for i, x in enumerate(ball_1):
        X1 = s1_indices[i]
        for j in x:
            X2 = s2_indices[j]
            counter += 1
            CC_labels[counter][0] = X1
            CC_labels[counter][1] = X2
            CC_labels[counter][2] = k
            CC_labels[counter][3] = contact_id

    # Discard identical pairs
    idems = CC_labels[:, 0] == CC_labels[:, 1]
    CC_labels = CC_labels[~idems]

    # Discard pairs whose distance is greater than the sum of the VdW radii
    vdw_radii_x1 = vdw_radii[CC_labels[:, 0]]
    vdw_radii_x2 = vdw_radii[CC_labels[:, 1]]
    vdw_sum = (vdw_radii_x1 + vdw_radii_x2) / 10
    dists = calc_dist(xyz[CC_labels[:, 0]], xyz[CC_labels[:, 1]])
    vdw_contact = dists <= vdw_sum
    CC_labels = CC_labels[vdw_contact]

    return CC_labels


@njit(parallel=False)
def single_contacts(xyz, k, contact_id, s1_indices, s2_indices, cut=0.3):
    """
    Args:
        xyz (ndarray): Coordinates of the atoms in the frame.
        k (int): Frame identifier.
        contact_id: Identifier for the contact type.
        s1_indices (ndarray): Indices of the atoms in s1.
        s2_indices (ndarray): Indices of the atoms in s2.
        cut: Cutoff distance for the tree query.

    Returns:
        CC_labels (list): List of labels for each close contact.
    """

    # Create & query the trees
    s2_tree = nckd(xyz[s2_indices])
    ball_1 = s2_tree.query_radius(xyz[s1_indices], cut)

    # Find the CC
    n_doubles = sum([len(x) for x in ball_1])
    CC_labels = np.zeros((n_doubles, 4), dtype=np.int32)

    counter = -1
    for i, x in enumerate(ball_1):
        X1 = s1_indices[i]
        for j in x:
            X2 = s2_indices[j]
            counter += 1
            CC_labels[counter][0] = X1
            CC_labels[counter][1] = X2
            CC_labels[counter][2] = k
            CC_labels[counter][3] = contact_id

    idems = CC_labels[:, 0] == CC_labels[:, 1]
    CC_labels = CC_labels[~idems]
    return CC_labels


@njit(parallel=False)
def double_contacts(xyz, k, contact_id,
                    s1_donors, s1_hydros, s1_acc,
                    s2_donors, s2_hydros, s2_acc,
                    cut_HA, cut_DA, cut_DHA):
    """
    Args:
        xyz: Coordinates of the atoms in the frame.
        k (int): Frame identifier.
        contact_id: Identifier for the contact type.
        s1_donors: Indices of the donor atoms in s1.
        s1_hydros: Indices of the hydrogen atoms in s1.
        s1_acc: Indices of the acceptor atoms in s1.
        s2_donors: Indices of the donor atoms in s2.
        s2_hydros: Indices of the hydrogen atoms in s2.
        s2_acc: Indices of the acceptor atoms in s2.
        cut_HA: Cutoff distance between Hydrogen and Acceptor.
        cut_DA: Cuttoff distance between Donor and Acceptor.
        cut_DHA: Cutoff angle between Donor, Hydrogen and Acceptor.

    Returns:
        DHA_labels (list): List of DHA labels for each hydrogen bond.
    """

    # Create & query the trees
    s2_A_tree = nckd(xyz[s2_acc])
    s1_A_tree = nckd(xyz[s1_acc])

    ball_1 = s2_A_tree.query_radius(xyz[s1_hydros], cut_HA)
    ball_2 = s1_A_tree.query_radius(xyz[s2_hydros], cut_HA)

    # Find the DHA candidates
    n_triples = 0
    for x in ball_1:
        n_triples += x.size
    for x in ball_2:
        n_triples += x.size
    DHA_labels = np.empty((n_triples, 4), dtype=np.int32)
    DHA_coords = np.empty((n_triples, 3, 3), dtype=np.float32)

    counter = -1
    for i in range(len(ball_1)):
        D = s1_donors[i]
        H = s1_hydros[i]
        for j in ball_1[i]:
            A = s2_acc[j]
            counter += 1
            DHA_labels[counter] = [D, A, k, 0]
            DHA_coords[counter][0] = xyz[D]
            DHA_coords[counter][1] = xyz[H]
            DHA_coords[counter][2] = xyz[A]

    for i in range(len(ball_2)):
        D = s2_donors[i]
        H = s2_hydros[i]
        for j in ball_2[i]:
            A = s1_acc[j]
            counter += 1
            DHA_labels[counter] = [D, A, k, 0]
            DHA_coords[counter][0] = xyz[D]
            DHA_coords[counter][1] = xyz[H]
            DHA_coords[counter][2] = xyz[A]

    # Filter the DHA candidates
    DA_dist_correct = calc_dist(DHA_coords[:, 0],
                                DHA_coords[:, 2]) <= cut_DA
    DHA_angle_correct = calc_angle(DHA_coords[:, 0, :],
                                   DHA_coords[:, 1, :],
                                   DHA_coords[:, 2, :]) > cut_DHA

    correct_DHA = DHA_labels[DA_dist_correct & DHA_angle_correct]

    return correct_DHA.astype(np.int32)


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
    CC_labels = np.zeros((n_doubles, 4), dtype=np.int32)

    counter = -1
    for i, x in enumerate(ball_1):
        X1 = s1_indices[i]
        for j in x:
            X2 = s2_indices[j]
            counter += 1
            CC_labels[counter][0] = X1
            CC_labels[counter][1] = X2
            CC_labels[counter][2] = k
            CC_labels[counter][3] = 1

    idems = CC_labels[:, 0] == CC_labels[:, 1]
    CC_labels = CC_labels[~idems]
    return CC_labels


@njit(parallel=False)
def inters(xyz, frame, s1_donors, s1_hydros, s1_acc, s2_donors,
           s2_hydros, s2_acc,
           s1_indices, s2_indices):
    """
    Find the interactions between two selections in a frame.

    Args:
        frame (int): Frame identifier.
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
    list_inters = []

    for i in range(num_frames):
        hbs = hbonds(xyz[i], i + frame,
                     s1_donors, s1_hydros, s1_acc, s2_donors,
                     s2_hydros, s2_acc)
        ccs = close_contacts(xyz[i], i + frame, s1_indices, s2_indices)
        inters_i = np.concatenate((hbs, ccs))
        list_inters.append(inters_i)

    return list_inters


def get_intermap(inter_list, labels, inters_types, n, prev_cutoff=1):
    """
    Get the intermap DataFrame from the interactions list and the labels.

    Args:
        inter_list (ndarray): Array of interactions for each frame.
        labels (ndarray): Labels of the atoms in the selection.
        inters_types (ndarray): Types of interactions.
        n (int): Number of frames in the trajectory.
        prev_cutoff (float): Prevalence cutoff for the interactions.

    Returns:
        df (DataFrame): DataFrame with the interactions.
    """
    df = pd.DataFrame(inter_list, columns=['r1', 'r2', 'frame', 'type'])
    df = df.groupby(['r1', 'r2', 'type'])['frame'].unique()
    df = df.reset_index()
    df['prevalence'] = df['frame'].map(lambda x: round(len(x) / n * 100, 2))
    df.reset_index(inplace=True)

    df = df[df['prevalence'] > prev_cutoff]
    df['r1'] = df['r1'].map(lambda x: labels[x])
    df['r2'] = df['r2'].map(lambda x: labels[x])
    df['type'] = df['type'].map(lambda x: inters_types[x])
    df[['resname1', 'resid1', 'name1']] = df['r1'].str.split('-', expand=True)
    df[['resname2', 'resid2', 'name2']] = df['r2'].str.split('-', expand=True)
    df = df[
        ['resname1', 'resid1', 'name1', 'resname2', 'resid2', 'name2', 'type',
         'prevalence', 'frame']]
    df.sort_values(by='prevalence', ascending=True, inplace=True)

    return df


@njit(parallel=False)
def compute_centroids(rings, xyz):
    centroids = np.zeros((len(rings), 3), dtype=float)
    for i, ring in enumerate(rings):
        atoms = ring[:ring[-1]]
        centroids[i] = xyz[atoms].sum(axis=0) / len(atoms)
    return centroids


@njit(parallel=False)
def calc_min_dist(coords1, coords2):
    """
    Get the minimumm distance between two sets of coordinates

    Args:
        coords1: coordinates of the first residue
        coords2: coordinates of the second residue

    Returns:
        The minimum distance between two sets of coordinates
    """
    # Constants
    n1 = coords1.shape[0]
    n2 = coords2.shape[0]

    # Find minimum distance using square values to save time
    min_dist_squared = np.inf
    for i in range(n1):
        for j in range(n2):
            dist_squared = (
                    (coords1[i][0] - coords2[j][0]) ** 2
                    + (coords1[i][1] - coords2[j][1]) ** 2
                    + (coords1[i][2] - coords2[j][2]) ** 2
            )
            if dist_squared < min_dist_squared:
                min_dist_squared = dist_squared
    return np.sqrt(min_dist_squared)


@njit(parallel=False)
def calc_normal_vector(p1, p2, p3):
    """
    Calculate the normal vector of a plane defined by three points

    Args:
        p1 (np.ndarray): First point
        p2 (np.ndarray): Second point
        p3 (np.ndarray): Third point

    Returns:
        np.ndarray: Normal vector of the plane defined by the three points
    """
    # Calculate vectors from points
    v1 = p2 - p1
    v2 = p3 - p1

    # Calculate the normal vector
    normal = np.cross(v1, v2)
    norm = np.linalg.norm(normal)

    if norm == 0:
        return np.zeros(3, dtype=np.float32)

    return (normal / norm).astype(np.float32)


@njit(parallel=False)
def pi_stacking(xyz, k, contact_id,
                sel1_rings, sel2_rings,
                cut_ctd_dist, cut_min_dist, cut_min_angle, cut_max_angle):
    """
    Detect pi-stacking interactions between aromatic rings

    Returns:
        rings_finals (np.ndarray): Array with the indices of the pi-stacking interactions

    """
    sel1_centroids = compute_centroids(sel1_rings, xyz)
    sel2_centroids = compute_centroids(sel2_rings, xyz)
    s2_tree = nckd(sel2_centroids)
    ball_1 = s2_tree.query_radius(sel1_centroids, cut_ctd_dist)

    # Find the DHA candidates
    rings_list = []

    for i, ball in enumerate(ball_1):
        if ball.size == 0:
            continue
        # Detect close rings
        ring1 = sel1_rings[i][:sel1_rings[i][-1]]
        r1_1 = xyz[ring1[0]]
        r1_3 = xyz[ring1[2]]
        r1_5 = xyz[ring1[4]]

        for j in ball:
            ring2 = sel2_rings[j][:sel2_rings[j][-1]]

            # Ignore adjacent ones
            if np.intersect1d(ring1, ring2).size:
                continue

            min_ring_dist = calc_min_dist(xyz[ring1], xyz[ring2])
            if min_ring_dist <= cut_min_dist:
                # Detect angle-compliant rings
                r2_1 = xyz[ring2[0]]
                r2_3 = xyz[ring2[2]]
                r2_5 = xyz[ring2[4]]
                normal_ring1 = calc_normal_vector(r1_1, r1_3, r1_5)
                normal_ring2 = calc_normal_vector(r2_1, r2_3, r2_5)
                dot_product = np.dot(normal_ring1, normal_ring2)
                dot_64 = np.float64(dot_product)

                # Manual clipping
                if dot_64 < -1.0:
                    dot_64 = -1.0
                elif dot_64 > 1.0:
                    dot_64 = 1.0
                angle = np.arccos(dot_64)

                angle_degrees = np.degrees(angle)
                if cut_min_angle < angle_degrees < cut_max_angle:
                    rings_list.append((ring1[0], ring2[0]))

    rings_finals = np.empty((len(rings_list), 4), dtype=np.int32)
    for i, (r1, r2) in enumerate(rings_list):
        rings_finals[i] = r1, r2, k, contact_id
    return rings_finals

# %% ==========================================================================
# Debugging area
# =============================================================================
# import MDAnalysis as mda
#
# topo = '/home/rglez/RoyHub/oxo-8/data/raw/A1/8oxoGA1_1_dry.prmtop'
# traj = '/home/rglez/RoyHub/oxo-8/data/raw/A1/8oxoGA1_1_dry.nc'
# u = mda.Universe(topo, traj)
#
# aroms = u.guess_TopologyAttrs(to_guess=['aromaticities'])
#
#
# u.guess_TopologyAttrs('aromaticities')
#
#
# a = mda.guesser.default_guesser.DefaultGuesser('aromaticities')
# u.select_atoms("charge < 0")

# import mdtraj as md
# from intermap.indexman import IndexManager as iman
# from intermap import topo_trajs as tt
# import numpy as np
# import time
# import pandas as pd
# from numba import set_num_threads
#
#
# sel1 = "resname =~ '8OG'"
# sel2 = "protein"
#
# chunk_size = 2000
# n_cores = 8

# %% ==========================================================================

# Yield chunks of the trajectory
# chunks = tt.get_traj_chunks(topo, [traj], chunk_size=chunk_size)
# chunk = next(chunks)
# xyz = chunk.xyz
# N = len(xyz)
# prev_cutoff = 2
# # Get the selections
# master_traj = md.load_frame(traj, top=topo, index=0)
# seles = iman(sel1, sel2, master_traj)
# labels = seles.labels
# inters_types = np.asarray(['hb', 'cc'])
# idxs = seles.indices
# s1_indices = seles.s1_idx
# s2_indices = seles.s2_idx
# s1_hydros = idxs['hbonds']['s1_hydros']
# s1_donors = idxs['hbonds']['s1_donors']
# s1_acc = idxs['hbonds']['s1_acc']
# s2_hydros = idxs['hbonds']['s2_hydros']
# s2_donors = idxs['hbonds']['s2_donors']
# s2_acc = idxs['hbonds']['s2_acc']

# # %% Run the functions in parallel ============================================
# first_timer = time.time()
# set_num_threads(n_cores)
# cc1 = close_contacts(xyz[0], 0, s1_indices, s2_indices)
# hb1 = hbonds(xyz[0], 0, s1_donors, s1_hydros, s1_acc, s2_donors, s2_hydros,
#              s2_acc)
# estimate = int((len(cc1) + len(hb1)) * len(xyz) * 1.5)
# frame = 0
# inter_list = inters(xyz, frame,
#                     s1_donors, s1_hydros, s1_acc,
#                     s2_donors, s2_hydros, s2_acc,
#                     s1_indices, s2_indices)
# inter_list = np.concatenate(inter_list)
#
# print(f"Elapsed time: {time.time() - first_timer:.2f} s")
#
# df_map = get_intermap(inter_list, labels, inters_types, N)

# %% ==========================================================================
