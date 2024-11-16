# Created by rglez at 11/16/24
import numpy as np
from scipy.spatial import cKDTree as ckd

from intermap.geometric import calc_angle, calc_dist


def find_hbonds(xyz, s1_donors, s1_hydros, s1_acc,
                s2_donors, s2_hydros, s2_acc):
    """
    Args:
        xyz: Coordinates of the atoms in the frame.
        s1_donors: Indices of the donor atoms in s1.
        s1_hydros: Indices of the hydrogen atoms in s1.
        s1_acc: Indices of the acceptor atoms in s1.
        s2_donors: Indices of the donor atoms in s2.
        s2_hydros: Indices of the hydrogen atoms in s2.
        s2_acc: Indices of the acceptor atoms in s2.

    Returns:
        neihgbors1 (list): List of neighbors for each hydrogen in s1.
        ball_2 (list): List of neighbors for each hydrogen in s2.
    """
    # Hardcoded cutoffs
    HA_cut = 0.25  # distance between Hydrogen and Acceptor cutoff
    DA_cut = 0.39  # distance between Donor and Acceptor cutoff
    DHA_cut = 90  # angle between Donor, Hydrogen and Acceptor cutoff

    # Create & query the trees
    s1_H_tree = ckd(xyz[s1_hydros])
    s2_A_tree = ckd(xyz[s2_acc])
    s2_H_tree = ckd(xyz[s2_hydros])
    s1_A_tree = ckd(xyz[s1_acc])

    ball_1 = s1_H_tree.query_ball_tree(s2_A_tree, HA_cut)
    ball_2 = s2_H_tree.query_ball_tree(s1_A_tree, HA_cut)
    del s1_H_tree, s2_A_tree, s2_H_tree, s1_A_tree

    # Find the DHA candidates
    n_triples = sum([len(x) for x in ball_1] + [len(x) for x in ball_2])
    DHA_labels = np.zeros((n_triples, 3), dtype=np.int64)
    DHA_coords = np.zeros((n_triples, 3, 3), dtype=float)

    counter = -1
    for i, x in enumerate(ball_1):
        D = s1_donors[i]
        H = s1_hydros[i]
        for j in x:
            A = s2_acc[j]
            counter += 1
            DHA_labels[counter] = [D, H, A]
            DHA_coords[counter] = [xyz[D], xyz[H], xyz[A]]

    for i, x in enumerate(ball_2):
        D = s2_donors[i]
        H = s2_hydros[i]
        for j in x:
            A = s1_acc[j]
            counter += 1
            DHA_labels[counter] = [D, H, A]
            DHA_coords[counter] = [xyz[D], xyz[H], xyz[A]]

    # Filter the DHA candidates
    DA_dist_correct = calc_dist(DHA_coords[:, 0],
                                DHA_coords[:, 2]) <= DA_cut
    DHA_angle_correct = calc_angle(DHA_coords[:, 0, :],
                                   DHA_coords[:, 1, :],
                                   DHA_coords[:, 2, :]) > DHA_cut
    return DHA_labels[DA_dist_correct & DHA_angle_correct]
