# Created by rglez at 12/10/24
import numpy as np
from numba import njit, prange
from numba_kdtree import KDTree as nckd
from numpy import concatenate as concat

import intermap.commons as cmn


# todo: Why loading all atoms instead of selected indices? (in nowat trajs will speed up)
# todo: remove concatenation of arrays and use slicing of preallocated arrays instead


@njit(parallel=False, cache=True)
def detect_vdw(inter_name, vdw_radii, row1, row2, dists, to_compute_others):
    """
    Detect the Van der Waals interactions

    Args:
        inter_name (str): Name of the interaction
        vdw_radii (ndarray): Van der Waals radii of the atoms
        row1 (ndarray): Indices of the atoms in the first selection
        row2 (ndarray): Indices of the atoms in the second selection
        dists (ndarray): Distances between the atoms in the first and second
                         selections
        to_compute_others (ndarray): Interactions to compute

    Returns:
        inter_idx (int): Index of the interaction in the to_compute_others array
        passing_dist (ndarray): Boolean array with the atoms that pass the
                                distance criterion
    """
    inter_idx = cmn.indices(to_compute_others, [inter_name])[0]

    if row1.size == 0:
        return inter_idx, np.zeros(dists.size, dtype=np.bool_)

    row1_vdw = vdw_radii[row1]
    row2_vdw = vdw_radii[row2]
    vdw_sum = (row1_vdw + row2_vdw)
    passing_dist = dists <= vdw_sum
    return inter_idx, passing_dist


@njit(parallel=False, cache=True)
def detect_1d(inter_name, dists, row1, type1, row2, type2, cutoffs_others,
              to_compute_others):
    """
    Detect the 1D interactions (only one distance involved)

    Args:
        inter_name (str): Name of the interaction
        dists (ndarray): Distances between the atoms in the first and second
        row1 (ndarray): Indices of the atoms in the first selection
        type1 (ndarray): Indices of the atoms of the first type
        row2 (ndarray): Indices of the atoms in the second selection
        type2 (ndarray): Indices of the atoms of the second type

    Returns:
        inter_idx (int): Index of the interaction in the to_compute_others array
        passing_dist (ndarray): Boolean array with the atoms that pass the
                                distance criterion
    """
    inter_idx = cmn.indices(to_compute_others, [inter_name])[0]

    if type1.size == 0 or type2.size == 0:
        return inter_idx, np.zeros(dists.size, dtype=np.bool_)

    dist_cutoff = cutoffs_others[0, inter_idx]
    passing_dists = dists <= dist_cutoff

    if inter_name == 'CloseContacts':
        return inter_idx, passing_dists

    else:
        s1_is_type = cmn.isin(row1, type1)
        s2_is_type = cmn.isin(row2, type2)
        are_type = s1_is_type & s2_is_type
        return inter_idx, passing_dists & are_type


@njit(parallel=False, cache=True)
def detect_2d1a(inter_name, dists, xyz, row1, row2, hb_acc, hb_hydros,
                hb_donors,
                ha_cut, da_cut, min_ang, max_ang, to_compute_others):
    """
    Detect the 2D1A interactions (two distances and one angle involved)

    Args:
        inter_name (str): Name of the interaction
        dists (ndarray): Distances between the atoms in the first and second
        row1 (ndarray): Indices of the atoms in the first selection
        row2 (ndarray): Indices of the atoms in the second selection
        hb_acc (ndarray): Indices of the hydrogen bond acceptors
        hb_hydros (ndarray): Indices of the hydrogen bond hydrogens
        hb_donors (ndarray): Indices of the hydrogen bond donors
        ha_cut (float): Cutoff distance for the hydrogen bond acceptor
        da_cut (float): Cutoff distance for the hydrogen bond donor-acceptor
        min_ang (float): Minimum angle for the hydrogen bond
        max_ang (float): Maximum angle for the hydrogen bond

    Returns:
        inter_idx (int): Index of the interaction in the to_compute_others array
        hbonds (ndarray): Boolean array with the atoms that pass the distance
    """
    idx_name = cmn.indices(to_compute_others, [inter_name])[0]

    if 'Acceptor' in inter_name:
        are_acceptors = cmn.isin(row1, hb_acc)
        are_hydros = cmn.isin(row2, hb_hydros)
        passing_ha = are_hydros & are_acceptors & (dists <= ha_cut)
        acceptors = row1[passing_ha]
        hydros = row2[passing_ha]

    elif 'Donor' in inter_name:
        are_hydros = cmn.isin(row1, hb_hydros)
        are_acceptors = cmn.isin(row2, hb_acc)
        passing_ha = are_hydros & are_acceptors & (dists <= ha_cut)
        hydros = row1[passing_ha]
        acceptors = row2[passing_ha]
    else:
        raise ValueError(f"Invalid interaction name: {inter_name}")

    if passing_ha.size == 0:
        return idx_name, np.zeros(dists.size, dtype=np.bool_)

    ext_idx = cmn.indices(hb_hydros, hydros)
    donors = hb_donors[ext_idx]

    ha_vectors = xyz[acceptors] - xyz[hydros]
    hd_vectors = xyz[donors] - xyz[hydros]
    da_vectors = xyz[donors] - xyz[acceptors]

    da_dists = np.sqrt((da_vectors ** 2).sum(axis=1))
    passing_da = da_dists <= da_cut

    if not passing_da.any():
        return idx_name, np.zeros(dists.size, dtype=np.bool_)

    angles = cmn.calc_angles_2v(hd_vectors[passing_da], ha_vectors[passing_da])
    passing_angles = (angles >= min_ang) & (angles <= max_ang)

    if passing_angles.size == 0:
        return idx_name, np.zeros(dists.size, dtype=np.bool_)

    hbonds = np.zeros(dists.size, dtype=np.bool_)
    hbonds[passing_da] = passing_angles
    return idx_name, hbonds


@njit(parallel=False, cache=True)
def stackings(inter_name, dists, mindists, s1_normals, s2_normals, cutoffs_aro,
              to_compute_aro):
    """
    Helper function to compute the pi-stacking interactions

    """

    # Parse the cutoffs
    idx = cmn.indices(to_compute_aro, [inter_name])[0]
    dist_cut = cutoffs_aro[0, idx]
    min_dist = cutoffs_aro[1, idx]
    min_ang = cutoffs_aro[2, idx]
    max_ang = cutoffs_aro[3, idx]

    # Apply restraints
    passing_dist1 = dists <= dist_cut
    passing_dist2 = mindists <= min_dist
    angles = cmn.calc_angles_2v(s1_normals, s2_normals)
    passing_angles = ((angles >= min_ang) & (angles <= max_ang))
    stacking = passing_dist1 & passing_dist2 & passing_angles
    return idx, stacking


@njit(parallel=False, cache=True)
def pications(inter_name, ijf, dists, xyz2, rings_normals, rings_idx, cat_idx,
              cutoffs_aro, to_compute_aro):
    """
    Helper function to compute the pi-cation // cation-pi interactions

    """
    # Parse the cutoffs
    idx = cmn.indices(to_compute_aro, [inter_name])[0]
    dist_cut = cutoffs_aro[0, idx]
    min_ang = cutoffs_aro[2, idx]
    max_ang = cutoffs_aro[3, idx]

    # Select the pairs
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]
    if inter_name == 'PiCation':
        s1_is_type = cmn.isin(row1, rings_idx)
        s2_is_type = cmn.isin(row2, cat_idx)
    elif inter_name == 'CationPi':
        s1_is_type = cmn.isin(row1, cat_idx)
        s2_is_type = cmn.isin(row2, rings_idx)
    else:
        raise ValueError(f"Invalid interaction name: {inter_name}")
    pairs = s1_is_type & s2_is_type

    if pairs.any():
        # Calculate angles between normals and vectors
        row1_pairs = row1[pairs]
        row2_pairs = row2[pairs]
        vector_ctr_cat = xyz2[row1_pairs] - xyz2[row2_pairs]
        normals = rings_normals[cmn.indices(rings_idx, row1_pairs)]
        angles = cmn.calc_angles_2v(normals, vector_ctr_cat)

        # Apply restraints
        passing_dist = dists[pairs] <= dist_cut
        passing_angles = ((angles >= min_ang) & (angles <= max_ang))
        return idx, pairs, passing_dist & passing_angles
    else:
        return idx, pairs, pairs


@njit(parallel=False, cache=True)
def not_aro(xyz, k, s1_indices_raw, s2_indices_raw, anions, cations,
            hydrophobes, metal_donors, metal_acceptors, vdw_radii,
            hb_hydrogens, hb_donors, hb_acceptors, xb_halogens, xb_donors,
            xb_acceptors, cutoffs_others, to_compute_others):
    """
    Compute the interactions not related to aromatic rings

    Args:
        xyz (ndarray): Coordinates of the atoms in the system
        k (int): Index of the frame in the trajectory
        s1_indices_raw (ndarray): Indices of the atoms in the first selection
        s2_indices_raw (ndarray): Indices of the atoms in the second selection
        anions (ndarray): Indices of the anions
        cations (ndarray): Indices of the cations
        hydrophobes (ndarray): Indices of the hydrophobic atoms
        metal_donors (ndarray): Indices of the metal donors
        metal_acceptors (ndarray): Indices of the metal acceptors
        vdw_radii (ndarray): Van der Waals radii of the atoms
        hb_hydrogens (ndarray): Indices of the hydrogen bond hydrogens
        hb_donors (ndarray): Indices of the hydrogen bond donors
        hb_acceptors (ndarray): Indices of the hydrogen bond acceptors
        xb_acceptors (ndarray): Indices of the halogen bond acceptors
        xb_donors (ndarray): Indices of the halogen bond donors
        xb_halogens (ndarray): Indices of the halogen atoms
        cutoffs_others (ndarray): Cutoff distances for the interactions not
        to_compute_others (ndarray): Interactions to compute

    Returns:
        ijf (ndarray): Indices of the atoms in the first and second selections
        interactions (ndarray): Container for the interactions
    """

    if to_compute_others.size == 0:
        ijf = np.zeros((0, 3), dtype=np.int32)
        inters = np.zeros((0, len(to_compute_others)), dtype=np.bool_)
        return ijf, inters

    # Create & query the trees
    dist_cut = cutoffs_others[:2].max()
    s2_tree = nckd(xyz[s2_indices_raw])
    ball_1 = s2_tree.query_radius(xyz[s1_indices_raw], dist_cut)

    # Get containers
    ijf, dists, interactions = cmn.get_containers(
        xyz, k, np.arange(xyz.shape[0]), ball_1, s1_indices_raw,
        s2_indices_raw, to_compute_others)

    # Start detecting interactions
    to_detect = list(to_compute_others)
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]

    # Contacts
    inter_name = 'VdWContact'
    if inter_name in to_detect:
        inter_idx, vdw = detect_vdw(inter_name, vdw_radii, row1, row2, dists,
                                    to_compute_others)
        interactions[:, inter_idx] = vdw

    inter_name = 'CloseContacts'
    if inter_name in to_detect:
        inter_idx, close_contacts = detect_1d(inter_name, dists, row1, row1,
                                              row2, row2, cutoffs_others,
                                              to_compute_others)
        interactions[:, inter_idx] = close_contacts

    inter_name = 'Hydrophobic'
    if inter_name in to_detect:
        inter_idx, hydrophobic = detect_1d(inter_name, dists, row1,
                                           hydrophobes, row2, hydrophobes,
                                           cutoffs_others, to_compute_others)
        interactions[:, inter_idx] = hydrophobic

    # Salt bridges
    exist_all = anions.size > 0 and cations.size > 0
    find_cationic = ('Cationic' in to_detect) and exist_all
    find_anionic = ('Anionic' in to_detect) and exist_all

    if find_cationic:
        inter_name = 'Cationic'
        inter_idx, cationic = detect_1d(inter_name, dists, row1, cations, row2,
                                        anions, cutoffs_others,
                                        to_compute_others)
        interactions[:, inter_idx] = cationic

    if find_anionic in to_detect:
        inter_name = 'Anionic'
        inter_idx, anionic = detect_1d(inter_name, dists, row1, anions, row2,
                                       cations, cutoffs_others,
                                       to_compute_others)
        interactions[:, inter_idx] = anionic

    # Metallic
    exist_all = metal_donors.size > 0 and metal_acceptors.size > 0
    find_metal_don = ('MetalDonor' in to_detect) and exist_all
    find_metal_acc = ('MetalAcceptor' in to_detect) and exist_all

    if find_metal_don in to_detect:
        inter_name = 'MetalDonor'
        inter_idx, metal_donor = detect_1d(inter_name, dists, row1,
                                           metal_donors, row2, metal_acceptors,
                                           cutoffs_others, to_compute_others)
        interactions[:, inter_idx] = metal_donor

    if find_metal_acc in to_detect:
        inter_name = 'MetalAcceptor'
        inter_idx, metal_acceptors = detect_1d(inter_name, dists, row1,
                                               metal_acceptors, row2,
                                               metal_donors, cutoffs_others,
                                               to_compute_others)
        interactions[:, inter_idx] = metal_acceptors

    # HBonds
    exists_all = hb_acceptors.size > 0 and hb_donors.size > 0
    find_hb_acc = ('HBAcceptor' in to_detect) and exists_all
    find_hb_donor = ('HBDonor' in to_detect) and exists_all
    if find_hb_acc or find_hb_donor:
        idx_both = cmn.indices(to_compute_others, ['HBAcceptor'])[0]
        da_cut, ha_cut, min_ang, max_ang = cutoffs_others[:4, idx_both]

        if find_hb_acc:
            inter_name = 'HBAcceptor'
            inter_idx, hb_a = detect_2d1a(inter_name, dists, xyz, row1, row2,
                                          hb_acceptors, hb_hydrogens,
                                          hb_donors, ha_cut, da_cut, min_ang,
                                          max_ang, to_compute_others)
            interactions[:, inter_idx] = hb_a

        if find_hb_donor:
            inter_name = 'HBDonor'
            inter_idx, hb_d = detect_2d1a(inter_name, dists, xyz, row1, row2,
                                          hb_acceptors, hb_hydrogens,
                                          hb_donors, ha_cut,
                                          da_cut, min_ang, max_ang,
                                          to_compute_others)
            interactions[:, inter_idx] = hb_d

    # XBonds
    exists_all = xb_acceptors.size > 0 and xb_donors.size > 0
    find_xb_acc = ('XBAcceptor' in to_detect) and exists_all
    find_xb_donor = ('XBDonor' in to_detect) and exists_all

    if find_xb_acc or find_xb_donor:
        idx_both = cmn.indices(to_compute_others, ['XBAcceptor'])[0]
        da_cut, ha_cut, min_ang, max_ang = cutoffs_others[:4, idx_both]

        if find_xb_acc:
            inter_name = 'XBAcceptor'
            inter_idx, xb_a = detect_2d1a(inter_name, dists, xyz, row1, row2,
                                          xb_acceptors, xb_halogens, xb_donors,
                                          ha_cut, da_cut, min_ang, max_ang,
                                          to_compute_others)
            interactions[:, inter_idx] = xb_a

        if find_xb_donor:
            inter_name = 'XBDonor'
            inter_idx, xb_d = detect_2d1a(inter_name, dists, xyz, row1, row2,
                                          xb_acceptors, xb_halogens, xb_donors,
                                          ha_cut, da_cut, min_ang, max_ang,
                                          to_compute_others)
            interactions[:, inter_idx] = xb_d

    mask = cmn.get_compress_mask(interactions)
    return ijf[mask], interactions[mask]


@njit(parallel=False, cache=True)
def aro(xyz, k, s1_indices_raw, s2_indices_raw, cations, rings, cutoffs_aro,
        to_compute_aro):
    """
    Compute the aromatic interactions

    Args:
        xyz (ndarray): Coordinates of the atoms in the system
        k (int): Index of the frame in the trajectory
        s1_indices_raw (ndarray): Indices of the atoms in the first selection
        s2_indices_raw (ndarray): Indices of the atoms in the second selection
        cations (ndarray): Indices of the cations
        rings (ndarray): Indices of the aromatic rings
        cutoffs_aro (ndarray): Cutoff distances for the aromatic interactions
        to_compute_aro (ndarray): Interactions to compute

    Returns:
        ijf (ndarray): Indices of the atoms in the first and second selections
        interactions (ndarray): Container for the interactions
    """
    # =========================================================================
    # STEP I: Find all pair of atoms/centroids within the max cutoff distance
    # =========================================================================
    if to_compute_aro.size == 0:
        ijf = np.zeros((0, 3), dtype=np.int32)
        inters = np.zeros((0, len(to_compute_aro)), dtype=np.bool_)
        return ijf, inters

    # Get cations
    s1_cat = s1_indices_raw[cmn.isin(s1_indices_raw, cations)]
    s2_cat = s2_indices_raw[cmn.isin(s2_indices_raw, cations)]

    # Get the aromatic rings
    s1_rings = rings[cmn.isin(rings[:, 0], s1_indices_raw)]
    s2_rings = rings[cmn.isin(rings[:, 0], s2_indices_raw)]

    # Compute the centroids
    s1_centr = cmn.calc_centroids(s1_rings, xyz)
    s2_centr = cmn.calc_centroids(s2_rings, xyz)

    # Compute the normal vectors
    s1_at1, s1_at3, s1_at5 = s1_rings[:, 0], s1_rings[:, 2], s1_rings[:, 4]
    s2_at1, s2_at3, s2_at5 = s2_rings[:, 0], s2_rings[:, 2], s2_rings[:, 4]
    s1_norm = cmn.calc_normal_vector(xyz[s1_at1], xyz[s1_at3], xyz[s1_at5])
    s2_norm = cmn.calc_normal_vector(xyz[s2_at1], xyz[s2_at3], xyz[s2_at5])

    # Create a new xyz array with the cations & centroids only
    xyz2 = concat((xyz[s1_cat], xyz[s2_cat], s1_centr, s2_centr), axis=0)
    ext_idx = concat((s1_cat, s2_cat, s1_rings[:, 0], s2_rings[:, 0]))

    # Internal indexing for xyz2 coordinates
    n0 = s1_cat.size + s2_cat.size
    n1 = n0 + s1_centr.shape[0]
    n2 = n1 + s2_centr.shape[0]
    s1_cat_idx = np.arange(0, s1_cat.size)
    s2_cat_idx = np.arange(s1_cat.size, n0)
    s1_rings_idx = np.arange(n0, n1)
    s2_rings_idx = np.arange(n1, n2)

    # Create & query the trees
    s1_indices = concat((s1_cat_idx, s1_rings_idx))
    s2_indices = concat((s2_cat_idx, s2_rings_idx))
    dist_cut_aro = cutoffs_aro[:2].max()

    if s2_indices.size == 0 or s1_indices.size == 0:
        ijf = np.zeros((0, 3), dtype=np.int32)
        inters = np.zeros((0, len(to_compute_aro)), dtype=np.bool_)
        return ijf, inters

    s2_tree = nckd(xyz2[s2_indices])
    ball_1 = s2_tree.query_radius_parallel(xyz2[s1_indices], dist_cut_aro)
    ijf, dists, inters = cmn.get_containers(
        xyz2, k, ext_idx, ball_1, s1_indices, s2_indices, to_compute_aro)

    # =========================================================================
    # STEP II: Compute the aromatic interactions
    # =========================================================================
    set_aro = list(to_compute_aro)

    if 'PiCation' in set_aro:
        idx, pairs, pi_cat = pications('PiCation', ijf, dists, xyz2, s1_norm,
                                       s1_rings_idx, s2_cat_idx,
                                       cutoffs_aro, to_compute_aro)
        if pairs.any():
            inters[pairs, idx] = pi_cat

    if 'CationPi' in set_aro:
        idx, pairs, cat_pi = pications('CationPi', ijf, dists, xyz2, s2_norm,
                                       s2_rings_idx, s1_cat_idx,
                                       cutoffs_aro, to_compute_aro)
        if pairs.any():
            inters[pairs, idx] = cat_pi

    # Stackings
    find_PiStacking = 'PiStacking' in set_aro
    find_EdgeToFace = 'EdgeToFace' in set_aro
    find_FaceToFace = 'FaceToFace' in set_aro
    if find_PiStacking or find_EdgeToFace or find_FaceToFace:

        # Get the ring pairs
        row1, row2 = ijf[:, 0], ijf[:, 1]
        s1_is_ctr = cmn.isin(row1, s1_rings_idx)
        s2_is_ctr = cmn.isin(row2, s2_rings_idx)
        pairs = s1_is_ctr & s2_is_ctr
        if pairs.any():
            ring_pairs = ijf[pairs]
            ring_dists = dists[pairs]
            s1_normals = s1_norm[cmn.indices(s1_rings_idx, ring_pairs[:, 0])]
            s2_normals = s2_norm[cmn.indices(s2_rings_idx, ring_pairs[:, 1])]

            # Compute the minimum distance between the rings
            num_pairs = ring_pairs.shape[0]
            mindists = np.zeros(num_pairs, dtype=np.float32)
            for i in range(num_pairs):
                s1_ring = \
                    s1_rings[cmn.indices(s1_rings_idx, [ring_pairs[i, 0]])][0]
                s1_ring_idx = s1_ring[:s1_ring[-1]]

                s2_ring = \
                    s2_rings[cmn.indices(s2_rings_idx, [ring_pairs[i, 1]])][0]
                s2_ring_idx = s2_ring[:s2_ring[-1]]
                mindists[i] = cmn.calc_min_dist(xyz[s1_ring_idx],
                                                xyz[s2_ring_idx])

            if find_PiStacking:
                idx, pi_stacking = stackings('PiStacking', ring_dists,
                                             mindists, s1_normals, s2_normals,
                                             cutoffs_aro, to_compute_aro)
                inters[pairs, idx] = pi_stacking

            if find_EdgeToFace:
                idx, etf_stacking = stackings('EdgeToFace', ring_dists,
                                              mindists,
                                              s1_normals, s2_normals,
                                              cutoffs_aro, to_compute_aro)
                inters[pairs, idx] = etf_stacking

            if find_FaceToFace:
                idx, ftf_stacking = stackings('FaceToFace', ring_dists,
                                              mindists,
                                              s1_normals, s2_normals,
                                              cutoffs_aro, to_compute_aro)
                inters[pairs, idx] = ftf_stacking

    mask = cmn.get_compress_mask(inters)
    return ijf[mask], inters[mask]


def get_estimation(xyz_all, n_samples, s1_indices_raw, s2_indices_raw, cations,
                   rings, cutoffs_aro, to_compute_aro, anions, hydrophobes,
                   metal_donors, metal_acceptors, vdw_radii, hb_hydrogens,
                   hb_donors, hb_acceptors, xb_halogens, xb_donors,
                   xb_acceptors, cutoffs_others, to_compute_others,
                   factor=1.5):
    # Preallocate the arrays
    n_frames = xyz_all.shape[0]
    samples = xyz_all[::n_frames // n_samples]
    num_detected = np.zeros(n_samples, dtype=np.int32)

    # Use parallel loop to fill
    N = samples.shape[0]
    for i in range(N):
        xyz = samples[i]

        ijf_aro, inters_aro = aro(
            xyz, i, s1_indices_raw, s2_indices_raw, cations, rings,
            cutoffs_aro, to_compute_aro)

        ijf_others, inters_others = not_aro(
            xyz, i, s1_indices_raw, s2_indices_raw, anions, cations,
            hydrophobes, metal_donors, metal_acceptors, vdw_radii,
            hb_hydrogens, hb_donors, hb_acceptors, xb_halogens, xb_donors,
            xb_acceptors, cutoffs_others, to_compute_others)
        num_detected[i] = ijf_aro.shape[0] + ijf_others.shape[0]

    # Estimate the number of contacts
    ijf_vert = np.int32(num_detected.max() * factor)
    if ijf_vert.size == 0:
        raise ValueError("No interactions detected in the sampled trajectory")

    inters_hori = inters_others.shape[1] + inters_aro.shape[1]
    ijf_template = np.empty((n_frames, ijf_vert, 3), dtype=np.int32)
    inters_template = np.empty((n_frames, ijf_vert, inters_hori),
                               dtype=np.bool_)
    ijf_template.fill(0)
    inters_template.fill(False)
    return ijf_template, inters_template


@njit(parallel=True, cache=True)
def run_parallel(xyz_all, ijf_template, inters_template, len_others, len_aro,
                 s1_indices_raw, s2_indices_raw, anions, cations,
                 hydrophobes, metal_donors, metal_acceptors, vdw_radii,
                 hb_hydrogens, hb_donors, hb_acceptors, xb_halogens, xb_donors,
                 xb_acceptors, rings, cutoffs_others, to_compute_others,
                 cutoffs_aro, to_compute_aro):
    num_frames = xyz_all.shape[0]
    limits = np.zeros(num_frames, dtype=np.int32)

    # Use parallel loop to fill
    for i in prange(num_frames):
        xyz = xyz_all[i]

        # Compute the non-aromatic interactions
        ijf_others, inters_others = not_aro(
            xyz, i, s1_indices_raw, s2_indices_raw, anions, cations,
            hydrophobes, metal_donors, metal_acceptors, vdw_radii,
            hb_hydrogens, hb_donors, hb_acceptors, xb_halogens, xb_donors,
            xb_acceptors, cutoffs_others, to_compute_others)

        # Compute the aromatic interactions
        ijf_aro, inters_aro = aro(
            xyz, i, s1_indices_raw, s2_indices_raw, cations, rings,
            cutoffs_aro, to_compute_aro)

        # Fill the templates
        num_others = ijf_others.shape[0]
        num_aro = ijf_aro.shape[0]
        limits[i] = num_others + num_aro

        ijf_template[i, :num_others] = ijf_others
        ijf_template[i, num_others:limits[i]] = ijf_aro
        inters_template[i, :num_others, :len_others] = inters_others
        inters_template[i, num_others:limits[i], len_others:] = inters_aro

    # Compress the chunks to the final size
    num_pairs = limits.sum()
    ijf_final = np.empty((num_pairs, 3), dtype=np.int32)
    inters_final = np.empty((num_pairs, len_others + len_aro), dtype=np.bool_)

    for i in range(num_frames):
        start = limits[:i].sum()
        end = limits[:i + 1].sum()
        ijf_final[start:end] = ijf_template[i, :limits[i]]
        inters_final[start:end] = inters_template[i, :limits[i]]
    return ijf_final, inters_final


# =============================================================================
# Debuging area
# =============================================================================
# import intermap.config as conf
# from argparse import Namespace
# import intermap.topo_trajs as tt
# from intermap.indices import IndexManager
# import intermap.cutoffs as cf
# from numba import set_num_threads
# from intermap import interdict as idt
#
# # Load the configuration
# config_path = '/home/rglez/RoyHub/intermap/tests/imap.cfg'
# config = conf.InterMapConfig(config_path, conf.allowed_parameters)
# args = Namespace(**config.config_args)
#
# # Parsing the interactions & cutoffs
# all_inters, all_cutoffs = cf.get_inters_cutoffs(args.cutoffs)
# if isinstance(args.interactions, str) and args.interactions == 'all':
#     to_compute = all_inters
# else:
#     to_compute = args.interactions
# selected_aro, selected_others, cutoffs_aro, cutoffs_others = \
#     cmn.get_cutoffs_and_inters(to_compute, all_inters, all_cutoffs)
#
# len_others, len_aro = len(selected_others), len(selected_aro)
# cutoffs_str = {x: args.cutoffs[x] for x in args.cutoffs if x in to_compute}
#
# # Load the necessary indices for detecting each interaction
# iman = IndexManager(args.topology, args.trajectory, args.selection_1,
#                     args.selection_2, args.interactions)
# vdw_radii = iman.radii
# hydrophobes = iman.hydroph
# anions, cations = iman.anions, iman.cations
# metal_donors, metal_acceptors = iman.metal_don, iman.metal_acc
# hb_hydrogens, hb_donors, hb_acceptors = iman.hb_H, iman.hb_D, iman.hb_A
# xb_halogens, xb_donors, xb_acceptors = iman.xb_H, iman.xb_D, iman.xb_A
# rings = iman.rings
#
# # Load & trim the trajectory
# n_frames = iman.n_frames
# last = tt.parse_last_param(args.last, n_frames)
# traj_frames = np.arange(args.start, last, args.stride)
#
# # Naming all atoms and interactions
# universe = iman.universe
# sel_idx = iman.sel_idx
# atnames = universe.atoms.names[sel_idx]
# resnames = universe.atoms.resnames[sel_idx]
# resids = universe.atoms.resids[sel_idx]
# n_sel_atoms = sel_idx.size
# names = [f"{resnames[i]}_{resids[i]}_{atnames[i]}" for i in
#          range(n_sel_atoms)]
# inters = np.asarray(selected_others.tolist() + selected_aro.tolist())
#
# # TESTING
# k = 0
# set_num_threads(12)
# s1_indices = iman.sel1_idx
# s2_indices = iman.sel2_idx
# self = idt.InterDict(
#     args.format, args.min_prevalence, traj_frames, names, inters)
#
# max_allocated = 0
# chunks = tt.split_in_chunks(traj_frames, args.chunk_size)
# for i, frames_chunk in enumerate(chunks):
#     xyz_chunk = tt.get_coordinates(universe, frames_chunk, sel_idx,
#                                    n_sel_atoms)
#
#     # Estimating the number of interactions per chunk to allocate memory
#     if i == 0:
#         print('Estimating the number of interactions per chunk')
#         ijf_template, inters_template = get_estimation(
#             xyz_chunk, 5, s1_indices, s2_indices, cations, rings,
#             cutoffs_aro, selected_aro, anions, hydrophobes, metal_donors,
#             metal_acceptors, vdw_radii, hb_hydrogens, hb_donors,
#             hb_acceptors, xb_halogens, xb_donors, xb_acceptors,
#             cutoffs_others, selected_others)
#
#     # Parallel computing of the interactions
#     ijf_chunk, inters_chunk = run_parallel(
#         xyz_chunk, ijf_template, inters_template, len_others, len_aro,
#         s1_indices, s2_indices, anions, cations, hydrophobes,
#         metal_donors, metal_acceptors, vdw_radii, hb_hydrogens, hb_donors,
#         hb_acceptors, xb_halogens, xb_donors, xb_acceptors, rings,
#         cutoffs_others, selected_others, cutoffs_aro, selected_aro)
#
#     # Raise if not enough space has been allocated
#     if (occupancy := ijf_chunk.shape[0] / max_allocated) >= 0.98:
#         raise ValueError(f"Chunk {i} occupancy: {round(occupancy, 2)}")
#     elif occupancy >= 0.90:
#         print('WARNING: Allocations almost full')
#
#     # Filling the interaction dictionary
#     ijf_chunk[:, 2] = frames_chunk[ijf_chunk[:, 2]]
#     self.fill(ijf_chunk, inters_chunk, traj_frames)
