# Created by gonzalezroy at 2/6/25
"""
Test the geometric functions to compute distances and angles
"""
import mdtraj as md
import numpy as np

import intermap.commons as cmn
from tests.conftest import mdtrajectory


def test_distance(mdtrajectory):
    """
    Test the distance function
    """

    trj = mdtrajectory
    n_atoms = trj.n_atoms

    # get random indices for the tests
    idx1 = np.random.randint(0, n_atoms, 500)
    idx2 = np.random.randint(0, n_atoms, 500)

    # compute distances with mdtraj
    dists_mdtraj = md.compute_distances(trj, zip(idx1, idx2), periodic=False)

    # compute distances with imap
    xyz1 = trj.xyz[:, idx1][0]
    xyz2 = trj.xyz[:, idx2][0]
    dists_imap = cmn.calc_dist(xyz1, xyz2)

    assert np.allclose(dists_mdtraj, dists_imap, atol=1e-5)


def test_angles_3p(mdtrajectory):
    trj = mdtrajectory
    n_atoms = trj.n_atoms

    idx1 = np.random.randint(0, n_atoms, 500)
    idx2 = np.random.randint(0, n_atoms, 500)
    idx3 = np.random.randint(0, n_atoms, 500)

    abc = np.stack([idx1, idx2, idx3], axis=1)
    no_equals = (abc[:, 0] != abc[:, 1]) & (abc[:, 1] != abc[:, 2]) & (
            abc[:, 0] != abc[:, 2])
    trio = abc[no_equals]

    # compute angles with mdtraj
    angles_mdtraj_radians = md.compute_angles(trj, trio, periodic=False)
    angles_mdtraj = np.degrees(angles_mdtraj_radians)

    # compute angles with imap
    xyz1 = trj.xyz[:, trio[:, 0]][0]
    xyz2 = trj.xyz[:, trio[:, 1]][0]
    xyz3 = trj.xyz[:, trio[:, 2]][0]
    angles_imap = cmn.calc_angle(xyz1, xyz2, xyz3)

    assert np.allclose(angles_mdtraj, angles_imap, atol=1e-5)


def test_angles_2v(mdtrajectory):
    trj = mdtrajectory
    n_atoms = trj.n_atoms

    idx1 = np.random.randint(0, n_atoms, 500)
    idx2 = np.random.randint(0, n_atoms, 500)
    idx3 = np.random.randint(0, n_atoms, 500)

    abc = np.stack([idx1, idx2, idx3], axis=1)
    no_equals = (abc[:, 0] != abc[:, 1]) & (abc[:, 1] != abc[:, 2]) & (
            abc[:, 0] != abc[:, 2])
    trio = abc[no_equals]

    # compute angles with mdtraj
    angles_mdtraj_radians = md.compute_angles(trj, trio, periodic=False)
    angles_mdtraj = np.degrees(angles_mdtraj_radians)[0]

    # compute angles with imap 2v
    xyz1 = trj.xyz[:, trio[:, 0]][0]
    xyz2 = trj.xyz[:, trio[:, 1]][0]
    xyz3 = trj.xyz[:, trio[:, 2]][0]

    v1 = xyz1 - xyz2
    v2 = xyz3 - xyz2
    angles_imap_2v = cmn.calc_angles_2v(v1, v2)
    angles_imap_3p = cmn.calc_angle(xyz1, xyz2, xyz3)

    assert np.allclose(angles_mdtraj, angles_imap_2v, atol=1e-3)
    assert np.allclose(angles_mdtraj, angles_imap_3p, atol=1e-3)
    assert np.allclose(angles_imap_2v, angles_imap_3p, atol=1e-3)
