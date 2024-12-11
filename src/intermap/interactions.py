# Created by rglez at 12/10/24

import numpy as np
from numba import njit
from numba_kdtree import KDTree as nckd

inter_identifiers = {'close_contacts': 0,
                     'anionic': 1,
                     'cationic': 2,
                     'hydrophobic': 3,
                     'metal_donor': 4}


@njit(parallel=False)
def close_contacts(xyz, k, contact_id, s1_indices, s2_indices, dist_cut):
    """
    Find the close contacts between two selections.

    Args:
        xyz (ndarray): Coordinates of the atoms in the frame.
        k (int): Frame identifier.
        contact_id: Identifier for the contact type.
        s1_indices (ndarray): Indices of the atoms in s1.
        s2_indices (ndarray): Indices of the atoms in s2.
        dist_cut: Cutoff distance for the tree query.

    Returns:
        CC_labels (list): List of labels for each close contact.
    """

    # Create & query the trees
    s2_tree = nckd(xyz[s2_indices])
    ball_1 = s2_tree.query_radius(xyz[s1_indices], dist_cut)

    # Find the close contacts
    n_contacts = sum([len(x) for x in ball_1])
    contact_indices = np.zeros((n_contacts, 4), dtype=np.int32)

    # Fill the labels
    counter = -1
    for i, x in enumerate(ball_1):
        X1 = s1_indices[i]
        for j in x:
            X2 = s2_indices[j]
            counter += 1
            contact_indices[counter][0] = X1
            contact_indices[counter][1] = X2
            contact_indices[counter][2] = k
            contact_indices[counter][3] = contact_id

    # Remove idems (self-contacts appearing if both selections overlap)
    idems = contact_indices[:, 0] == contact_indices[:, 1]
    contact_indices = contact_indices[~idems]
    return contact_indices


# %% ==========================================================================
# Instantiate the IndexManager class
# =============================================================================
# todo: implement pi-cation
# todo: keep an eye on selections for contacts calculation (directed vs undirected)
# todo: watch angular definitions for the interactions
# todo: assign pi-type based on the angle between the planes
# todo: test interactions
# todo: reorganize code
# todo: parallelize
# todo: benchmark

import time
from intermap.indices import IndexManager

start_time = time.time()
topo = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1.prmtop'
traj = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1_sk100.nc'
sel1 = "nucleic or resname 8OG"
sel2 = "protein"
inters = 'all'

iman = IndexManager(topo, traj, sel1, sel2, inters)
load_time = time.time() - start_time
print(f"Until setting indices via SMARTS: {load_time:.2f} s")

# =============================================================================
# Load coordinates
# =============================================================================
import intermap.topo_trajs as tt

# Trajectory load
start = 0
last = 1
stride = 1
chunk_size = (last - start) // 5
chunks = tt.get_traj_chunks(topo, [traj],
                            start=start, last=last, stride=stride,
                            chunk_size=chunk_size)
chunk = next(chunks)
xyz = chunk.xyz[0]
k = 0

chunk_time = time.time() - start_time
print(f"Until loading the chunk: {chunk_time:.2f} s")

# =============================================================================
# %% Start computing interactions
# =============================================================================
start = time.time()
import intermap.cutoffs as cf

# ==== close contacts =========================================================
s1_indices = iman.sel1_idx
s2_indices = iman.sel2_idx
dist_cut = cf.close_contacts
cc_id = inter_identifiers['close_contacts']
if s1_indices.size and s2_indices.size:
    cc = close_contacts(xyz, k, cc_id, s1_indices, s2_indices, dist_cut)

# ==== anionic ================================================================
s1_ani = np.intersect1d(iman.sel1_idx, iman.anions)
s2_cat = np.intersect1d(iman.sel2_idx, iman.cations)
anionic_id = inter_identifiers['anionic']
if s1_ani.size and s2_cat.size:
    ani = close_contacts(xyz, k, anionic_id, s1_ani, s2_cat, cf.ionic)

# ==== cationic ===============================================================
s1_cat = np.intersect1d(iman.sel1_idx, iman.cations)
s2_ani = np.intersect1d(iman.sel2_idx, iman.anions)
cationic_id = inter_identifiers['cationic']
if s1_cat.size and s2_ani.size:
    cat = close_contacts(xyz, k, cationic_id, s1_cat, s2_ani, cf.ionic)

# ==== hydrophobic ============================================================
s1_hp = np.intersect1d(iman.sel1_idx, iman.hydroph)
s2_hp = np.intersect1d(iman.sel2_idx, iman.hydroph)
hydrop_id = inter_identifiers['hydrophobic']
if s1_hp.size and s2_hp.size:
    hp = close_contacts(xyz, k, hydrop_id, s1_hp, s2_hp, cf.hydrophobic)

# ==== metal donor ============================================================
s1_met = np.intersect1d(iman.sel1_idx, iman.metal_don)
s2_acc = np.intersect1d(iman.sel2_idx, iman.metal_acc)
metd_id = inter_identifiers['metal_donor']
if s1_met.size and s2_acc.size:
    metd = close_contacts(xyz, k, metd_id, s1_met, s2_acc, cf.metallic)

#
print(f"Computing time: {time.time() - start:.2f} s")
# Selections
# s1_donors = np.intersect1d(self.s1_idx, self.hb_D)
# s1_donors_idx = npi.indices(self.hb_D, s1_donors)
# s1_hydros = self.hb_H[s1_donors_idx]
# s1_acc = np.intersect1d(self.s1_idx, self.hb_A)
#
# s2_donors = np.intersect1d(self.s2_idx, self.hb_D)
# s2_donors_idx = npi.indices(self.hb_D, s2_donors)
# s2_hydros = self.hb_H[s2_donors_idx]
# s2_acc = np.intersect1d(self.s2_idx, self.hb_A)
#
# s1_hydroph = np.intersect1d(self.s1_idx, self.hydroph)
# s2_hydroph = np.intersect1d(self.s2_idx, self.hydroph)
#
# s1_xdonors = np.intersect1d(self.s1_idx, self.xb_D)
# s1_xdonors_idx = npi.indices(self.xb_D, s1_xdonors)
# s1_xhydros = self.xb_H[s1_xdonors_idx]
# s1_xacc = np.intersect1d(self.s1_idx, self.xb_A)
#
# s2_xdonors = np.intersect1d(self.s2_idx, self.xb_D)
# s2_xdonors_idx = npi.indices(self.xb_D, s2_xdonors)
# s2_xhydros = self.xb_H[s2_xdonors_idx]
# s2_xacc = np.intersect1d(self.s2_idx, self.xb_A)
#
# max_vdw = self.get_max_vdw_dist()
# vdw_radii = self.radii
#
# padded_rings = self.rings
# contact_id = 3
# ctd_dist = 0.6
# min_dist = 0.38
#
# sel1_rings = padded_rings[np.isin(padded_rings[:, 0], self.s1_idx)]
# sel2_rings = padded_rings[np.isin(padded_rings[:, 0], self.s2_idx)]
#

#


# if self.s1_idx.size and self.s2_idx.size:
#     vdw = fnd.vdw_contacts(xyz, k, 1,
#                            self.s1_idx, self.s2_idx,
#                            max_vdw, vdw_radii)
#

#
# if s1_donors.size and s1_hydros.size and s1_acc.size and \
#         s2_donors.size and s2_hydros.size and s2_acc.size:
#     hb = fnd.double_contacts(xyz, k, 0,
#                              s1_donors, s1_hydros, s1_acc,
#                              s2_donors, s2_hydros, s2_acc,
#                              cf.HA_cut, cf.DA_cut, cf.DHA_cut)
#
# if self.xb_D.size and self.xb_H.size and self.xb_A.size:
#     xb = fnd.double_contacts(xyz, k, 1,
#                              self.xb_D, self.xb_H, self.xb_A,
#                              self.xb_D, self.xb_H, self.xb_A,
#                              cf.HA_cut, cf.DA_cut, cf.DHA_cut)
#
# if sel1_rings.size and sel2_rings.size:
#     pipi = fnd.pi_stacking(xyz, k, 7, sel1_rings, sel2_rings, 0.6, 0.38, 0, 90)
#
# print(f"Computing time: {time.time() - start:.2f} s")
