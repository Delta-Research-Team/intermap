# Created by gonzalezroy at 11/14/24
import intermap.descriptors.topo_trajs as tt
import mdtraj as md
from intermap.indexman import IndexManager as iman
from scipy.spatial import cKDTree as ckd

# =============================================================================
# User defined variables
# =============================================================================
topo = '/media/gonzalezroy/Expansion/RoyData/oxo-8/raw/water/A1/8oxoGA1_1_hmr.prmtop'
traj = '/media/gonzalezroy/Expansion/RoyData/oxo-8/raw/water/A1/8oxoGA1_1_sk100.nc'
sel1 = "(resname =~ '(5|3)?D([ATGC])|(8OG){1}(3|5)?$')"
sel2 = "water"

# Get the selections
master_traj = md.load_frame(traj, top=topo, index=0)
seles = iman(sel1, sel2, master_traj)
s1_donors, s1_hydros, s1_acc = seles.s1_donors, seles.s1_hydros, seles.s1_acc
s2_donors, s2_hydros, s2_acc = seles.s2_donors, seles.s2_hydros, seles.s2_acc

# Get the trajectory xyz chunks
chunks = tt.get_traj_chunks(topo, [traj])
chunk_parsed = next(chunks)
chunk_xyz = chunk_parsed.xyz
xyz = chunk_xyz[0]

# =============================================================================
# Find hbonds
# =============================================================================
DA_cut = 0.39  # distance between Donor and Acceptor cutoff
HA_cut = 0.25  # distance between Hydrogen and Acceptor cutoff
DHA_cut = 90  # angle between Donor, Hydrogen and Acceptor cutoff

s1_D_tree = ckd(xyz[s1_donors])
s2_A_tree = ckd(xyz[s2_acc])
s2_D_tree = ckd(xyz[s2_donors])
s1_A_tree = ckd(xyz[s1_acc])

neihgbors1 = s1_D_tree.query_ball_tree(s2_A_tree, DA_cut)
D1 = s1_donors
H1 = s1_hydros
A1 = [s2_acc[x] for x in neihgbors1]
DHA1 = []
for i, x in enumerate(D1):
    D = D1[i]
    H = H1[i]
    for A in A1[i]:
        DHA1.append((D, H, A))


neighbors2 = s2_D_tree.query_ball_tree(s1_A_tree, DA_cut)

# =============================================================================
# Testing
# =============================================================================
import numpy as np


def calc_angle(d, h, a):
    """
    Computes the angle between three atoms

    Args:
        d (donor): Coordinates of the first atom (x, y, z)
        h (hydrogen): Coordinates of the second atom (x, y, z).
        a (acceptor): Coordinates of the third atom (x, y, z).

    Returns:
        angle_deg: the angle in degrees
    """
    # Compute vectors
    dh = d - h
    ah = a - h

    # Compute dot & norms
    dot_product = np.dot(dh, ah)
    dh_norm = np.linalg.norm(dh)
    ah_norm = np.linalg.norm(ah)

    # Compute angle
    angle_rad = np.arccos(dot_product / (dh_norm * ah_norm))
    angle_deg = np.rad2deg(angle_rad)
    return angle_deg


# Convert to a universal function
vec_angle = np.vectorize(calc_angle, signature='(n),(n),(n)->()')

# Example usage with multiple sets of points
points = np.array([
    [[0, 0], [1, 0], [0, 1]],  # First triangle
    [[0, 0], [1, 1], [1, 0]]  # Second triangle
])

# Split points into separate arrays
A = points[:, 0]
B = points[:, 1]
C = points[:, 2]

# Calculate angles
angles = vec_angle(A, B, C)

# Convert to degrees
angles_degrees = np.degrees(angles)

print("Angles in degrees:", angles_degrees)
