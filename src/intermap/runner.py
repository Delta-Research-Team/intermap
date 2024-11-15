# Created by gonzalezroy at 11/14/24
import intermap.descriptors.topo_trajs as tt
import intermap.descriptors.main as mm
topo = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A1/8oxoGA1_1_hmr.prmtop'
traj = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A1/8oxoGA1_1_sk100.nc'

heavies = ['S', 'N', 'O', 'P']

start = 0
stride = 1
last = 200
chunk_size = 250

DA_cut = 0.39  # distance between Donor and Acceptor cutoff
HA_cut = 0.25  # distance between Hydrogen and Acceptor cutoff
DHA_cut = 90  # angle between Donor, Hydrogen and Acceptor cutoff

sel1 = "(resname =~ '(5|3)?D([ATGC])|(8OG){1}(3|5)?$')"
sel2 = "water"

mini_traj, resids_to_atoms, resids_to_noh, donors, hydros, acceptors = \
    tt.prepare_datastructures(topo, traj, heavies)

# Get the trajectory xyz chunks
chunks = tt.get_traj_chunks(topo, [traj])
for chunk in chunks:
    xyz = chunk.xyz


ave_min_dist, occ_nb, cp, occ_sb, occ_hb, occ_int, mi, gc = \
    mm.compute_descriptors(mini_traj, trajs, arg, resids_to_atoms,
                           resids_to_noh, calphas, oxy, nitro, donors,
                           hydros, acceptors, corr_indices, first_timer)
# invert CP matrix
