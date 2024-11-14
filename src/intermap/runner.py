# Created by gonzalezroy at 11/14/24
import intermap.topo_trajs as tt

top = '/home/gonzalezroy/Dropbox/RoyHub/oxo-8/data/raw/water/A1/8oxoGA1_1_hmr.prmtop'
traj = '/home/gonzalezroy/Dropbox/RoyHub/oxo-8/data/raw/water/A1/8oxoGA1_1_sk100.nc'


sel1 = "(resname =~ '(5|3)?D([ATGC])|(8OG){1}(3|5)?$')"
sel2 = "water"
topo_df


# Get the trajectory xyz chunks
chunks = tt.get_traj_chunks(top, [traj])
for chunk in chunks:
    xyz = chunk.xyz



DA_cut = 0.39  # distance between Donor and Acceptor cutoff
HA_cut = 0.25  # distance between Hydrogen and Acceptor cutoff
DHA_cut = 90  # angle between Donor, Hydrogen and Acceptor cutoff
heavy = ['S', 'N', 'O', 'P']

start = 0
stride = 1
last = 200

chunk_size = 250
