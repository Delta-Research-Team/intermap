# Created by rglez at 5/17/25
import mdtraj as md

# IgG3
topo = '/media/rglez/Roy2TB/Dropbox/RoyData/intermap/100_frames_trajs/02_IgG3/M1_IgG3.pdb'
traj = '/media/rglez/Roy2TB/Dropbox/RoyData/intermap/100_frames_trajs/02_IgG3/IgG3.xtc'

parsed = md.load(topo)
all = parsed.topology.select('all')
prot1 = parsed.topology.select('chainid 0 1')
prot2 = parsed.topology.select('chainid 2 3 4 5')
print(len(all))
print(len(prot1))
print(len(prot2))

import matplotlib.pyplot as plt

a = [2, 5, 6, 4, 7, 88, 8, 2, 3, 6, 9, 10, 12, 17]
plt.boxplot(a)
plt.show()
