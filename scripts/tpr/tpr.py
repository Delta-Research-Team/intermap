# Created by rglez at 5/8/25
from os.path import join

import MDAnalysis as mda

# =============================================================================
#
# =============================================================================
topo = '/home/rglez/test_data/molecular_dynamics.tpr'
traj = '/home/rglez/test_data/molecular_dynamics_sample.xtc'
out_dir = '/home/rglez/RoyHub/intermap/scripts/tpr'

# =============================================================================
#
# =============================================================================
out_top = join(out_dir, 'mpro.pdb')
u = mda.Universe(topo, traj)
ag = u.select_atoms("all")
ag.write(out_top)

# =============================================================================
#
# =============================================================================
