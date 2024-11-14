import cupy as cp
import cupyx.scipy.spatial
import mdtraj as md
from cupyx.scipy.spatial import KDTree

top = '/media/rglez/Expansion/RoyData/oxo-8/raw/A1/8oxoGA1_1_dry.prmtop'
traj = '/media/rglez/Expansion/RoyData/oxo-8/raw/A1/8oxoGA1_1_dry.nc'

parsed = md.load_netcdf(traj, top=top)
xyz = parsed.xyz

# use cupy to build kdtree

import cupy as cp
import numpy as np

# Sample XYZ data (replace with your actual loading)
xyz_data = np.random.rand(1000, 3)  # 1000 atoms, 3 coordinates each

# Transfer to CuPy
xyz_cupy = cp.asarray(xyz)
