# Created by rglez at 3/15/25


mdscan_trajs = (
    'arc', 'dcd', 'binpos', 'xtc', 'trr', 'hdf5', 'h5', 'ncdf', 'netcdf', 'nc',
    'pdb.gz', 'pdb', 'lh5', 'crd', 'mdcrd', 'inpcrd', 'restrt', 'rst7',
    'ncrst', 'lammpstrj', 'dtr', 'stk', 'gro', 'xyz.gz', 'xyz', 'tng', 'xml',
    'mol2', 'hoomdxml', 'gsd')

mdscan_topos = (
    'pdb', 'pdb.gz', 'h5', 'lh5', 'prmtop', 'parm7', 'prm7', 'psf', 'mol2',
    'hoomdxml', 'gro', 'arc', 'hdf5', 'gsd')

mdanalysis_topos = (
    'psf', 'pdb', 'ent', 'pqr', 'pdbqt', 'gro', 'top', 'prmtop', 'parm7',
    'dms', 'tpr', 'itp', 'mol2', 'data', 'lammpsdump', 'xyz', 'txyz', 'arc',
    'gms', 'log', 'config', 'history', 'xml', 'gsd', 'mmtf', 'in')

mdanalysis_trajs = (
    'dcd', 'data', 'lammpsdump', 'xtc', 'trr', 'xyz', 'txyz', 'arc', 'gsd',
    'gms', 'log', 'out', 'trj', 'mdcrd', 'inpcrd', 'restrt', 'ncdf', 'nc',
    'pdb', 'ent', 'pdbqt', 'pqr', 'gro', 'crd', 'dms', 'trz', 'mol2', 'config',
    'history', 'mmtf', 'coor', 'namdbin', 'in')

imap_topos = set.intersection(set(mdscan_topos), set(mdanalysis_topos))
imap_trajs = set.intersection(set(mdscan_trajs), set(mdanalysis_trajs))
