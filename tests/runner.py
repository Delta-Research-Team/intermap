# Created by rglez at 1/19/25
from os.path import join

import MDAnalysis as mda
import prolif as plf

import intermap.config as conf


def run_prolif(out_dir):
    """
    Run prolif on a test trajectory and save the results to a pickle file.

    Args:
        out_dir (str): The directory to save the output files
    """
    topo = join(conf.proj_dir, 'tests', 'data', 'traj1', 'top.pdb')
    traj = join(conf.proj_dir, 'tests', 'data', 'traj1', 'traj.xtc')
    cases = {'lig_prot': {'sel1': 'resname LIG', 'sel2': 'protein'},
             'prot_prot': {'sel1': 'protein', 'sel2': 'protein'}}

    for case in cases:
        u = mda.Universe(topo, traj)
        s1 = u.select_atoms(cases[case]['sel1'])
        s2 = u.select_atoms(cases[case]['sel2'])

        fp = plf.Fingerprint(['Anionic',
                              'CationPi',
                              'Cationic',
                              'EdgeToFace',
                              'FaceToFace',
                              'HBAcceptor',
                              'HBDonor',
                              'Hydrophobic',
                              'MetalAcceptor',
                              'MetalDonor',
                              'PiCation',
                              'PiStacking',
                              'VdWContact',
                              'XBAcceptor',
                              'XBDonor'], count=True)
        fp.run(u.trajectory[::10], s1, s2)
        out_name = join(out_dir, f'{case}_ifp.pkl')
        fp.to_pickle(out_name)


# =============================================================================
#
# =============================================================================
out_dir = join(conf.proj_dir, 'tests', 'outputs')
run_prolif(out_dir)
