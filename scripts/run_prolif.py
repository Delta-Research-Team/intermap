# Created by gonzalezroy at 3/4/25


import MDAnalysis as mda
import pandas as pd
import prolif as plf
from rgpack import generals as gnl

# =============================================================================
#
# =============================================================================
out_name = '/media/rglez/Expansion1/RoyData/intermap/correctness/prolif/lig-prot.pkl'

# load topology and trajectory
u = mda.Universe(plf.datafiles.TOP, plf.datafiles.TRAJ)

# create selections for the ligand and protein
ligand_selection = u.select_atoms("resname LIG")
protein_selection = u.select_atoms("protein")

# use default interactions
fp = plf.Fingerprint(
    ['Anionic', 'CationPi', 'Cationic', 'EdgeToFace', 'FaceToFace',
     'HBAcceptor', 'HBDonor', 'Hydrophobic', 'MetalAcceptor', 'MetalDonor',
     'PiCation', 'PiStacking', 'VdWContact', 'XBAcceptor', 'XBDonor'],
    count=True, parameters={"VdWContact": {"preset": "rdkit"}})
fp.run(u.trajectory, ligand_selection, protein_selection)
fp.to_pickle(out_name)
df = fp.to_dataframe()
gnl.pickle_to_file(df, out_name)

u.trajectory[0]
prot = plf.Molecule.from_mda(protein_selection)
lig = plf.Molecule.from_mda(ligand_selection)
data = fp.generate(lig, prot, metadata=True)
records = []
for residue_pair, interactions in data.items():
    ligand_residue, protein_residue = residue_pair
    ligand_residue_str = f"{ligand_residue}"
    protein_residue_str = f"{protein_residue}"
    for interaction, details in interactions.items():
        for detail in details:
            record = {
                'interaction': interaction,
                'ligand_residue': ligand_residue_str,
                'protein_residue': protein_residue_str,
                'ligand_indices': detail['indices']['ligand'],
                'protein_indices': detail['indices']['protein'],
                'ligand_parent_indices': detail['parent_indices']['ligand'],
                'protein_parent_indices': detail['parent_indices']['protein'],
                'distance': detail['distance'],
                'DHA_angle': detail.get('DHA_angle'),
                'plane_angle': detail.get('plane_angle'),
                'normal_to_centroid_angle': detail.get('normal_to_centroid_angle'),
                'intersect_distance': detail.get('intersect_distance')
            }
            records.append(record)
df = pd.DataFrame(records)
df.to_csv('metadata.csv', index=False, na_rep='NaN')

import numpy as np
liggy = np.asarray(range(4988, 5067))
proty = np.asarray(range(0, 4988))
