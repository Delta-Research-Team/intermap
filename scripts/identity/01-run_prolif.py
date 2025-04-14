# Created by gonzalezroy at 4/14/25
import os
from os.path import join

import MDAnalysis as mda
import numpy as np
import prolif as plf
import rdkit
import rgpack.generals as gnl


def guess_from_name(name, mass, pt_symbols, pt_masses, real_names):
    """
    Guess the element from the name of the atom
    """

    if len(name) > 1:
        name1 = name.strip().lower()[0]
        name2 = name.strip().lower()[0:2]
    else:
        name1 = name.strip().lower()[0]
        if real_name := real_names[name1]:
            return real_name
        else:
            raise ValueError(f"Unknown element: {name}")

    if name2 in real_names:
        mass2 = round(pt_masses[pt_symbols[real_names[name2]]])
        if mass2 == mass:
            return real_names[name2]
    return real_names[name1]


def get_periodic_table_info():
    """
    Get the periodic table information
    """
    pt = rdkit.Chem.GetPeriodicTable()
    pt_symbols = {pt.GetElementSymbol(x): x for x in range(1, 119)}
    pt_masses = {x: pt.GetAtomicWeight(x) for x in range(1, 119)}
    real_names = {x.lower(): x for x in pt_symbols}
    return pt_symbols, pt_masses, real_names


def guess_elements(universe):
    masses = [round(x) for x in universe.atoms.masses]
    names = universe.atoms.names

    # Ensure all elements are present
    try:
        elements = universe.atoms.elements
    except Exception:
        pt_symbols, pt_masses, real_names = get_periodic_table_info()
        elements = []
        for i, name in enumerate(names):
            element = guess_from_name(
                name, masses[i], pt_symbols, pt_masses, real_names)
            if not element:
                raise ValueError(
                    f"Fail to guess element from {name}."
                    f" Please check the topology.")
            elements.append(element)
        elements = np.asarray(elements)
    return elements


# =============================================================================
# TopoTrajs
# =============================================================================
user = os.getlogin()

topo_trajs = {
    'MPRO':
        {
            'topo': f'/media/{user}/Expansion/romie/TRAJECTORIES_INPUTS_DATA_mpro_wt_variants_amarolab/a173v/a173v_mpro_chainA_rep123.pr5.aligned_CA.not_waters_or_ions.psf',
            'traj': f'/media/{user}/Expansion/romie/TRAJECTORIES_INPUTS_DATA_mpro_wt_variants_amarolab/a173v/a173v_mpro_chainA_rep1.pr5.aligned_CA.not_waters.dcd'}}

# =============================================================================
#
# =============================================================================
out_dir = f'/home/{user}/RoyHub/intermap/scripts/identity/runs/mpro/prolif'
n_samples = 100
sel1 = 'protein'
sel2 = 'protein'

topo = topo_trajs['MPRO']['topo']
traj = topo_trajs['MPRO']['traj']
u = mda.Universe(topo, traj)
elements = guess_elements(u)
u.add_TopologyAttr('elements', elements)

# create selections for the ligand and protein
ligand_selection = u.select_atoms(sel1)
protein_selection = u.select_atoms(sel2)

# use default interactions
fp = plf.Fingerprint(
    ['Anionic', 'CationPi', 'Cationic', 'EdgeToFace', 'FaceToFace',
     'HBAcceptor', 'HBDonor', 'Hydrophobic', 'MetalAcceptor', 'MetalDonor',
     'PiCation', 'PiStacking', 'VdWContact', 'XBAcceptor', 'XBDonor'],
    count=True, parameters={"VdWContact": {"preset": "rdkit"}})
fp.run(u.trajectory[:n_samples], ligand_selection, protein_selection)
# fp.to_pickle(join(out_dir, 'prolif.pkl'))
df = fp.to_dataframe()
gnl.pickle_to_file(df, join(out_dir, 'prolif.pkl'))
