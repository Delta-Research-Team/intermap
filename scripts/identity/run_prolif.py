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
topo_traj_path = '/home/gonzalezroy/RoyHub/intermap/data/100_frames_trajs/'
out_dir = f'/home/gonzalezroy/RoyHub/intermap/scripts/identity/runs/prolif'

# =============================================================================

# Find the topology and trajectory files
topo_trajs = gnl.recursive_defaultdict()
for case in os.listdir(topo_traj_path):
    case_path = os.path.join(topo_traj_path, case)
    for file in os.listdir(case_path):
        file_path = os.path.join(case_path, file)
        if file_path.endswith('.dcd'):
            topo_trajs[case]['traj'] = file_path
        else:
            topo_trajs[case]['topo'] = file_path

# Define the selections
selections = {
    '01_MPRO_TRAJECTORIES_INPUTS_DATA_mpro_wt_variants_amarolab': {
        'sel1': 'protein',
        'sel2': 'protein'},
    '02_IgG3': {
        'sel1': 'resid 291:1718',
        'sel2': 'resid 1:290'},
    '03_p53_fl_p53_Anton': {
        'sel1': 'nucleic',
        'sel2': 'protein'},
    '04_NSP13_DESRES-Trajectory_sarscov2-13795965-no-water': {
        'sel1': 'nucleic',
        'sel2': 'protein'},
    '05_PAO1': {
        'sel1': 'protein',
        'sel2': 'protein'},
    '06_spike_TRAJECTORIES_spike_open_prot_glyc_amarolab': {
        'sel1': 'protein',
        'sel2': 'not protein'},
    '07_ace2_TRAJECTORIES_ace2_rbd_prot_glyc_memb_amarolab': {
        'sel1': 'protein',
        'sel2': 'not protein'},
    '08_Nucleosome': {
        'sel1': 'nucleic',
        'sel2': 'protein'},
}

if __name__ == '__main__':
    for case in topo_trajs:
        print(f'Processing {case}')
        topo = topo_trajs[case]['topo']
        traj = topo_trajs[case]['traj']
        sel1 = selections[case]['sel1']
        sel2 = selections[case]['sel2']

        # create the output directory
        out_dir_case = join(out_dir, case)
        os.makedirs(out_dir_case, exist_ok=True)

        # load the topology and trajectory
        u = mda.Universe(topo, traj)
        elements = guess_elements(u)
        u.add_TopologyAttr('elements', elements)

        # create selections for the ligand and protein
        ligand_selection = u.select_atoms(sel1)
        protein_selection = u.select_atoms(sel2)

        # use default interactions
        fp = plf.Fingerprint(
            ['Anionic', 'CationPi', 'Cationic', 'EdgeToFace', 'FaceToFace',
             'HBAcceptor', 'HBDonor', 'Hydrophobic', 'MetalAcceptor',
             'MetalDonor',
             'PiCation', 'PiStacking', 'VdWContact', 'XBAcceptor', 'XBDonor'],
            count=True, parameters={"VdWContact": {"preset": "rdkit"}})
        fp.run(u.trajectory, ligand_selection, protein_selection)
        df = fp.to_dataframe()
        gnl.pickle_to_file(df, join(out_dir, 'prolif.pkl'))
