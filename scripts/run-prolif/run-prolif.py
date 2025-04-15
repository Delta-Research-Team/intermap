# Created by gonzalezroy at 4/15/25


import MDAnalysis as mda
import numpy as np
import prolif as plf
import rdkit


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
#
# =============================================================================
topo = '/home/gonzalezroy/RoyHub/intermap/data/100_frames_trajs/01_MPRO_TRAJECTORIES_INPUTS_DATA_mpro_wt_variants_amarolab/wt_mpro_chainA_rep123.pr5.aligned_CA.not_waters_or_ions.psf'
traj = '/home/gonzalezroy/RoyHub/intermap/data/100_frames_trajs/01_MPRO_TRAJECTORIES_INPUTS_DATA_mpro_wt_variants_amarolab/mpro.dcd'

sel1 = 'protein'
sel2 = 'protein'

u = mda.Universe(topo, traj)
elements = guess_elements(u)
u.add_TopologyAttr('elements', elements)

ligand_selection = u.select_atoms(sel1)
protein_selection = u.select_atoms(sel2)

fp = plf.Fingerprint(count=True)
fp.run(u.trajectory, ligand_selection, protein_selection)
