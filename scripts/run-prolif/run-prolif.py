# Created by gonzalezroy at 4/15/25
from collections import defaultdict
from os.path import join

import MDAnalysis as mda
import numpy as np
import prolif as plf
import rdkit
from rgpack import generals as gnl


def guess_bonds(universe, topo, traj):
    # Guess bonds if not present
    pt = rdkit.Chem.GetPeriodicTable()
    radii = {e: pt.GetRvdw(e) for e in elements}

    try:
        any_bond = universe.bonds[0]
    except:
        universe = mda.Universe(topo, traj, guess_bonds=True,
                                vdwradii=radii)
        any_bond = universe.bonds[0]
    if any_bond is None:
        raise ValueError(
            "No bonds found in topology and MDAnalysis could not guess them.")
    return universe


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


def fix_hh_and_hvalence(universe):
    """
    Delete H-H bonds from the universe

    Args:
        universe (MDAnalysis.Universe): Universe object
    """
    bonds = universe.bonds
    neighbors = defaultdict(list)

    hh = []
    for a1, a2 in bonds:
        a1_idx = a1.index
        a2_idx = a2.index
        if (a1.element == 'H') and (a2.element == 'H'):
            hh.append((a1_idx, a2_idx))
            continue
        if a1.element == 'H':
            neighbors[a1_idx].append(a2_idx)
        if a2.element == 'H':
            neighbors[a2_idx].append(a1_idx)

    # Remove H-H bonds
    universe.delete_bonds(hh)

    offending = {k: v for k, v in neighbors.items() if len(v) > 1}
    hv = []
    for k, v in offending.items():
        distances = []
        for j in v:
            distance = calc_dist(universe.atoms[k].position,
                                 universe.atoms[j].position)
            distances.append(distance)
        argmin = np.argmin(distances)
        v.pop(argmin)
        for x in v:
            hv.append((k, x))
    universe.delete_bonds(hv)
    return universe


def any_hh_bonds(universe):
    """
    Get the hydrogen-hydrogen bonds in the Universe

    Returns:
        bonds (list): List of hydrogen-hydrogen bonds
    """
    bonds = universe.bonds
    for bond in bonds:
        atom1, atom2 = bond
        if (atom1.element == 'H') and (atom2.element == 'H'):
            return True
    return False


def load_traj(topo, trajs):
    """
    Load the trajectory into a Universe

    Returns:
        universe (mda.Universe): Universe object
    """

    # Load the trajectory
    universe = mda.Universe(*([topo] + trajs))
    masses = [round(x) for x in universe.atoms.masses]
    names = universe.atoms.names

    # Ensure all elements are present
    try:
        elements = universe.atoms.elements.lolita
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

    # Guess bonds if not present
    pt = rdkit.Chem.GetPeriodicTable()
    radii = {e: pt.GetRvdw(e) for e in elements}

    try:
        any_bond = universe.bonds[0]
    except:
        universe = mda.Universe(topo, traj, guess_bonds=True,
                                vdwradii=radii)
        any_bond = universe.bonds[0]
    if any_bond is None:
        raise ValueError(
            "No bonds found in topology and MDAnalysis could not guess them.")

    # Remove the hydrogen-hydrogen bonds if any
    universe.add_TopologyAttr('elements', elements)

    universe = fix_hh_and_hvalence(universe)

    return universe


# =============================================================================
#
# =============================================================================

topo = plf.datafiles.TOP
traj = plf.datafiles.TRAJ

sel1 = 'protein'
sel2 = 'resname LIG'

# trajs = [traj]
# u = load_traj(topo, trajs)
u = mda.Universe(topo, traj)

ligand_selection = u.select_atoms(sel1)
protein_selection = u.select_atoms(sel2)

# fp = plf.Fingerprint(count=True)
fp = plf.Fingerprint(count=True,
                     parameters={"VdWContact": {"preset": "rdkit"}})
fp.run(u.trajectory, ligand_selection, protein_selection)
df = fp.to_dataframe()
out_name = join(
    '/media/gonzalezroy/Roy2TB/RoyData/intermap/IDENTITY-FINAL/prolif/toy-250/',
    'fingerprint.dfpickle')
gnl.pickle_to_file(df, out_name)
