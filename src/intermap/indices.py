# Created by rglez at 12/8/24
import itertools as it
import time
from collections import defaultdict

import MDAnalysis as mda
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import rdchem

import intermap.smarts as smarts


def get_hh_bonds(universe):
    """
    Get the hydrogen-hydrogen bonds in the Universe

    Returns:
        bonds (list): List of hydrogen-hydrogen bonds
    """
    bonds = universe.bonds
    for bond in bonds:
        atom1, atom2 = bond
        if (atom1.element == 'H') and (atom2.element == 'H'):
            yield bond


def get_uniques(universe):
    """
    Get the unique residues in the universe

    Args:
        universe (mda.Universe): Universe object

    Returns:
        unique_mda_res (dict): Dictionary with the unique residues
        unique_rdmols (dict): Dictionary with the unique RDKit molecules
        unique_idx (dict): Dictionary with the indices of the unique residues
    """
    by_atnames = universe.residues.names
    by_resnames = universe.residues.resnames

    # Find the unique residues by atom names
    uniques = defaultdict(list)
    frozen_labels = {}
    for i, names in enumerate(by_atnames):
        frozen_atnames = frozenset(names)
        uniques[frozen_atnames].append(i)
        frozen_labels[frozen_atnames] = by_resnames[i]
    unique_indices = list(uniques.values())

    # Parse the unique residues
    unique_mda_res = {}
    unique_rdmols = {}
    unique_idx = {}
    for i, indices in enumerate(unique_indices, 1):
        # Get the label of representative residues
        representative = indices[0]
        res = universe.residues[representative]
        label = f'{res.resname}_{i}'

        # Get the unique mda residues
        unique_mda_res[label] = res

        # Convert the residue to an RDKit molecule
        try:
            unique_rdmols[label] = res.atoms.convert_to("RDKIT", force=True)
        except rdchem.AtomValenceException:
            raise rdchem.AtomValenceException()

        # Get the indices of the unique residues
        unique_idx[label] = indices
    return unique_mda_res, unique_rdmols, unique_idx


class IndexManager:
    """
    Class to manage the indices of the selections in a trajectory.
    """

    def __init__(self, topo, traj, sel1, sel2, interactions):
        # Parse the arguments
        self.topo = topo
        self.traj = traj
        self.sel1 = sel1
        self.sel2 = sel2
        self.interactions = interactions
        self.smarts_patterns = smarts.ProlifSmarts().interactions

        # Initialize the trajectory attributes
        self.unique_residues = None
        self.unique_rdmols = None
        self.universe, self.renamed_universe = self.load_traj()
        self.radii = self.get_vdw_radii()

        # Get the indices of the selections
        self.sel1_idx, self.sel2_idx = self.get_selections_indices()

        self.hydroph = self.get_singles('hydroph')
        self.cations = self.get_singles('cations')
        self.anions = self.get_singles('anions')
        self.metal_acc = self.get_singles('metal_acc')
        self.metal_don = self.get_singles('metal_don')

        self.hb_D, self.hb_H, self.hb_A = self.get_doubles('hb_don', 'hb_acc')
        self.xb_D, self.xb_H, self.xb_A = self.get_doubles('xb_don', 'xb_acc')

        self.rings = self.get_rings()

    def load_traj(self):
        """
        Load the trajectory into a Universe

        Returns:
            universe (mda.Universe): Universe object
        """
        # Load the trajectory
        universe = mda.Universe(self.topo, self.traj)

        # Remove the hydrogen-hydrogen bonds
        t1 = time.time()
        universe.delete_bonds(get_hh_bonds(universe))
        print(f"Just deleting bonds: {time.time() - t1:.2f} s")

        # Get the unique residues, RDKit molecules and indices
        unique_residues, unique_rdmols, unique_idx = get_uniques(universe)
        self.unique_residues = unique_residues
        self.unique_rdmols = unique_rdmols

        # Change the resnames of the universe to the unique labels
        old_resnames = universe.residues.resnames
        new_resnames = old_resnames.copy()
        new_universe = universe.copy()
        for label in unique_idx:
            new_resnames[unique_idx[label]] = label
        new_universe.add_TopologyAttr('resname', new_resnames)
        return universe, new_universe

    def get_selections_indices(self):
        """
        Get the indices of the selections in the trajectory.

        Returns:
            s1_idx: array with the indices of the atoms in selection 1
            s2_idx: array with the indices of the atoms in selection 2
        """
        s1_idx = self.universe.select_atoms(self.sel1).indices.astype(np.int32)
        s2_idx = self.universe.select_atoms(self.sel2).indices.astype(np.int32)

        if len(s1_idx) == 0:
            raise ValueError("No atoms found for selection 1")
        if len(s2_idx) == 0:
            raise ValueError("No atoms found for selection 2")
        return s1_idx, s2_idx

    def get_singles(self, identifier):
        """
        Get the indices associated with the single interactions

        Args:
            identifier (str): Identifier of the interaction

        Returns:
            singles: array with the indices of the single atoms
        """
        smart = self.smarts_patterns[identifier]
        singles = []
        for case in self.unique_rdmols:
            res = self.unique_residues[case]
            atom_names = res.atoms.names
            mol = self.unique_rdmols[case]

            query = Chem.MolFromSmarts(smart)
            match = [list(x) for x in mol.GetSubstructMatches(query)]
            match_names = [' '.join(atom_names[x]) for x in match]
            for sel in match_names:
                hydrophs = 'resname {} and name {}'.format(case, sel)
                idx = self.renamed_universe.select_atoms(hydrophs).indices
                singles.extend(idx)
        return np.array(singles).astype(np.int32)

    def get_doubles(self, donor_identifier, acceptor_identifier):
        """
        Get the indices associated with the hydrogen bonds

        Returns:
            hb_D: array with the indices of the donors
            hb_H: array with the indices of the hydrogens
            hb_A: array with the indices of the acceptors
        """
        smarts = self.smarts_patterns
        smart_dx = smarts[donor_identifier]
        smart_a = smarts[acceptor_identifier]

        hx_D = []
        hx_H = []
        hx_A = []
        for case in self.unique_rdmols:
            res = self.unique_residues[case]
            atom_names = res.atoms.names
            mol = self.unique_rdmols[case]

            # Find the donors and hydrogens
            query_dh = Chem.MolFromSmarts(smart_dx)
            match_dh = [list(x) for x in mol.GetSubstructMatches(query_dh)]
            match_names_dh = [' '.join(atom_names[x]) for x in match_dh]
            for sel in match_names_dh:
                d, x = sel.split()
                donors = 'resname {} and name {}'.format(case, d)
                hydros = 'resname {} and name {}'.format(case, x)
                idx1 = self.renamed_universe.select_atoms(donors).indices
                idx2 = self.renamed_universe.select_atoms(hydros).indices
                assert len(idx1) == len(idx2)
                hx_D.extend(idx1)
                hx_H.extend(idx2)

            # Find the acceptors
            query_a = Chem.MolFromSmarts(smart_a)
            match_a = [list(x) for x in mol.GetSubstructMatches(query_a)]
            match_names_a = [' '.join(atom_names[x]) for x in match_a]
            for a in match_names_a:
                acc = 'resname {} and name {}'.format(case, a)
                idx = self.renamed_universe.select_atoms(acc).indices
                hx_A.extend(idx)

        hx_D = np.array(hx_D).astype(np.int32)
        hx_H = np.array(hx_H).astype(np.int32)
        hx_A = np.array(hx_A).astype(np.int32)
        return hx_D, hx_H, hx_A

    def get_rings(self):
        """
        Get the indices associated with the aromatic rings

        Returns:
            rings: List with the indices of the aromatic rings
        """
        patterns = [self.smarts_patterns['rings5'],
                    self.smarts_patterns['rings6']]

        rings = []
        for case in self.unique_rdmols:
            res = self.unique_residues[case]
            atom_names = res.atoms.names
            mol = self.unique_rdmols[case]

            for smart in patterns:
                query = Chem.MolFromSmarts(smart)
                match = [list(x) for x in mol.GetSubstructMatches(query)]
                match_names = [' '.join(atom_names[x]) for x in match]
                for sel in match_names:
                    ring = 'resname {} and name {}'.format(case, sel)
                    idx = self.renamed_universe.select_atoms(ring).indices
                    rings.extend(idx.reshape(-1, len(sel.split())))

        padded_rings = np.full((len(rings), 7), dtype=int, fill_value=-1)
        for i, ring in enumerate(rings):
            padded_rings[i, :len(ring)] = ring
            padded_rings[i, -1] = len(ring)
        return padded_rings.astype(np.int32)

    def get_max_vdw_dist(self):
        """
        Get the maximum van der Waals distance between the atoms in the
         selections.

        Returns:
            max_vdw (float): Maximum van der Waals distance between the atoms
        """

        s1_elements = set(self.universe.atoms[self.sel1_idx].elements)
        s2_elements = set(self.universe.atoms[self.sel2_idx].elements)

        product = it.product(s1_elements, s2_elements)
        unique_pairs = set(tuple(sorted((a, b))) for a, b in product)
        unique_elements = set([y for x in unique_pairs for y in x])

        pt = rdkit.Chem.GetPeriodicTable()
        radii = {x: pt.GetRvdw(x) for x in unique_elements}

        radiis = [radii[pair[0]] + radii[pair[1]] for pair in unique_pairs]
        max_vdw = max(radiis) / 10
        return np.float32(max_vdw)

    def get_vdw_radii(self):
        """
        Get the van der Waals radii of the atoms in the universe

        Returns:
            radii (ndarray): Array with the van der Waals radii of the atoms
        """
        elements = self.universe.atoms.elements
        pt = rdkit.Chem.GetPeriodicTable()
        radii = np.array([pt.GetRvdw(e) for e in elements])
        return radii.astype(np.float32)

# %% ==========================================================================
#
# =============================================================================
# import time
#
# topo = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1.prmtop'
# traj = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1_sk100.nc'
#
# sel1 = "nucleic"
# sel2 = "resname WAT"
# interactions = 'all'
#
# start = time.time()
# self = IndexManager(topo, traj, sel1, sel2, interactions)
# print(f"Elapsed time: {time.time() - start:.2f} s")
