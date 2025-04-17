# Created by rglez at 12/8/24
import itertools as it
import logging
import shutil
import time
from collections import defaultdict
from pprint import pformat

import MDAnalysis as mda
import numpy as np
import numpy_indexed as npi
import rdkit
from rdkit import Chem

import intermap.commons as cmn
import intermap.interactions.geometry as geom
import intermap.managers.cutoffs as cf

logger = logging.getLogger('InterMapLogger')


def get_periodic_table_info():
    """
    Get the periodic table information
    """
    pt = rdkit.Chem.GetPeriodicTable()
    pt_symbols = {pt.GetElementSymbol(x): x for x in range(1, 119)}
    pt_masses = {x: pt.GetAtomicWeight(x) for x in range(1, 119)}
    real_names = {x.lower(): x for x in pt_symbols}
    return pt_symbols, pt_masses, real_names


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


def ag2rdkit(ag):
    """
    Convert an AtomGroup to an RDKit molecule

    Args:
        ag: AtomGroup object from MDAnalysis

    Returns:
        rdkit_mol: RDKit molecule
    """
    try:
        return ag.convert_to("RDKIT", force=False)
    except AttributeError:
        return ag.convert_to("RDKIT", force=True)


def match_rings(mol):
    aromatic_atoms = set([atom.GetIdx() for atom in mol.GetAtoms() if
                          atom.GetIsAromatic()])

    aromatic_rings = []
    for ring in mol.GetRingInfo().AtomRings():
        if all(atom in aromatic_atoms for atom in ring):
            aromatic_rings.append(list(ring))
    return aromatic_rings


def get_uniques_triads(universe):
    """
    Get the unique residues in the universe taking into account the connected
    residues

    Args:
        universe (mda.Universe): Universe object

    Returns:
        unique_mda_res (dict): Dictionary with the unique residues
        unique_rdmols (dict): Dictionary with the unique RDKit molecules
        unique_idx (dict): Dictionary with the indices of the unique residues
    """
    # Load universe information
    stamp = time.time()
    by_resnames = universe.residues.resnames
    by_resindex = universe.residues.resindices
    bonds = universe.bonds

    # Get connected residues
    connected = defaultdict(list)
    for a1, a2 in bonds:
        r1 = a1.resindex
        r2 = a2.resindex
        if r1 != r2:
            connected[r1].append(r2)
            connected[r2].append(r1)

    # Convert connected triads to rdkit format
    mda_connected = {}
    rdk_connected = {}
    for residue in connected:
        neighbors = {residue}
        neighbors.update(connected[residue])
        triad = universe.residues[list(neighbors)]
        mda_connected[residue] = triad
        rdk_connected[residue] = ag2rdkit(triad.atoms)

    # Get unique disconnected residues
    connected_set = set(connected.keys())
    disconnected = set.difference(set(by_resindex), connected_set)
    uniq_disconnected = defaultdict(list)
    for x in disconnected:
        uniq_disconnected[by_resnames[x]].append(x)

    # Convert disconnected monomers to rdkit format
    mda_disconnected = {}
    rdk_disconnected = {}
    for residue in uniq_disconnected:
        mono = universe.residues[uniq_disconnected[residue][0]]
        mda_disconnected[residue] = mono
        rdk_disconnected[residue] = ag2rdkit(mono.atoms)

    rdkit_ = time.time() - stamp
    logger.info(f"Residues converted to Rdkit format in {rdkit_:.2f} s")
    return (mda_connected, rdk_connected, mda_disconnected, rdk_disconnected,
            uniq_disconnected)


class IndexManager:
    """
    Class to manage the indices of the selections in a trajectory.
    """
    smarts = {
        'hydroph': '[c,s,Br,I,S&v2&H0,$([D3,D4;#6])&!$([#6]~[#7,#8,#9])&!$([#6&X4&H0]);+0]',
        'cations': '[+{1-},$([N&X3&!$([N&X3]-O)]-C=[N&X3&+])]',
        'anions': '[-{1-},$(O=[C,S,P]-[O&-])]',
        'metal_acc': '[O,#7&!$([n&X3])&!$([N&X3]-*=[!#6])&!$([N&X3]-a)&!$([N&X4]),-{1-};!+{1-}]',
        'metal_don': '[#20,#48,#27,#29,#26,#12,#25,#28,#30]',
        'hb_acc': '[#7&!$([nX3])&!$([NX3]-*=[O,N,P,S])&!$([NX3]-[a])&!$([Nv4&+1]),O&!$([OX2](C)C=O)&!$(O(~a)~a)&!$(O=N-*)&!$([O-]-N=O),o+0,F&$(F-[#6])&!$(F-[#6][F,Cl,Br,I])]',
        'xb_acc': '[$([O,S;+0]),$([N;v3,v4&+1]),n+0]-[H]',
        'hb_don': '[$([O,S;+0]),$([N;v2,v3,v4&+1]),n+0]-[H]',
        'xb_don': '[#6,#7,#14,F,Cl,Br,I]-[Cl,Br,I,#85]',
        'rings5': '[a&r]1:[a&r]:[a&r]:[a&r]:[a&r]:1',
        'rings6': '[a&r]1:[a&r]:[a&r]:[a&r]:[a&r]:[a&r]:1',
        'water': ' [OH2]'}

    def __init__(self, args):
        logger.info("Loading trajectory and selections."
                    " Using the following SMARTS patterns:\n\n"
                    f"{pformat(self.smarts)}\n")

        # Get the arguments from the parser
        self.args = args
        self.topo = args.topology
        self.traj = args.trajectory
        self.sel1 = args.selection_1
        self.sel2 = args.selection_2
        self.raw_inters = args.interactions
        self.last_frame = args.last

        # Load the trajectory as a Universe
        (self.universe, self.traj_frames, self.n_atoms,
         self.n_frames) = self.load_traj()

        # Get annotations
        self.annotations = self.get_annotations()

        # Get indices of the selections
        (self.sel_idx, self.s1_idx, self.s2_idx,
         self.overlap) = self.get_selections_indices()

        # Get the names of the atoms
        (self.resconv, self.resid_names, self.atom_names, self.resid_notes,
         self.atom_notes) = self.get_resids_and_names()

        # Get triads (connected) / monomers (disconnected) of residues
        (self.mda_connected, self.rdk_connected, self.mda_disconnected,
         self.rdk_disconnected, self.uniq_disconnected) = get_uniques_triads(
            self.universe)

        # Get VDW radii
        self.vdw_radii = self.get_vdw_radii()

        # Get indices of the 1D interactions
        self.hydroph = self.get_singles('hydroph')
        self.cations = self.get_singles('cations')
        self.anions = self.get_singles('anions')
        self.met_acc = self.get_singles('metal_acc')
        self.met_don = self.get_singles('metal_don')

        # Get indices of the 2D1A interactions
        self.hb_don, self.hb_hydro, self.hb_acc = self.get_doubles(
            'hb_don', 'hb_acc')
        self.xb_don, self.xb_hal, self.xb_acc = self.get_doubles(
            'xb_don', 'xb_acc')

        # Get indices of waters
        self.waters, self.raw_waters = self.get_waters()

        # Get indices of the aromatic interactions
        (self.rings, self.s1_cat, self.s2_cat, self.s1_cat_idx,
         self.s2_cat_idx, self.s1_rings, self.s2_rings, self.s1_rings_idx,
         self.s2_rings_idx, self.s1_aro_idx, self.s2_aro_idx,
         self.xyz_aro_idx) = self.get_aro()

        # Report the interactions detected
        self.inters_requested = self.report()

    def load_traj(self):
        """
        Load the trajectory into a Universe

        Returns:
            universe (mda.Universe): Universe object
        """

        # Load the trajectory
        stamp0 = time.time()
        trajs = [cmn.check_path(x.strip()) for x in self.traj.split(',')]
        universe = mda.Universe(*([self.topo] + trajs))
        masses = [round(x) for x in universe.atoms.masses]
        names = universe.atoms.names

        # Ensure all elements are present
        try:
            elements = universe.atoms.elements.lolita
        except Exception:
            logger.warning("Element information not found in the topology."
                           " Using the atoms name to guess it.")
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
            logger.warning(f'The passed topology does not contain bonds. '
                           f'MDAnalysis will guess them automatically.')
            guessing = time.time()
            universe = mda.Universe(self.topo, self.traj, guess_bonds=True,
                                    vdwradii=radii)
            n_bonds = len(universe.bonds)
            logger.info(
                f" {n_bonds} bonds guessed in {time.time() - guessing:.2f} s")
            any_bond = universe.bonds[0]
        if any_bond is None:
            raise ValueError(
                "No bonds found in topology and MDAnalysis could not guess them.")

        # Remove the hydrogen-hydrogen bonds if any
        universe.add_TopologyAttr('elements', elements)
        stamp1 = time.time()
        are_hh = any_hh_bonds(universe)
        universe.delete_bonds(get_hh_bonds(universe))
        del_time = time.time() - stamp1
        if are_hh:
            logger.warning(
                f"This universe contained H-H bonds. Removed in {del_time:.2f} s")

        # Output load time
        loading = time.time() - stamp0

        # Chunks of frames to analyze
        last = cmn.parse_last_param(self.last_frame, len(universe.trajectory))
        traj_frames = np.arange(self.args.start, last, self.args.stride)

        logger.info(
            f"Trajectory loaded in {loading:.2f} s.\n"
            f" Number of frames to consider (start:last:stride): "
            f"{traj_frames.size} ({self.args.start}:{last}:{self.args.stride})")

        n_atoms = universe.atoms.n_atoms
        n_frames = len(universe.trajectory)

        # Copy the topology to the output dir
        shutil.copy(self.topo, self.args.output_dir)
        return universe, traj_frames, n_atoms, n_frames

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

        uniques = sorted(set(s1_idx).union(set(s2_idx)))
        sel_idx = np.asarray(uniques, dtype=np.int32)
        s1_idx = npi.indices(sel_idx, s1_idx).astype(np.int32)
        s2_idx = npi.indices(sel_idx, s2_idx).astype(np.int32)
        overlap = np.intersect1d(s1_idx, s2_idx).size > 0
        return sel_idx, s1_idx, s2_idx, overlap

    def get_singles(self, identifier):
        """
        Get the indices associated with the single interactions

        Args:
            identifier (str): Identifier of the interaction

        Returns:
            singles: array with the indices of the single atoms
        """

        query = Chem.MolFromSmarts(self.smarts[identifier])
        singles = []

        # Look for the single atoms in the connected residues
        for case in self.rdk_connected:
            tri_mol = self.rdk_connected[case]
            match = [y for x in tri_mol.GetSubstructMatches(query) for y in x]
            if match:
                tri_res = self.mda_connected[case]
                where = tri_res.atoms.resindices[match] == case
                selected = tri_res.atoms.indices[match][where]
                singles.extend(selected)

        # Look for the single atoms in the disconnected residues
        for case in self.rdk_disconnected:
            mono_mol = self.rdk_disconnected[case]
            match = [y for x in mono_mol.GetSubstructMatches(query) for y in x]
            if match:
                for similar in self.uniq_disconnected[case]:
                    mono_res = self.universe.residues[similar]
                    selected = mono_res.atoms.indices[match]
                    singles.extend(selected)

        selected_raw = npi.indices(
            self.sel_idx, np.asarray(singles), missing=-1)
        selected = selected_raw[selected_raw != -1].astype(np.int32)
        return selected

    def get_doubles(self, donor_identifier, acceptor_identifier):
        """
        Get the indices associated with the hydrogen bonds

        Returns:
            hb_D: array with the indices of the donors
            hb_H: array with the indices of the hydrogens
            hb_A: array with the indices of the acceptors
        """

        smart_dx = self.smarts[donor_identifier]
        smart_a = self.smarts[acceptor_identifier]

        hx_A = []
        hx_D = []
        hx_H = []
        query_dh = Chem.MolFromSmarts(smart_dx)
        query_a = Chem.MolFromSmarts(smart_a)

        # Look for D, H, A in the connected residues
        for case in self.rdk_connected:
            tri_mol = self.rdk_connected[case]
            tri_res = self.mda_connected[case]

            match_dh = [x for x in tri_mol.GetSubstructMatches(query_dh)]
            if match_dh:
                where_dh = tri_res.atoms.resindices[match_dh] == case
                selected_dh = tri_res.atoms.indices[match_dh][
                    where_dh.all(axis=1)]

                for i, x in enumerate(selected_dh):
                    hx_D.append(x[0])
                    hx_H.append(x[1])

            match_a = [x for x in tri_mol.GetSubstructMatches(query_a)]
            if match_a:
                where_a = tri_res.atoms.resindices[match_a] == case
                selected_a = tri_res.atoms.indices[match_a][where_a]
                hx_A.extend(selected_a)

        # Look for D, H, A in the disconnected residues
        for case in self.rdk_disconnected:
            mono_mol = self.rdk_disconnected[case]
            match_dh = [x for x in mono_mol.GetSubstructMatches(query_dh)]
            if match_dh:
                for similar in self.uniq_disconnected[case]:
                    mono_res = self.universe.residues[similar]
                    selected_dh = mono_res.atoms.indices[match_dh]
                    for i, x in enumerate(selected_dh):
                        hx_D.append(x[0])
                        hx_H.append(x[1])

            match_a = [x for x in mono_mol.GetSubstructMatches(query_a)]
            if match_a:
                for similar in self.uniq_disconnected[case]:
                    mono_res = self.universe.residues[similar]
                    selected_a = mono_res.atoms.indices[match_a]
                    hx_A.extend(selected_a[:, 0])

        # Filter the indices
        hx_D_raw = npi.indices(self.sel_idx, np.asarray(hx_D), missing=-1)
        hx_H_raw = npi.indices(self.sel_idx, np.asarray(hx_H), missing=-1)
        hx_A_raw = npi.indices(self.sel_idx, np.asarray(hx_A), missing=-1)

        hx_D = hx_D_raw[hx_D_raw != -1].astype(np.int32)
        hx_H = hx_H_raw[hx_H_raw != -1].astype(np.int32)
        hx_A = np.unique(hx_A_raw[hx_A_raw != -1]).astype(np.int32)
        assert len(hx_D) == len(
            hx_H), f"Donors ({hx_D.size}) and Hydrogens ({hx_H.size}) do not match"
        return hx_D, hx_H, hx_A

    def get_rings(self):
        """
        Get the indices associated with the aromatic rings

        Returns:
            rings: List with the indices of the aromatic rings
        """
        patterns = [self.smarts['rings5'], self.smarts['rings6']]

        rings = []
        # Look for rings in the connected residues
        for case in self.rdk_connected:
            tri_mol = self.rdk_connected[case]
            tri_res = self.mda_connected[case]

            for smart in patterns:
                query = Chem.MolFromSmarts(smart)
                match = [list(x) for x in tri_mol.GetSubstructMatches(query)]
                if match:
                    where_dh = tri_res.atoms.resindices[match] == case
                    selected_dh = tri_res.atoms.indices[match][where_dh]
                    if selected_dh.any():
                        rings.append(selected_dh)

        # Look for rings in the disconnected residues
        for case in self.rdk_disconnected:
            mono_mol = self.rdk_disconnected[case]
            for smart in patterns:
                query = Chem.MolFromSmarts(smart)
                match = [list(x) for x in mono_mol.GetSubstructMatches(query)]
                if match:
                    for similar in self.uniq_disconnected[case]:
                        mono_res = self.universe.residues[similar]
                        selected = mono_res.atoms.indices[match]
                        if selected.any():
                            for ring in selected:
                                rings.append(ring)

        # Pad rings with -1
        padded_rings = np.full((len(rings), 7), dtype=np.int32, fill_value=-1)
        for i, ring in enumerate(rings):
            padded_rings[i, :len(ring)] = ring
            padded_rings[i, -1] = len(ring)

        sel_rings = padded_rings.copy()
        for i, ring in enumerate(padded_rings):
            r = ring[:ring[-1]]
            sel_rings[i, :ring[-1]] = npi.indices(self.sel_idx, r, missing=-1)
        return sel_rings

    def get_aro(self):
        """
        Get the indices associated with the aromatic interactions

        Returns:
            s1_cat (ndarray): Indices of the atoms in selection 1 that are cations
            s2_cat (ndarray): Indices of the atoms in selection 2 that are cations

            s1_rings (ndarray): Indices of the atoms in selection 1 that are in rings
            s2_rings (ndarray): Indices of the atoms in selection 2 that are in rings
            s1_rings_idx (ndarray): Indices of the atoms in selection 1 that are in rings
            s2_rings_idx (ndarray): Indices of the atoms in selection 2 that are in rings
            s1_aro_idx (ndarray): Indices of the atoms in selection 1 that are in aromatic interactions
            s2_aro_idx (ndarray): Indices of the atoms in selection 2 that are in aromatic interactions
            xyz_aro_idx (ndarray): Indices of the atoms in the universe that are in aromatic interactions
        """
        rings_raw = self.get_rings()
        rings = rings_raw[rings_raw[:, 0] != -1]
        s1_cat = self.s1_idx[geom.isin(self.s1_idx, self.cations)]
        s2_cat = self.s2_idx[geom.isin(self.s2_idx, self.cations)]
        s1_rings = rings[geom.isin(rings[:, 0], self.s1_idx)]
        s2_rings = rings[geom.isin(rings[:, 0], self.s2_idx)]

        n0 = s1_cat.size + s2_cat.size
        n1 = n0 + s1_rings.shape[0]
        n2 = n1 + s2_rings.shape[0]

        s1_cat_idx = np.arange(0, s1_cat.size, dtype=np.int32)
        s2_cat_idx = np.arange(s1_cat.size, n0, dtype=np.int32)
        s1_rings_idx = np.arange(n0, n1, dtype=np.int32)
        s2_rings_idx = np.arange(n1, n2, dtype=np.int32)
        s1_aro_idx = np.concatenate((s1_cat_idx, s1_rings_idx)).astype(
            np.int32)
        s2_aro_idx = np.concatenate((s2_cat_idx, s2_rings_idx)).astype(
            np.int32)
        xyz_aro_idx = np.concatenate((s1_cat, s2_cat, s1_rings[:, 0],
                                      s2_rings[:, 0])).astype(np.int32)

        return (
            rings, s1_cat, s2_cat, s1_cat_idx, s2_cat_idx, s1_rings, s2_rings,
            s1_rings_idx, s2_rings_idx, s1_aro_idx, s2_aro_idx, xyz_aro_idx)

    def get_waters(self):
        """
        Get the indices associated with the water molecules

        Returns:
            waters: List with the indices of the water molecules
        """
        query = Chem.MolFromSmarts(self.smarts['water'])
        query = Chem.AddHs(query)
        waters = []

        # Look for the water molecules in the disconnected residues
        for case in self.rdk_disconnected:
            mono_mol = self.rdk_disconnected[case]
            match = [y for x in mono_mol.GetSubstructMatches(query) for y in x]
            if match:
                for similar in self.uniq_disconnected[case]:
                    mono_res = self.universe.residues[similar]
                    selected = mono_res.atoms.indices[match]
                    waters.extend(selected)

        selected_raw = npi.indices(
            self.sel_idx, np.asarray(waters), missing=-1)
        selected = selected_raw[selected_raw != -1].astype(np.int32)

        return selected, waters

    def get_max_vdw_dist(self):
        """
        Get the maximum van der Waals distance between the atoms in the
         selections.

        Returns:
            max_vdw (float): Maximum van der Waals distance between the atoms
        """

        s1_elements = set(self.universe.atoms[self.s1_idx].elements)
        s2_elements = set(self.universe.atoms[self.s2_idx].elements)

        product = it.product(s1_elements, s2_elements)
        unique_pairs = set(tuple(sorted((a, b))) for a, b in product)
        unique_elements = set([y for x in unique_pairs for y in x])

        pt = rdkit.Chem.GetPeriodicTable()
        radii = {x: pt.GetRvdw(x) for x in unique_elements}

        radiis = [radii[pair[0]] + radii[pair[1]] for pair in unique_pairs]
        max_vdw = max(radiis)
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
        all_radii = radii.astype(np.float32)
        return all_radii[self.sel_idx]

    def get_interactions(self):
        """
        Get the possible interactions between the selections in the trajectory

        Returns:

        """

        # Check the lenght of the selections
        len_s1, len_s2 = len(self.s1_idx), len(self.s2_idx)
        len_hp = len(self.hydroph)
        len_an, len_ca = len(self.anions), len(self.cations)
        len_ma, len_md = len(self.met_acc), len(self.met_don)
        len_hbd, len_hba, len_hbh = len(self.hb_don), len(self.hb_acc), len(
            self.hb_hydro)
        len_xbd, len_xba, len_xbh = len(self.xb_don), len(self.xb_acc), len(
            self.xb_hal)
        len_rings = len(self.rings)

        # Gather the interactions constraints
        is_possible = {
            'CloseContact': (len_s1 > 0) and (len_s2 > 0),
            'VdWContact': (len_s1 > 0) and (len_s2 > 0),
            'Hydrophobic': len_hp > 0,
            'Anionic': (len_an > 0) and (len_ca > 0),
            'Cationic': (len_an > 0) and (len_ca > 0),
            'MetalAcceptor': (len_ma > 0) and (len_md > 0),
            'MetalDonor': (len_ma > 0) and (len_md > 0),
            'HBDonor': (len_hbd > 0) and (len_hba > 0) and (len_hbh > 0),
            'HBAcceptor': (len_hbd > 0) and (len_hba > 0) and (len_hbh > 0),
            'XBDonor': (len_xbd > 0) and (len_xba > 0) and (len_xbh > 0),
            'XBAcceptor': (len_xbd > 0) and (len_xba > 0) and (len_xbh > 0),
            'PiStacking': len_rings > 0,
            'PiCation': (len_rings > 0) and (len_ca > 0),
            'CationPi': (len_rings > 0) and (len_ca > 0),
            'FaceToFace': len_rings > 0,
            'EdgeToFace': len_rings > 0,
            'WaterBridge': len(self.waters) > 0,
        }

        raw = self.raw_inters
        if isinstance(raw, str):
            requested = cf.interactions
        else:
            requested = raw

        to_compute = []
        to_skip = []
        for x in requested:
            if not is_possible[x]:
                to_skip.append(x)
            else:
                to_compute.append(x)

        if to_skip:
            logger.warning(
                f"Skipping the following interactions because there "
                f"are no atoms to compute them under the current "
                f"selections for the passed system:\n"
                f" {pformat(to_skip)}")
        return to_compute

    def get_resids_and_names(self):
        """
        Get the residue IDs and names of the atoms in the selections

        Returns:
            resids: array with the residue IDs
            names: list with the names of the atoms
        """
        atnames = self.universe.atoms.names[self.sel_idx]
        resnames = self.universe.atoms.resnames[self.sel_idx]
        resids = self.universe.atoms.resids[self.sel_idx]
        atindex = self.universe.atoms.indices[self.sel_idx]
        resindex = self.universe.atoms.resindices[self.sel_idx].astype(
            np.int32)

        annotations = self.annotations
        if annotations:
            at_annots = {x: k for k, v in annotations.items() for x in v}
            res_annots = {
                x: at_annots.get(i, '') for i, x in
                enumerate(self.universe.atoms.resindices)}
        else:
            at_annots = {}
            res_annots = {}

        atom_names = {
            i: f"{resnames[i]}_{resids[i]}_{resindex[i]}_{atnames[i]}_{atindex[i]}"
            for i, x in enumerate(self.sel_idx)}

        resid_names = {resindex[i]: f"{resnames[i]}_{resids[i]}_{resindex[i]}"
                       for i, x in
                       enumerate(self.sel_idx)}

        return resindex, resid_names, atom_names, res_annots, at_annots

    def report(self):
        """
        Report the interactions detected

        Returns:
            inters_requested: List with the interactions requested
        """
        logger.info(
            f'The following selected atoms were detected and classified:\n\n'
            f' Total number of atoms: {self.n_atoms}\n'
            f" Number of selected atoms: {self.sel_idx.size}\n"
            f"   In Selection 1 ({self.sel1}): {len(self.s1_idx)}\n"
            f"   In Selection 2 ({self.sel2}): {len(self.s2_idx)}\n"
            f"   Hydrophobic: {len(self.hydroph)}\n"
            f"   Cations: {len(self.cations)}\n"
            f"   Anions: {len(self.anions)}\n"
            f"   Metal acceptors: {len(self.met_acc)}\n"
            f"   Metal donors: {len(self.met_don)}\n"
            f"   Hydrogen bond donors: {len(self.hb_don)}\n"
            f"   Hydrogen bond hydrogens: {len(self.hb_hydro)}\n"
            f"   Hydrogen bond acceptors: {len(self.hb_acc)}\n"
            f"   Halogen bond donors: {len(self.xb_don)}\n"
            f"   Halogens: {len(self.xb_hal)}\n"
            f"   Halogen bond acceptors: {len(self.xb_acc)}\n"
            f"   Aromatic rings: {len(self.rings)}\n"
            f"   Water molecules: {len(self.waters)}\n")

        # Possible interactions
        inters_requested = self.get_interactions()
        return inters_requested

    def get_annotations(self):
        """"
        Get the atoms for annotations specified

        Returns:
            annotations: dictionary with the annotations
        """
        annot_file = self.args.annotations
        universe = self.universe

        annotations = defaultdict(set)
        cumulative = set()

        if annot_file:
            with open(annot_file) as f:
                lines = f.readlines()
            lines = [x.strip() for x in lines if x not in ['', '\n', '#']]

            for line in lines:
                if line.startswith('#'):
                    continue

                try:
                    name, value = line.split('=')
                except ValueError:
                    raise ValueError(
                        f'Error parsing the annotations file. The line: {line}'
                        f' is not valid. The expected format is: name = value.'
                        f' While the name can be any string, the value must be'
                        f' a valid MDAnalysis selection string.')

                indices = universe.select_atoms(value).indices
                annotations[name.strip()].update(indices)
                if len(indices) == 0:
                    logger.warning(
                        f"Selection {value} for the key {name} returned no"
                        f" atoms in the topology {self.args.topology}. "
                        f"Please check the selection.")

                intersect = set(indices).intersection(cumulative)
                cumulative.update(indices)
                if len(intersect) > 0:
                    raise ValueError(
                        f"Error parsing the annotations file. The selection "
                        f"{value} for the key {name} overlaps with one or"
                        f" more of the previously defined selections.")

        return annotations

# =============================================================================
# import intermap.managers.config as conf
# from argparse import Namespace
#
# config = conf.ConfigManager(mode='debug')
# args = Namespace(**config.config_args)
# self = IndexManager(args)
