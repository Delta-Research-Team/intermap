# Created by rglez at 12/8/24
import itertools as it
import logging
import time
from collections import defaultdict
from pprint import pformat

import MDAnalysis as mda
import numpy as np
import numpy_indexed as npi
import rdkit
from rdkit import Chem

import intermap.cutoffs as cf

logger = logging.getLogger('InterMapLogger')


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
    try:
        return ag.convert_to("RDKIT", force=False)
    except AttributeError:
        return ag.convert_to("RDKIT", force=True)


def get_uniques_triads(universe):
    """
    Get the unique residues in the universe

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
    logger.info(f"Time to convert residues to Rdkit format: {rdkit_:.2f} s")
    return (mda_connected, rdk_connected, mda_disconnected, rdk_disconnected,
            uniq_disconnected)


class IndexManager:
    """
    Class to manage the indices of the selections in a trajectory.
    """
    smarts = {
        'hydroph': '[c,s,Br,I,S&H0&v2,$([D3,D4;#6])&!$([#6]~[#7,#8,#9])&!$([#6&X4&H0]);+0]',
        'cations': '[+{1-},$([N&X3&!$([N&X3]-O)]-C=[N&X3&+])]',
        'anions': '[-{1-},$(O=[C,S,P]-[O&-])]',
        'metal_acc': '[O,#7&!$([n&X3])&!$([N&X3]-*=[!#6])&!$([N&X3]-a)&!$([N&X4]),-{1-};!+{1-}]',
        'metal_don': '[#20,#48,#27,#29,#26,#12,#25,#28,#30]',
        'hb_acc': '[#7&!$([n&X3])&!$([N&X3]-*=[O,N,P,S])&!$([N&X3]-a)&!$([N&v4&+]),O&!$([O&X2](C)C=O)&!$(O(~a)~a)&!$(O=N-*)&!$([O&-]-N=O),o&+0,F&$(F-[#6])&!$(F-[#6][F,Cl,Br,I])]',
        'xb_acc': '[#7,#8,P,S,#34,#52,a;!+{1-}]!#*',
        'hb_don': '[$([O,S;+0]),$([N;v2,v3,v4&+1]),n+0]-[H]',
        'xb_don': '[#6,#7,#14,F,Cl,Br,I]-[Cl,Br,I,#85]',
        'rings5': '[a&r]1:[a&r]:[a&r]:[a&r]:[a&r]:1',
        'rings6': '[a&r]1:[a&r]:[a&r]:[a&r]:[a&r]:[a&r]:1'}

    def __init__(self, topo, traj, sel1, sel2, interactions):

        # Parse the arguments
        self.topo = topo
        self.traj = traj
        self.sel1 = sel1
        self.sel2 = sel2
        self.raw_inters = interactions

        # Load the trajectory as a Universe
        self.universe = self.load_traj()
        self.n_atoms = self.universe.atoms.n_atoms
        self.n_frames = len(self.universe.trajectory)
        logger.info(f'Total number of atoms: {self.n_atoms}\n'
                    f'Total number of frames: {self.n_frames}')

        # Make the general selections
        sel1_idx, sel2_idx = self.get_selections_indices()
        uniques = sorted(set(sel1_idx).union(set(sel2_idx)))
        self.sel_idx = np.asarray(uniques, dtype=np.int32)
        self.sel1_idx = npi.indices(self.sel_idx, sel1_idx).astype(np.int32)
        self.sel2_idx = npi.indices(self.sel_idx, sel2_idx).astype(np.int32)

        # Get VDW radii
        all_radii = self.get_vdw_radii()
        self.radii = all_radii[self.sel_idx]

        # Get triads (connected) / monomers (disconnected) of residues
        (self.mda_connected, self.rdk_connected, self.mda_disconnected,
         self.rdk_disconnected, self.uniq_disconnected) = get_uniques_triads(
            self.universe)

        # Single interactions
        self.hydroph = self.get_singles('hydroph')
        self.cations = self.get_singles('cations')
        self.anions = self.get_singles('anions')
        self.metal_acc = self.get_singles('metal_acc')
        self.metal_don = self.get_singles('metal_don')

        # Double interactions
        self.hb_D, self.hb_H, self.hb_A = self.get_doubles('hb_don', 'hb_acc')
        self.xb_D, self.xb_H, self.xb_A = self.get_doubles('xb_don', 'xb_acc')

        # Rings
        rings = self.get_rings()
        sel_rings = rings.copy()
        for i, ring in enumerate(rings):
            r = ring[:ring[-1]]
            sel_rings[i, :ring[-1]] = npi.indices(self.sel_idx, r, missing=-1)
        self.rings = sel_rings[sel_rings[:, 0] != -1]
        logger.debug(f"Detected atom types:\n"
                     f"In Selection 1 ({self.sel1}): {len(self.sel1_idx)}\n"
                     f"In Selection 2 ({self.sel2}): {len(self.sel2_idx)}\n"
                     f"Hydrophobic: {len(self.hydroph)}\n"
                     f"Cations: {len(self.cations)}\n"
                     f"Anions: {len(self.anions)}\n"
                     f"Metal acceptors: {len(self.metal_acc)}\n"
                     f"Metal donors: {len(self.metal_don)}\n"
                     f"Hydrogen bond donors: {len(self.hb_D)}\n"
                     f"Hydrogen bond hydrogens: {len(self.hb_H)}\n"
                     f"Hydrogen bond acceptors: {len(self.hb_A)}\n"
                     f"Halogen bond donors: {len(self.xb_D)}\n"
                     f"Halogens: {len(self.xb_H)}\n"
                     f"Halogen bond acceptors: {len(self.xb_A)}\n"
                     f"Aromatic rings: {len(self.rings)}")
        # Possible interactions
        self.interactions = self.get_interactions()

        # Check the SMARTS patterns
        logger.debug(f"Using the following SMARTS patterns:\n"
                     f" {pformat(self.smarts)}")

    def load_traj(self):
        """
        Load the trajectory into a Universe

        Returns:
            universe (mda.Universe): Universe object
        """

        # Load the trajectory & Guess the bonds if not present
        stamp0 = time.time()
        universe = mda.Universe(self.topo, self.traj)
        try:
            any_bond = universe.bonds[0]
        except:
            logger.warning(f'The passed topology does not contain bonds. '
                           f'MDTraj will guess them automatically.')
            universe = mda.Universe(self.topo, self.traj, guess_bonds=True)
            any_bond = universe.bonds[0]
        if any_bond is None:
            raise ValueError(
                "No bonds found in topology and MDAnalysis could not guess them.")

        # Remove the hydrogen-hydrogen bonds if any
        stamp1 = time.time()
        are_hh = any_hh_bonds(universe)
        universe.delete_bonds(get_hh_bonds(universe))
        del_time = time.time() - stamp1
        if are_hh:
            logger.warning(f"This universe contains H-H bonds.\n"
                           f"Time to remove H-H bonds: {del_time:.2f} s")

        # Output load time
        loading = time.time() - stamp0
        logger.info(f"Time to load the trajectory: {loading:.2f} s")
        return universe

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
                selected_dh = tri_res.atoms.indices[match_dh][where_dh]

                for i, x in enumerate(selected_dh):
                    if i % 2 == 0:
                        hx_D.append(x)
                    else:
                        hx_H.append(x)

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
                        if i % 2 == 0:
                            hx_D.append(x[0])
                        else:
                            hx_H.append(x[1])
            match_a = [x for x in mono_mol.GetSubstructMatches(query_a)]
            if match_a:
                for similar in self.uniq_disconnected[case]:
                    mono_res = self.universe.residues[similar]
                    selected_a = mono_res.atoms.indices[match_a]
                    hx_A.extend(selected_a[0])

        # Filter the indices
        hx_D_raw = npi.indices(self.sel_idx, np.asarray(hx_D), missing=-1)
        hx_H_raw = npi.indices(self.sel_idx, np.asarray(hx_H), missing=-1)
        hx_A_raw = npi.indices(self.sel_idx, np.asarray(hx_A), missing=-1)

        hx_D = hx_D_raw[hx_D_raw != -1].astype(np.int32)
        hx_H = hx_H_raw[hx_H_raw != -1].astype(np.int32)
        hx_A = np.unique(hx_A_raw[hx_A_raw != -1]).astype(np.int32)
        assert len(hx_D) == len(hx_H), "Donors and Hydrogens do not match"
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
        return padded_rings

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
        return radii.astype(np.float32)

    def get_interactions(self):
        """
        Get the possible interactions between the selections in the trajectory

        Returns:

        """

        len_s1, len_s2 = len(self.sel1_idx), len(self.sel2_idx)
        len_hp = len(self.hydroph)
        len_an, len_ca = len(self.anions), len(self.cations)
        len_ma, len_md = len(self.metal_acc), len(self.metal_don)
        len_hbd, len_hba, len_hbh = len(self.hb_D), len(self.hb_A), len(
            self.hb_H)
        len_xbd, len_xba, len_xbh = len(self.xb_D), len(self.xb_A), len(
            self.xb_H)
        len_rings = len(self.rings)

        is_possible = {
            'CloseContacts': (len_s1 > 0) and (len_s2 > 0),
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
        }

        raw = self.raw_inters
        for x in raw:
            if x == 'all':
                requested = list(cf.parse_cutoffs().keys())
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
                f"selections for the given system:\n"
                f"{to_skip} ")
        return to_compute

# =============================================================================
#
# =============================================================================
# import intermap.config as conf
# from argparse import Namespace
#
# config_path = '/home/gonzalezroy/RoyHub/intermap/tests/imaps/imap1.cfg'
# config = conf.ConfigManager(config_path, conf.allowed_parameters)
# args = Namespace(**config.config_args)
# self = IndexManager(args.topology, args.trajectory, args.selection_1,
#                     args.selection_2, args.interactions)
