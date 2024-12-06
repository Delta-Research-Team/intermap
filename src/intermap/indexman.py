# Created by gonzalezroy at 11/15/24
from collections import defaultdict

import networkx as nx
import numpy as np

# Constants and global variables
valid_interactions = ['hbonds', 'close_contacts']
heavies = ['S', 'N', 'O', 'P']


def parse_interactions(interactions):
    """
    Parse the interactions declared by the user.

    Args:
        interactions: str, 'all' or list of interactions to calculate

    Returns:
        interactions: list of interactions to calculate
    """

    if interactions == 'all':
        return valid_interactions
    elif isinstance(interactions, list):
        for interaction in interactions:
            if interaction not in valid_interactions:
                raise ValueError(f"Invalid interaction: {interaction}")
        return interactions
    else:
        raise ValueError(f"Invalid interactions declared: {interactions}")


class IndexManager:
    """
    Class to manage the indices of the selections in a trajectory.
    """

    def __init__(self, sel1, sel2, master_traj, interactions='all'):
        """
        Initialize the IndexManager class.

        Args:
            sel1: MDtraj selection string # 1
            sel2: MDtraj selection string # 2
            master_traj: MDtraj trajectory object (one frame only)
            interactions: str, 'all' or list of interactions to calculate
        """
        # Parse arguments
        self.s1 = sel1
        self.s2 = sel2
        self.traj = master_traj

        # Parse interactions
        self.inters = parse_interactions(interactions)

        # Get the topology, bonds and labels
        self.topo, bonds = self.traj.topology.to_dataframe()
        self.bonds = bonds.astype(int)
        self.coords = self.traj.xyz[0]
        self.labels = self.get_labels()

        # Get the indices of the selections
        self.s1_idx, self.s2_idx = self.get_selections_indices()

        # Get the indices of interactions selections
        self.donors, self.hydros, self.heavies = self.get_all_dha_idx()
        self.indices = self.get_interaction_indices()

    def get_labels(self):
        """
        Get the labels of the atoms in the topology.

        Returns:
            labels (list): labels of the atoms in the topology
        """
        df = self.topo
        labels = [f'{x[0]}-{x[1]}-{x[2]}'
                  for x in zip(df.resName, df.resSeq, df.name)]
        return np.asarray(labels)

    def get_selections_indices(self):
        """
        Get the indices of the selections in the trajectory.

        Returns:
            s1_idx: array with the indices of the atoms in selection 1
            s2_idx: array with the indices of the atoms in selection 2
        """
        s1_idx = self.traj.top.select(self.s1)
        s2_idx = self.traj.top.select(self.s2)

        if len(s1_idx) == 0:
            raise ValueError("No atoms found for selection 1")
        if len(s2_idx) == 0:
            raise ValueError("No atoms found for selection 2")
        return s1_idx, s2_idx

    def get_all_dha_idx(self):
        """
        Get the indices of all (D)onors, (H)ydrogens and (A)cceptor.

        Returns:
            donors: array with the indices of the donor atoms
            hydros: array with the indices of the hydrogen atoms
            all_heavies: array with the indices of the heavy atoms

        """
        topo = self.topo
        bonds = self.bonds
        all_hydro = topo[topo.element == "H"].index
        all_heavies = np.where(topo.element.isin(heavies))[0]

        ats1 = bonds[:, 0]
        ats2 = bonds[:, 1]

        at1_is_h = np.isin(ats1, all_hydro)
        at2_is_heavy = np.isin(ats2, all_heavies)
        at2_is_h = np.isin(ats2, all_hydro)
        at1_is_heavy = np.isin(ats1, all_heavies)

        hydros = np.concatenate([ats1[np.bitwise_and(at1_is_h, at2_is_heavy)],
                                 ats2[np.bitwise_and(at2_is_h, at1_is_heavy)]])
        donors = np.concatenate([ats2[np.bitwise_and(at1_is_h, at2_is_heavy)],
                                 ats1[np.bitwise_and(at2_is_h, at1_is_heavy)]])

        return donors, hydros, all_heavies

    def get_sel_dha_idx(self, sel_idx):
        """
        Get the indices of the (D)onors, (H)ydrogens and (A)cceptor in a selection

        Args:
            sel_idx: array with the indices of the atoms in the selection

        Returns:
            sel_donors: array with the indices of the donor atoms in the selection
            sel_hydros: array with the indices of the hydrogen atoms in the selection
            sel_heavies: array with the indices of the heavy atoms in the selection
        """
        donors_hydros = dict(zip(self.donors, self.hydros))
        sel_heavies = np.intersect1d(sel_idx, self.heavies)
        sel_donors = np.intersect1d(sel_idx, self.donors)
        sel_hydros = np.asarray([donors_hydros[donor] for donor in sel_donors])
        return sel_donors, sel_hydros, sel_heavies

    def get_interaction_indices(self):
        """'
        Get the indices of the atoms involved in the interactions.

        Returns:
            selection_indices (dict): dictionary with the indices of the atoms
                involved in the interactions
        """
        selection_indices = {}

        if ('hbonds' in self.inters) or ('all' in self.inters):
            s1_donors, s1_hydros, s1_acc = self.get_sel_dha_idx(self.s1_idx)
            s2_donors, s2_hydros, s2_acc = self.get_sel_dha_idx(self.s2_idx)
            selection_indices['hbonds'] = {'s1_donors': s1_donors,
                                           's1_hydros': s1_hydros,
                                           's1_acc': s1_acc,
                                           's2_donors': s2_donors,
                                           's2_hydros': s2_hydros,
                                           's2_acc': s2_acc}

        if ('close_contacts' in self.inters) or ('all' in self.inters):
            selection_indices['close_contacts'] = {'s1_idx': self.s1_idx,
                                                   's2_idx': self.s2_idx}

        if ('aromatic' in self.inters) or ('all' in self.inters):
            selection_indices['aromatic'] = {'s1_centroids': self.s1_idx,
                                             's2_centroids': self.s2_idx}

        return selection_indices

    # =========================================================================
    # Aromatics
    # =========================================================================
    def get_sel_aro_idx(self, sel_idx, tolerance=25):
        coordinates = self.coords
        at1_in_sel = np.isin(self.bonds[:, 0], sel_idx)
        at2_in_sel = np.isin(self.bonds[:, 1], sel_idx)
        bonds_in_sel = self.bonds[np.bitwise_or(at1_in_sel, at2_in_sel)]

        # Create a dictionary of bonds for each atom
        topo_bonds = defaultdict(list)
        for bond in bonds_in_sel:
            topo_bonds[bond[0]].append(bond[1])

        # Convert topo_bonds to a NetworkX graph
        G = nx.Graph(topo_bonds)
        all_cycles = list(nx.simple_cycles(G))

        # Find all cycles in the topo_bonds graph
        interval = range(5, 7)
        cycles = [x for x in all_cycles if len(x) in interval]

        for c in cycles:
            if 18947 in c:
                print(c)

        cycles = [[18945, 18947, 18935, 18936, 18938, 18943]]

        # Check if the cycles are planar
        planar_cycles = []
        for cycle in cycles:
            num_atoms = len(cycle)
            passing = 0
            for i in range(num_atoms):
                c1 = coordinates[cycle[i]]
                c2 = coordinates[cycle[(i + 1) % num_atoms]]
                c3 = coordinates[cycle[(i + 2) % num_atoms]]
                c4 = coordinates[cycle[(i + 3) % num_atoms]]
                dihedral_angle = self._calculate_dihedral_angle(c1, c2, c3, c4)
                print(i, dihedral_angle)
                angle1 = abs(dihedral_angle)
                angle2 = abs(dihedral_angle - 180)
                if (angle1 > tolerance) and (angle2 > tolerance):
                    continue
                else:
                    passing += 1
            if passing == num_atoms:
                planar_cycles.append(cycle)
        return planar_cycles

    def _calculate_dihedral_angle(self, coords1, coords2, coords3, coords4):
        """
        Calculate the dihedral angle between four points in space.

        Args:
            coords1, coords2, coords3, coords4: Coordinates of the four points.

        Returns:
            Dihedral angle in degrees.
        """
        b1 = coords2 - coords1
        b2 = coords3 - coords2
        b3 = coords4 - coords3

        b1xb2 = np.cross(b1, b2)
        b2xb3 = np.cross(b2, b3)

        b1xb2_x_b2xb3 = np.cross(b1xb2, b2xb3)

        y = np.dot(b1xb2_x_b2xb3, b2) * (1.0 / np.linalg.norm(b2))
        x = np.dot(b1xb2, b2xb3)

        return np.degrees(np.arctan2(y, x))


# =============================================================================
# Debugging area
# =============================================================================
# from argparse import Namespace
# import mdtraj as md
#
# topo = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A1/8oxoGA1_1_hmr.prmtop'
# traj = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A1/8oxoGA1_1_sk100.nc'
# sel1 = "(resname =~ '(5|3)?D([ATGC])|(8OG){1}(3|5)?$')"
# sel2 = "protein"
# nprocs = 8
# chunk_size = 150
# args = Namespace(topo=topo, traj=traj, sel1=sel1, sel2=sel2, nprocs=nprocs,
#                  chunk_size=chunk_size)
#
# master_traj = md.load_frame(args.traj, top=args.topo, index=0)
#
# self = IndexManager(sel1, sel2, master_traj)
