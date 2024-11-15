# Created by gonzalezroy at 11/15/24
import numpy as np

heavies = ['S', 'N', 'O', 'P']


class IndexManager:
    """
    Class to manage the indices of the selections in a trajectory.
    """

    def __init__(self, sel1, sel2, master_traj):
        """
        Initialize the IndexManager class.

        Args:
            sel1: MDtraj selection string # 1
            sel2: MDtraj selection string # 2
            master_traj: MDtraj trajectory object (one frame only)
        """
        # Parse arguments
        self.s1 = sel1
        self.s2 = sel2
        self.traj = master_traj

        # Get the topology and bonds
        self.topo, bonds = self.traj.topology.to_dataframe()
        self.bonds = bonds.astype(int)

        # Get the indices of the selections
        self.s1_idx, self.s2_idx = self.get_selections_indices()

        # [hbonds]: get the indices of (D)onors, (H)ydrogens and (A)cceptor
        self.donors, self.hydros, self.heavies = self.get_all_dha_idx()
        self.s1_donors, self.s1_hydros, self.s1_acc = self.get_sel_dha_idx(
            self.s1_idx)
        self.s2_donors, self.s2_hydros, self.s2_acc = self.get_sel_dha_idx(
            self.s2_idx)

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
