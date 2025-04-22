# Created by rglez at 4/20/25
import os
from collections import defaultdict
from os.path import dirname, join

import MDAnalysis as mda
import numpy as np
import pandas as pd
import rgpack.generals as gnl

errors = {
    'noTopo':
        'Topology not found. Please manually copy the topology file used to run '
        'InterMap to the same directory as the csv file being loaded.'}


class CSVFilter:
    """
    A class to filter the InterMap CSV file before loading it with Shiny
    """

    def __init__(self, csv):
        """
        Initialize the CSVFilter class

        Args:
            csv (str): path to the InterMap CSV file
        """
        # Parse arguments
        self.csv_path = gnl.check_path(csv)
        self.root_dir = dirname(self.csv_path)
        self.topo_path = self.check_topology()

        # Load the main objects
        self.universe = mda.Universe(self.topo_path)
        self.master, self.resolution = self.parse_csv()
        self.long = True if 'timeseries' in self.master.columns else False
        self.prevalences = self.master['prevalence'].astype(float)

        # Map the indices
        self.at_idx = self.universe.atoms.indices
        self.res_idx = self.universe.atoms.resindices
        self.at2res = self.at2res()
        self.res2at = self.res2at()
        self.notes2df = self.notes2df()
        self.inters2df = self.inters2df()

    """
    def check_topology(self):
        with open(self.csv_path, 'rt') as f:
            topo = next(f).split('#')[-1].strip()
            topo_path = join(self.root_dir, topo)
            if not os.path.isfile(topo_path):
                raise FileNotFoundError(errors['noTopo'].format(self.csv_path))
            return topo_path
    """


    def check_topology(self):
        TOPO_DIR = "/home/fajardo01/03_Fajardo_Hub/02_InterMap/visualizations/data/last_version/DrHU-Tails-8k/DrHU-Tails-8k/"

        with open(self.csv_path, 'rt') as f:
            topo = next(f).split('#')[-1].strip()
            topo_path = os.path.join(TOPO_DIR, topo)
            print(f"Looking for topology file at: {topo_path}")
            if not os.path.isfile(topo_path):
                raise FileNotFoundError(errors['noTopo'].format(self.csv_path))
            return topo_path


    def parse_csv(self):
        """
        Parse the InterMap CSV file

        Returns:
            csv (pd.DataFrame): DataFrame containing the InterMap data
        """
        csv = pd.read_csv(self.csv_path, header=1, na_values=('', ' '))
        raw = csv['sel1'][0]
        resolution = 'residue' if len(raw.split('_')) == 3 else 'atom'

        N = 2 if resolution == 'residue' else 4
        idx_1 = csv['sel1'].str.split('_', expand=True)[N].astype(int)
        idx_2 = csv['sel2'].str.split('_', expand=True)[N].astype(int)
        csv['idx1'] = idx_1
        csv['idx2'] = idx_2

        return csv, resolution

    # =========================================================================
    # Static mappings
    # =========================================================================
    def at2res(self):
        """
        Map the atom indices to residue indices

        Returns:
            at2res (dict): dictionary mapping atom indices to residue indices
        """
        at2res = {i: self.res_idx[i] for i in range(len(self.at_idx))}
        return at2res

    def res2at(self):
        """
        Map the residue indices to atom indices

        Returns:
            res2at (dict): dictionary mapping residue indices to atom indices
        """
        res2at = defaultdict(set)
        for i, res in enumerate(self.res_idx):
            res2at[res].add(self.at_idx[i])
        return res2at

    def notes2df(self):
        """
        Map the annotations to the DataFrame indices

        Returns:
            notes2df (dict): dictionary mapping annotations to DataFrame
             indices
        """
        notes2df = defaultdict(set)
        for x in ['note1', 'note2']:
            sub_df = self.master[x][self.master[x].isnull() == False]
            if not sub_df.empty:
                notes = sub_df.str.strip()
                for i, y in enumerate(notes):
                    if not pd.isna(y):
                        notes2df[y].add(i)
        return notes2df

    def inters2df(self):
        """
        Map the interaction names to the DataFrame indices

        Returns:
            inters2df (dict): dictionary mapping interaction names to DataFrame
             indices
        """
        inters2df = defaultdict(set)
        inters = self.master['interaction_name'].str.strip()
        for i, y in enumerate(inters):
            inters2df[y].add(i)
        return inters2df

    # =========================================================================
    # Dynamic filters
    # =========================================================================

    def by_mda(self, sele):
        """
        Filter the InterMap CSV file based on the MDAnalysis selection string

        Args:
            sele (str): selection string compliant with MDAnalysis syntax
        """
        idx = set()
        if self.resolution == 'atom':
            indices = set(self.universe.select_atoms(sele).indices)
        else:
            indices = set(self.universe.select_atoms(sele).resindices)

        for x in ['idx1', 'idx2']:
            where = np.where(self.master[x].isin(indices))[0]
            idx.update(where)
        status = 0 if len(indices) > 0 else -1
        return idx, status

    def by_prevalence(self, preval):
        """
        Map the prevalence to the DataFrame indices

        Args:
            preval (float): prevalence value to filter the DataFrame

        Returns:
            idx (np.ndarray): indices of the DataFrame with prevalence >= preval
            status (int): status of the filter. 0 if successful, -1 if no
            indices found
        """
        idx = np.asarray(self.prevalences[self.prevalences >= preval].index)
        status = 0 if len(idx) > 0 else -1
        return idx, status

    def by_inters(self, inters):
        """
        Map the interaction names to the DataFrame indices

        Args:
            inters (tuple): tuple of interaction names to filter the DataFrame

        Returns:

        """
        idx = set()
        if inters == 'all':
            inters = set(self.inters2df.keys())
        else:
            inters = inters

        for inter in inters:
            if inter in self.inters2df:
                idx.update(self.inters2df[inter])
        status = 0 if len(idx) > 0 else -1
        return idx, status

    def by_notes(self, notes):
        """
        Map the annotations to the DataFrame indices

        Args:
            notes (tuple): tuple of annotations to filter the DataFrame

        Returns:
            idx (np.ndarray): indices of the DataFrame with prevalence >= preval
            status (int): status of the filter. 0 if successful, -1 if no
            indices found
        """
        idx = set()
        if notes == 'all':
            notes = set(self.notes2df.keys())
        else:
            notes = notes

        for note in notes:
            if note in self.notes2df:
                idx.update(self.notes2df[note])
        status = 0 if len(idx) > 0 else -1
        return idx, status


# =============================================================================
#
# =============================================================================
full = '/home/fajardo01/03_Fajardo_Hub/02_InterMap/visualizations/data/last_version/DrHU-Tails-8k/DrHU-Tails-8k/dna-prot_InterMap_full.csv'
short = '/home/fajardo01/03_Fajardo_Hub/02_InterMap/visualizations/data/last_version/DrHU-Tails-8k/DrHU-Tails-8k/dna-prot_InterMap_short.csv'

self = CSVFilter(full)

#mda_sele, mda_status = self.by_mda('all')
#prevalence, preval_status = self.by_prevalence(95)
#interactions, inters_status = self.by_inters('all')
#annotations, notes_status = self.by_notes('all')

#df_idx = set.intersection(mda_sele, prevalence, interactions, annotations)
#df_status = 0 if len(df_idx) > 0 else -1
#print(df_status)
#self.master.iloc[260]
#self.master.iloc[list(mda_sele)].columns
#self.master.iloc[list(mda_sele)]['sel1']
#self.master.iloc[list(mda_sele)]['sel2']
#self.master.iloc[prevalence]['prevalence']
