# Created by rglez at 4/20/25
from collections import defaultdict
from os.path import dirname

import MDAnalysis as mda
import numpy as np
import pandas as pd
import rgpack.generals as gnl
from bitarray import bitarray as ba

pd.options.mode.chained_assignment = None

errors = {
    'noTopo':
        'Topology not found. Please manually copy the topology file used to run '
        'InterMap to the same directory as the csv file being loaded.'}


def compress_wb(df):
    """
    Compress the timeseries of the WaterBridge interactions.
    This is a destructive operation on the original df

    Args:
        df (pd.DataFrame): DataFrame containing the timeseries data

    Returns:
        df (pd.DataFrame): DataFrame with the timeseries data compressed
    """
    # Convert timeseries to bitarrays
    df['timeseries'] = df['timeseries'].apply(ba)

    # Detect wbs between two residues with exchanging waters
    wbs = df[df.interaction_name == 'WaterBridge'].groupby(
        ['sel1', 'sel2']).indices
    wbs = {k: v for k, v in wbs.items() if len(v) > 1}

    # Merge the timeseries of the waters
    to_del = []
    for case in wbs:
        indices = wbs[case]
        ref_idx = indices[0]
        bit_ref = df['timeseries'].loc[ref_idx].copy()
        for idx in indices:
            bit_ref |= df.loc[idx, 'timeseries']
            if idx != ref_idx:
                to_del.append(idx)
        df['timeseries'].iloc[ref_idx] = bit_ref

    # Drop the merged timeseries
    df.drop(to_del, inplace=True)

    # Convert bitarrays to numpy arrays for future plotting
    df['timeseries'] = df['timeseries'].apply(lambda x: np.asarray(x.tolist()))
    return df


class CSVFilter:
    """
    A class to filter the InterMap CSV file before loading it with Shiny
    """
    mda_topols = (
        'psf', 'pdb', 'ent', 'pqr', 'pdbqt', 'gro', 'top', 'prmtop', 'parm7',
        'dms', 'tpr', 'itp', 'mol2', 'data', 'lammpsdump', 'xyz', 'txyz',
        'arc', 'gms', 'log', 'config', 'history', 'xml', 'gsd', 'mmtf', 'in')

    def __init__(self, csv, topo):
        """
        Initialize the CSVFilter class

        Args:
            csv (str): path to the InterMap CSV file
            topo (str): path to the topology file
        """
        # Parse arguments
        self.csv_path = gnl.check_path(csv)
        if self.csv_path.split('.')[-1] != 'csv':
            raise ValueError('CSV file must have a .csv extension')

        # Validate topology file
        self.topo_path = gnl.check_path(topo)
        if self.topo_path.split('.')[-1] not in self.mda_topols:
            raise ValueError(
                'Topology file must be in a format compatible with MDAnalysis.'
                f' The following are supported: {self.mda_topols}')

        self.root_dir = dirname(self.csv_path)

        # Load the main objects
        try:
            self.universe = mda.Universe(self.topo_path)
        except Exception as e:
            raise ValueError(f"Error loading topology file: {str(e)}")

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

    def parse_csv(self):
        """
        Parse the InterMap CSV file

        Returns:
            csv (pd.DataFrame): DataFrame containing the InterMap data
        """
        df = pd.read_csv(self.csv_path, header=1, na_values=('', ' '),
                         dtype=str)
        raw = df['sel1'][0]
        resolution = 'residue' if len(raw.split('_')) == 3 else 'atom'

        N = 2 if resolution == 'residue' else 4
        idx_1 = df['sel1'].str.split('_', expand=True)[N].astype(int)
        idx_2 = df['sel2'].str.split('_', expand=True)[N].astype(int)
        df['idx1'] = idx_1
        df['idx2'] = idx_2
        df = compress_wb(df)
        return df, resolution

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
full_csv = '/home/gonzalezroy/RoyHub/intermap/tests/imaps/nuc-wat/8oxoGA2_1_InterMap_full.csv'
topology = '/home/gonzalezroy/RoyHub/intermap/tests/imaps/nuc-wat/ncp-OCS-nosalt.prmtop'
self = CSVFilter(full_csv, topology)
