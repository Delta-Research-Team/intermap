# Created by rglez at 4/20/25
import configparser
from collections import defaultdict
from os.path import dirname

import MDAnalysis as mda
import numpy as np
import pandas as pd
import rgpack.generals as gnl
from bitarray import bitarray as ba, util as bu

from intermap.shiny.app.css import all_interactions_colors

pd.options.mode.chained_assignment = None

errors = {
    'noTopo':
        'Topology not found. Please manually copy the topology file used to run '
        'InterMap to the same directory as the csv file being loaded.'}


def transpose(df):
    """
    Transpose the DataFrame; sel1 becomes sel2 and sel2 becomes sel1

    Args:
        df (pd.DataFrame): DataFrame to transposes
    """
    df.rename(
        columns={'sel1': 'sel2', 'sel2': 'sel1',
                 'note1': 'note2', 'note2': 'note1',
                 'idx1': 'idx2', 'idx2': 'idx1',
                 }, inplace=True)
    return df


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


def sortby(df, choice):
    """
    Sort the DataFrame by the selected column

    Args:
        df (pd.DataFrame): DataFrame to sort
        choice (str): column to sort by
    """
    choices = ['name', 'number', 'annotation']
    if choice not in choices:
        raise ValueError(
            f'Invalid choice: {choice}. Available choices: {choices}')

    if choice == 'name':
        sorted_df = df.sort_values(by=['sel1', 'idx1', 'note1',
                                       'sel2', 'idx2', 'note2'])
    elif choice == 'number':
        sorted_df = df.sort_values(by=['idx1', 'note1', 'sel1',
                                       'idx2', 'note2', 'sel2'])
    else:
        sorted_df = df.sort_values(by=['note1', 'idx1', 'sel1',
                                       'note2', 'idx2', 'sel2'])
    return sorted_df


def parse_pickle(pickle):
    """
    Parse the InterMap output pickle file

    Args:
        pickle (str): path to the InterMap output pickle file

    Returns:
        bit_dict (dict): dictionary containing the interactions
    """

    # Create the DataFrame from the pickle file
    pickle_path = gnl.check_path(pickle)
    bit_dict = gnl.unpickle_from_file(pickle_path)
    df = pd.DataFrame(bit_dict.keys(), columns=[
        'sel1', 'note1', 'sel2', 'note2', 'water', 'interaction_name'])

    # Extract the indices from the selections
    raw = df['sel1'][0]
    resolution = 'residue' if len(raw.split('_')) == 3 else 'atom'
    N = 2 if resolution == 'residue' else 4
    idx_1 = df['sel1'].str.split('_', expand=True)[N].astype(int)
    idx_2 = df['sel2'].str.split('_', expand=True)[N].astype(int)
    df['idx1'] = idx_1
    df['idx2'] = idx_2

    # Add the timeseries and the prevalence column
    df['timeseries'] = bit_dict.values()
    df['timeseries'] = df['timeseries'].apply(bu.sc_decode)
    df['prevalence'] = df['timeseries'].apply(
        lambda x: x.count() / len(x) * 100)
    df['timeseries'] = df['timeseries'].apply(ba.to01)

    return df, resolution



def parse_cfg(cfg):
    """
    Parse the InterMap configuration file

    Args:
        cfg (str): path to the InterMap configuration file

    Returns:
        config_obj (configparser.ConfigParser): ConfigParser object with the
            configuration settings
    """
    config_obj = configparser.ConfigParser(allow_no_value=True,
                                           inline_comment_prefixes='#')
    config_obj.optionxform = str
    config_obj.read(cfg)
    config_dict = {section: dict(config_obj.items(section))
                   for section in config_obj.sections()}

    axisx = config_dict['interactions']['selection_1']
    axisy = config_dict['interactions']['selection_2']

    return config_dict, axisx, axisy


class CSVFilter:
    """
    A class to filter the InterMap CSV file before loading it with Shiny
    """
    mda_topols = (
        'psf', 'pdb', 'ent', 'pqr', 'pdbqt', 'gro', 'top', 'prmtop', 'parm7',
        'dms', 'tpr', 'itp', 'mol2', 'data', 'lammpsdump', 'xyz', 'txyz',
        'arc', 'gms', 'log', 'config', 'history', 'xml', 'gsd', 'mmtf', 'in')

    def __init__(self, pickle, cfg):
        """
        Initialize the CSVFilter class

        Args:
            pickle (str): path to the InterMap output pickle file
            cfg (str): path to the InterMap configuration file
        """

        # Parse the objects from arguments
        self.root_dir = dirname(pickle)
        self.master, self.resolution = parse_pickle(pickle)
        self.cfg, self.axisx, self.axisy = parse_cfg(cfg)
        self.prevalences = self.master['prevalence']

        # Map the indices
        self.universe = self.load_universe()
        self.at_idx = self.universe.atoms.indices
        self.res_idx = self.universe.atoms.resindices
        self.at2res = self.at2res()
        self.res2at = self.res2at()
        self.notes2df = self.notes2df()
        self.inters2df = self.inters2df()

    def load_universe(self):
        """
        Load the MDAnalysis Universe from the configuration file

        Returns:
            universe (MDAnalysis.Universe): Universe object containing the
             topology and trajectory data
        """
        topo = self.cfg['topo-traj']['topology']
        traj = self.cfg['topo-traj']['trajectory']
        universe = mda.Universe(topo, traj)
        return universe

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
#pickle = '/media/rglez/Roy2TB/Dropbox/RoyData/intermap/tutorial-mayank/outputs/ligs-channel_InterMap.pickle'
#cfg = '/media/rglez/Roy2TB/Dropbox/RoyData/intermap/tutorial-mayank/outputs/ligs-channel_InterMap.cfg'
#self = CSVFilter(pickle, cfg)
#self.universe
#self.master.head()
