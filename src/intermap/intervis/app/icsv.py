# Created by rglez at 4/20/25
import configparser
from collections import defaultdict
from os.path import dirname

import MDAnalysis as mda
import numpy as np
import pandas as pd
import rgpack.generals as gnl
from bitarray import bitarray as ba, util as bu

from intermap.intervis.app.css import all_interactions_colors

pd.options.mode.chained_assignment = None

errors = {
    'noTopo':
        'Topology not found. Please manually copy the topology file used to run '
        'InterMap to the same directory as the csv file being loaded.'}

interaction_priority = {
    'Anionic': 1, 'Cationic': 2, 'HBDonor': 3, 'HBAcceptor': 4,
    'MetalDonor': 5, 'MetalAcceptor': 6, 'PiCation': 7, 'CationPi': 8,
    'PiAnion': 9, 'PiStacking': 10, 'FaceToFace': 11, 'EdgeToFace': 12,
    'XBDonor': 13, 'XBAcceptor': 14, 'Hydrophobic': 15, 'VdWContact': 16,
    'CloseContact': 17, 'WaterBridge': 18
}


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
    #df['timeseries'] = df['timeseries'].apply(bu.sc_decode)
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


# =========================================================================
# Funtions to process plots
# =========================================================================


def process_heatmap_data(df):
    interaction_priority = {
        'Anionic': 1, 'Cationic': 2, 'HBDonor': 3, 'HBAcceptor': 4,
        'MetalDonor': 5, 'MetalAcceptor': 6, 'PiCation': 7, 'CationPi': 8,
        'PiAnion': 9, 'PiStacking': 10, 'FaceToFace': 11, 'EdgeToFace': 12,
        'XBDonor': 13, 'XBAcceptor': 14, 'Hydrophobic': 15, 'VdWContact': 16,
        'CloseContact': 17, 'WaterBridge': 18
    }

    df['priority'] = df['interaction_name'].map(interaction_priority)

    priority_df = (df.sort_values(['sel1', 'sel2', 'priority', 'prevalence'],
                                  ascending=[True, True, True, False])
                   .groupby(['sel1', 'sel2']).first().reset_index())

    pivot_interaction = pd.pivot_table(priority_df, values='interaction_name',
                                       index='sel2', columns='sel1',
                                       aggfunc='first', fill_value='')

    pivot_prevalence = pd.pivot_table(priority_df, values='prevalence',
                                      index='sel2', columns='sel1',
                                      aggfunc='first', fill_value="")

    pivot_prevalence_rounded = pivot_prevalence.round(1).astype(str)

    pivot_note1 = pd.pivot_table(priority_df, values='note1',
                                 index='sel2', columns='sel1',
                                 aggfunc='first', fill_value="")

    pivot_note2 = pd.pivot_table(priority_df, values='note2',
                                 index='sel2', columns='sel1',
                                 aggfunc='first', fill_value="")

    present_interactions = sorted(priority_df['interaction_name'].unique(),
                                  key=lambda x: interaction_priority[x])

    return {
        'pivot_interaction': pivot_interaction,
        'pivot_prevalence': pivot_prevalence,
        'pivot_prevalence_rounded': pivot_prevalence_rounded,
        'pivot_note1': pivot_note1,
        'pivot_note2': pivot_note2,
        'present_interactions': present_interactions
    }


def process_prevalence_data(df, selection_column, batch_size=250):
    """Process data for interaction plots.

    Args:
        df: DataFrame with interaction data
        selection_column: Column to use for selection ('sel1' for ligand, 'sel2' for receptor)
        batch_size: Size of batches for processing

    Returns:
        List of dictionaries with processed data for plotting
    """
    df['annotation'] = df['note1'].fillna('') + ' ' + df['note2'].fillna('')
    df['annotation'] = df['annotation'].str.strip()

    sort_idx = np.lexsort((-df['prevalence'].to_numpy(),
                           df[selection_column].to_numpy()))
    sorted_df = df.iloc[sort_idx]

    batched_data = []
    legend_entries = set()
    unique_selections = sorted_df[selection_column].unique()

    for i in range(0, len(unique_selections), batch_size):
        batch_selections = unique_selections[i:i + batch_size]
        batch_data = sorted_df[
            sorted_df[selection_column].isin(batch_selections)]

        for _, row in batch_data.iterrows():
            show_legend = row['interaction_name'] not in legend_entries
            if show_legend:
                legend_entries.add(row['interaction_name'])

            batched_data.append({
                'sel1': row['sel1'],
                'sel2': row['sel2'],
                'prevalence': row['prevalence'],
                'interaction_name': row['interaction_name'],
                'show_legend': show_legend,
                'color': all_interactions_colors[row['interaction_name']],
                'annotation': row['annotation']
            })

    return batched_data


def process_prevalence_data2(df, selection_column, sort_by='note'):
    """Process data for interaction plots.

    Args:
        df: DataFrame with interaction data
        selection_column: Column to use for selection ('sel1' for ligand, 'sel2' for receptor)
        sort_by: Additional sorting column ('resname', 'resnum', 'idx', 'note')

    Returns:
        List of dictionaries with processed data for plotting
    """
    df['annotation'] = df['note1'].fillna('') + ' ' + df['note2'].fillna('')
    df['annotation'] = df['annotation'].str.strip()
    # Assign colors based on interaction names
    sorted_df = df.copy()
    inter_name = df['interaction_name']
    inter_colors = all_interactions_colors
    sorted_df['color'] = inter_name.map(inter_colors)

    # Show legend for the first occurrence of each interaction
    _, indices = np.unique(inter_name.to_numpy(), return_index=True)
    show_legend = np.zeros(len(sorted_df), dtype=bool)
    show_legend[indices] = True
    sorted_df['show_legend'] = show_legend

    # Sorting
    which = selection_column[-1]
    resname = f'resname{which}'
    resnum = f'resnum{which}'
    idx = f'idx{which}'
    note = f'note{which}'

    # Sort by the selected column
    if sort_by == 'resname':
        groups = sorted_df.groupby([resname, resnum, idx, note]).indices
    elif sort_by == 'resnum':
        groups = sorted_df.groupby([resnum, resname, idx, note]).indices
    elif sort_by == 'idx':
        groups = sorted_df.groupby([idx, resname, resnum, note]).indices
    elif sort_by == 'note':
        groups = sorted_df.groupby([note, resname, resnum, idx]).indices
    else:
        raise ValueError(
            f'Invalid sort_by value: {sort_by}. Available options: '
            'resname, resnum, idx, note')

    # Sort by prevalence
    for group in groups:
        indices = groups[group]
        sorted_indices = (sorted_df.iloc[indices].sort_values(
            by='prevalence', ascending=False).index)
        groups[group] = sorted_indices
    concat_indices = np.concatenate(list(groups.values()))
    sorted_df = sorted_df.iloc[concat_indices]

    # Create batched data
    batched_df = sorted_df[['sel1', 'sel2', 'prevalence', 'interaction_name',
                            'show_legend', 'color', 'annotation']]
    batched_data = batched_df.to_dict(orient='records')
    return batched_data


def process_time_series_data(df):
    """Process data for time series plot."""
    interaction_abbreviations = {
        'HBDonor': 'HBD', 'HBAcceptor': 'HBA', 'Cationic': 'Cat',
        'Anionic': 'Ani', 'WaterBridge': 'WB', 'PiStacking': 'πS',
        'PiCation': 'πC', 'CationPi': 'Cπ', 'PiAnion': 'πA',
        'FaceToFace': 'F2F',
        'EdgeToFace': 'E2F', 'MetalDonor': 'MD', 'MetalAcceptor': 'MA',
        'VdWContact': 'VdW', 'CloseContact': 'CC', 'Hydrophobic': 'Hyd',
        'XBAcceptor': 'XBA', 'XBDonor': 'XBD'
    }

    df = df.copy()
    df['selection_pair'] = (df['sel1'] + ' - ' + df['sel2'] + ' (' +
                            df['interaction_name'].map(
                                interaction_abbreviations) + ')')

    frame_interactions = []
    for _, row in df.iterrows():
        timeseries = np.array(list(row['timeseries']), dtype=int)
        frames_with_interaction = np.where(timeseries == 1)[0]
        for frame in frames_with_interaction:
            frame_interactions.append({
                'selection_pair': row['selection_pair'],
                'frame': frame,
                'prevalence': row['prevalence'],
                'interaction_name': row['interaction_name']
            })

    scatter_df = pd.DataFrame(frame_interactions)
    prevalence_data = df.groupby('selection_pair')['prevalence'].mean()

    interaction_colors = []
    for pair in prevalence_data.index:
        interaction_abbrev = pair.split('(')[1].rstrip(')')
        for full_name, abbrev in interaction_abbreviations.items():
            if abbrev == interaction_abbrev:
                interaction_colors.append(all_interactions_colors[full_name])
                break

    return {
        'scatter_df': scatter_df,
        'prevalence_data': prevalence_data,
        'interaction_colors': interaction_colors
    }


def process_lifetime_data(df):
    """Process data for lifetime violin plot with frame ranges."""
    import numpy as np

    # Diccionario de abreviaciones
    interaction_abbreviations = {
        'HBDonor': 'HBD', 'HBAcceptor': 'HBA', 'Cationic': 'Cat',
        'Anionic': 'Ani', 'WaterBridge': 'WB', 'PiStacking': 'πS',
        'PiCation': 'πC', 'PiAnion': 'πA', 'CationPi': 'Cπ',
        'FaceToFace': 'F2F',
        'EdgeToFace': 'E2F', 'MetalDonor': 'MD', 'MetalAcceptor': 'MA',
        'VdWContact': 'VdW', 'CloseContact': 'CC', 'Hydrophobic': 'Hyd',
        'XBAcceptor': 'XBA', 'XBDonor': 'XBD'
    }

    data_list = []

    for idx, row in df.iterrows():
        ts = row['timeseries']

        if isinstance(ts, tuple):
            bit_array = np.array(ts, dtype=int)

        elif isinstance(ts, (list, np.ndarray)):
            bit_array = np.array(ts, dtype=int)

        elif isinstance(ts, str) and set(ts.strip()) <= {'0', '1'}:
            bit_array = np.array([int(x) for x in ts.strip()])
        else:

            print(f"Formato no soportado en fila {idx}: {type(ts)} -> {ts}")
            continue

        one_indices = np.where(bit_array == 1)[0]
        if len(one_indices) > 0:
            intervals = []
            start = one_indices[0]
            prev = start
            for current in one_indices[1:]:
                if current != prev + 1:
                    intervals.append((start, prev + 1))
                    start = current
                prev = current
            intervals.append((start, prev + 1))

            abbrev = interaction_abbreviations.get(row['interaction_name'],
                                                   row['interaction_name'])
            pair_name = f"{row['sel1']} - {row['sel2']} ({abbrev})"
            for start, end in intervals:
                data_list.append({
                    'pair': pair_name,
                    'lifetime': end - start,
                    'interaction_name': row['interaction_name'],
                    'prevalence': row['prevalence'],
                    'start_frame': start,
                    'end_frame': end
                })

    return pd.DataFrame(data_list)

# =============================================================================
#
# =============================================================================
# pickle = '/media/rglez/Roy2TB/Dropbox/RoyData/intermap/tutorial-mayank/outputs/ligs-channel_InterMap.pickle'
# cfg = '/media/rglez/Roy2TB/Dropbox/RoyData/intermap/tutorial-mayank/outputs/ligs-channel_InterMap.cfg'
# self = CSVFilter(pickle, cfg)
# self.universe
# self.master.head()
