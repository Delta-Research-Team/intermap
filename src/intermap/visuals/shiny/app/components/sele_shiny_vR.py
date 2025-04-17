# Created by gonzalezroy at 3/25/25
from collections import defaultdict
from os.path import dirname, join

import MDAnalysis as mda
import numpy as np
import pandas as pd
import rgpack.generals as gnl


class TopoDF:
    """
    Class to shrink the InterMap CSV file based on a topology selection
    """

    def __init__(self, csv):
        """
        Initialize the TopoPandas class

        Args:
            csv (str): path to the InterMap CSV file
        """
        self.csv_path = gnl.check_path(csv)
        self.topo_path = self.read_topo_path()
        self.df_master, self.resolution = self.load_csv()
        self.universe = mda.Universe(self.topo_path)
        self.top2df = self.map_topology_to_dataframe()

    def read_topo_path(self):
        """
        Read the topology path from the InterMap CSV file

        Returns:
            str: Path to the topology file
        """
        with open(self.csv_path) as f:
            topo_line = f.readline()
        csv_dirname = dirname(self.csv_path)
        topo_path = join(csv_dirname, topo_line.split('#')[1].strip())

        return gnl.check_path(topo_path)

    def load_csv(self):
        """
        Load the InterMap CSV file into a pandas DataFrame

        Returns:
            pd.DataFrame: DataFrame containing the InterMap data
        """

        df = pd.read_csv(self.csv_path, comment='#')
        N = len(df.iloc[0, 0].split('_'))
        resolution = 'atom' if N == 5 else 'residue'
        return df, resolution

    def _get_res2atoms(self):
        """
        Get the residue indices from the universe

        Returns:
            list: List of residue indices
        """
        resindices = self.universe.atoms.resindices
        res2atoms = defaultdict(set)
        for i, resindex in enumerate(resindices):
            res2atoms[resindex].add(i)
        return res2atoms

    def map_topology_to_dataframe(self):
        """
        Map the topology indices to the DataFrame indices

        Returns:
            dict: Dictionary mapping topology indices to DataFrame indices
        """
        # Automatically determine the number of separators
        N = 4 if self.resolution == 'atom' else 2

        # Extract the indices from the DataFrame
        s1 = self.df_master.iloc[:, 0].str.split('_', expand=True)[N]
        s1 = s1.astype(np.int32).tolist()
        s2 = self.df_master.iloc[:, 1].str.split('_', expand=True)[N]
        s2 = s2.astype(np.int32).tolist()

        if self.df_master.iloc[:, 2].any():
            s3 = self.df_master.iloc[:, 2].str.split('_', expand=True)[N]
            s3 = s3.astype(np.int32).tolist()
        else:
            s3 = None

        # Map the topology indices to the DataFrame indices
        top2df = defaultdict(set)
        for df_index, s1_topindex in enumerate(s1):
            s2_topindex = s2[df_index]

            if self.resolution == 'atom':
                top2df[s1_topindex].add(df_index)
                top2df[s2_topindex].add(df_index)
                if s3:
                    s3_topindex = s3[df_index]
                    top2df[s3_topindex].add(df_index)

            else:
                res2atoms = self._get_res2atoms()
                s1_topindices = res2atoms[s1_topindex]
                for topindex in s1_topindices:
                    top2df[topindex].add(df_index)

                s2_topindices = res2atoms[s2_topindex]
                for topindex in s2_topindices:
                    top2df[topindex].add(df_index)

                if s3:
                    s3_topindex = s3[df_index]
                    top2df[s3_topindex].add(df_index)
        return top2df

    def reselect(self, selection):
        """
        Reselect the dataframe based on the given MDAnalysis selection

        Args:
            selection: MDAnalysis selection string

        Returns:
            filtered_df: DataFrame filtered based on the selection
        """
        # Select atoms based on the provided selection string
        selected_indices = self.universe.select_atoms(selection).indices

        # Get the DataFrame indices corresponding to the selected atoms
        try:
            df_lines = set.union(*(self.top2df[x] for x in selected_indices))
        except TypeError:
            raise TypeError(
                f"Selection '{selection}' returned no atoms. Please check the selection string."
            )

        # Filter the DataFrame based on the selected indices
        filtered_df = self.df_master.iloc[list(df_lines), :].copy()
        return filtered_df


# =============================================================================
#
# =============================================================================

csv = '/home/gonzalezroy/RoyHub/intermap/scripts/identity/runs/mpro/imap/mpro_prot-prot_1000_InterMap.csv'
topo = '/media/gonzalezroy/Expansion/romie/TRAJECTORIES_INPUTS_DATA_mpro_wt_variants_amarolab/a173v/a173v_mpro_chainA_rep123.pr5.aligned_CA.not_waters_or_ions.psf'

self = TopoDF(csv)
mini_df = self.reselect('resid 0')
