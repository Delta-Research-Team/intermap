# Created by gonzalezroy at 3/25/25

import MDAnalysis as mda
import pandas as pd


class CSVFilter:
    """
    A class to filter the InterMap CSV file before loading it with Shiny
    """

    def __init__(self, topo, sele, csv):
        """
        Initialize the CSVFilter class

        Args:
            topo (str): path to the topology file
            sele (str): selection string compliant with MDAnalysis syntax
            csv (str): path to the InterMap CSV file
        """
        self.sele = sele
        self.topo = topo
        self.csv = csv

        self.sel_idx = self.get_sel_idx()
        self.df = self.get_df()

    def get_sel_idx(self):
        """
        Get the indices of the selected atoms in the topology file

        Returns:
            indices (set): a set of indices of the selected atoms
        """
        universe = mda.Universe(self.topo)
        ag = universe.select_atoms(self.sele)
        at_names = [f'{x[0]}_{x[1]}_{x[2]}'
                    for x in zip(ag.resnames, ag.resids, ag.names)]
        return set(at_names)

    def filter(self):
        """
        Filter the InterMap CSV file based on the selected indices
        """

        with open(self.csv, 'rt') as f:
            next(f)
            for i, line in enumerate(f):
                s1, s2, s3, itype, prev, time = line.split(',')
                s1_in_sel = s1 in self.sel_idx
                s2_in_sel = s2 in self.sel_idx
                s3_in_sel = s3 in self.sel_idx
                if s1_in_sel or s2_in_sel or s3_in_sel:
                    yield line.split(',')

    def get_df(self):
        """
        Get the filtered CSV file as a pandas DataFrame
        """
        generator = self.filter()
        df = pd.DataFrame(generator,
                          columns=['sel1_atom', 'sel2_atom', 'water_atom',
                                   'interaction_name', 'prevalence',
                                   ' timeseries'])
        return df
