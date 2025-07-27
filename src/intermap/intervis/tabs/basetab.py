# Created by rglez at 7/27/25
"""
BaseTab class for the InterVis application.

This class serves as a base for all tabs in the InterMap application,
providing common functionality and properties.
"""
import pandas as pd
import plotly.graph_objects as go

from intermap.intervis.app.icsv import CSVFilter, interaction_priority


class BaseTab:
    """
    Base class for all tabs in the InterVis application.
    """

    def __init__(self, df, width=600, height=800, **kwargs):
        """
        Initializes the BaseTab with the provided DataFrame and plot dimensions.

        Args:
            df (pd.DataFrame): DataFrame containing interaction data.
            width (int): Width of the plot.
            height (int): Height of the plot.
            **kwargs: Additional keyword arguments for customization.
        """
        self.df = df
        self.width = width
        self.height = height
        self.processed_data = self.process_data()
        self.show_prevalence = kwargs.get('show_prevalence', None)

    def process_data(self):
        """
        Process the DataFrame to prepare it for visualization.

        This method should be overridden in subclasses to implement specific
        data processing logic.
        """
        raise NotImplementedError("Subclasses must implement this method.")


class HeatMap(BaseTab):
    """
    HeatMap plot for visualizing interaction data as a heatmap.
    """

    def process_data(self):
        """
        Process DataFrame for heatmap visualization.

        Returns:
            dict: Processed data containing pivot tables and interaction info
        """
        df = self.df
        df['priority'] = df['interaction_name'].map(interaction_priority)

        ordered_cols = ['sel1', 'sel2', 'priority', 'prevalence']
        ascending = [True, True, True, False]
        priority_df = (df.sort_values(ordered_cols, ascending=ascending)
                       .groupby(['sel1', 'sel2']).first().reset_index())

        pivot_interaction = pd.pivot_table(data=priority_df,
                                           values='interaction_name',
                                           index='sel2', columns='sel1',
                                           aggfunc="first", fill_value='')

        pivot_prevalence = pd.pivot_table(priority_df, values='prevalence',
                                          index='sel2', columns='sel1',
                                          aggfunc='first',
                                          fill_value=0)

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
            'pivot_note1': pivot_note1,
            'pivot_note2': pivot_note2,
            'present_interactions': present_interactions
        }

    def create_base_heatmap(self):
        """
        Create the base heatmap layer.

        Returns:
            go.Heatmap: Base heatmap trace
        """

        data = self.processed_data
        base_trace = go.Heatmap(
            z=[[0] * len(data['pivot_interaction'].columns)] * len(
                data['pivot_interaction'].index),
            x=data['pivot_interaction'].columns,
            y=data['pivot_interaction'].index,
            showscale=False,
            colorscale=[[0, '#FEFBF6'], [1, '#FEFBF6']],
            hoverongaps=False,
            hoverinfo='skip',
            showlegend=False
        )

        return base_trace


# =============================================================================
#
# =============================================================================

pickle = '/home/rglez/RoyHub/intermap/data/tutorial-mayank/outputs/ligs-channel_InterMap.pickle'
cfg = '/home/rglez/RoyHub/intermap/data/tutorial-mayank/outputs/ligs-channel_InterMap.cfg'
csv_obj = CSVFilter(pickle, cfg)
master_df = csv_obj.master

self = HeatMap(master_df, width=800, height=600, show_prevalence=True)
