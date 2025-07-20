"""
Tab 1: Heatmap Plot
Contains the HeatmapPlot class for creating interaction heatmaps.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from intermap.shiny.app.css import all_interactions_colors


def process_heatmap_data(df):
    """
    Process DataFrame for heatmap visualization.

    Args:
        df (pd.DataFrame): DataFrame containing interaction data

    Returns:
        dict: Processed data containing pivot tables and interaction info
    """
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
                                      aggfunc='first', fill_value=0)  # Cambiado de "" a 0


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


class HeatmapPlot:
    """
    A class to create interactive heatmap visualizations for molecular interactions.
    """

    def __init__(self, df, plot_size=(800, 600), show_prevalence=False):
        """
        Initialize the HeatmapPlot with the provided DataFrame.

        Args:
            df (pd.DataFrame): DataFrame containing interaction data
            plot_size (tuple): Size of the plot (width, height)
            show_prevalence (bool): Whether to show prevalence values on heatmap
        """
        self.df = df
        self.width = plot_size[0]
        self.height = plot_size[1]
        self.show_prevalence = show_prevalence
        self.processed_data = None

    def process_data(self):
        """
        Process the DataFrame for heatmap visualization.

        Returns:
            dict: Processed data ready for plotting
        """
        self.processed_data = process_heatmap_data(self.df)
        return self.processed_data

    def create_base_heatmap(self):
        """
        Create the base heatmap layer.

        Returns:
            go.Heatmap: Base heatmap trace
        """
        if not self.processed_data:
            self.process_data()

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

    def create_interaction_traces(self):
        """
        Create heatmap traces for each interaction type.

        Returns:
            list: List of plotly traces for interactions
        """
        if not self.processed_data:
            self.process_data()

        data = self.processed_data
        traces = []

        interaction_to_num = {inter: i for i, inter in
                              enumerate(data['present_interactions'])}

        for interaction in data['present_interactions']:
            # Create interaction heatmap trace
            heatmap_trace = self._create_single_interaction_heatmap(
                interaction, interaction_to_num, data
            )
            traces.append(heatmap_trace)

            # Create legend entry
            legend_trace = self._create_legend_entry(
                interaction, interaction_to_num
            )
            traces.append(legend_trace)

        return traces

    def _create_single_interaction_heatmap(self, interaction,
                                           interaction_to_num, data):
        """
        Create heatmap trace for a single interaction type.

        Args:
            interaction (str): Interaction name
            interaction_to_num (dict): Mapping of interactions to numbers
            data (dict): Processed data

        Returns:
            go.Heatmap: Heatmap trace for the interaction
        """
        mask = data['pivot_interaction'].values == interaction
        z_values = np.where(mask, interaction_to_num[interaction], None)

        # Fix para el redondeo - manejar valores numéricos correctamente
        prevalence_values = data['pivot_prevalence'].values
        text_matrix = np.where(
            mask,
            np.where(prevalence_values > 0, np.round(prevalence_values, 1).astype(str), ''),
            ''
        )

        # Create customdata for hover information
        customdata = self._create_customdata(mask, data, interaction)

        return go.Heatmap(
            z=z_values,
            x=data['pivot_interaction'].columns,
            y=data['pivot_interaction'].index,
            text=text_matrix,
            texttemplate="%{text}" if self.show_prevalence else "",
            textfont={"size": 9, "family": "Roboto", "color": "black"},
            showscale=False,
            colorscale=[[0, all_interactions_colors[interaction]],
                        [1, all_interactions_colors[interaction]]],
            hoverongaps=False,
            hoverlabel=dict(
                bgcolor=all_interactions_colors[interaction],
                bordercolor='#1a1a1a',
                font=dict(family="Roboto", size=15, color='rgb(26, 26, 26)')
            ),
            customdata=customdata,
            hovertemplate=(
                    "<b>Sel1:</b> %{customdata[0]}<br>" +
                    "<b>Note1:</b> %{customdata[1]}<br>" +
                    "<b>Sel2:</b> %{customdata[2]}<br>" +
                    "<b>Note2:</b> %{customdata[3]}<br>" +
                    f"<b>Interaction:</b> {interaction}<br>" +
                    "<b>Prevalence:</b> %{text}%<br>" +
                    "<extra></extra>"
            ),
            legendgroup=interaction,
            showlegend=False,
            visible=True,
            xgap=1,
            ygap=1
        )

    def _create_customdata(self, mask, data, interaction):
        """
        Create customdata matrix for hover information.

        Args:
            mask (np.ndarray): Boolean mask for current interaction
            data (dict): Processed data
            interaction (str): Current interaction name

        Returns:
            np.ndarray: Customdata array for hover information
        """
        customdata = []
        for i, sel2 in enumerate(data['pivot_interaction'].index):
            row = []
            for j, sel1 in enumerate(data['pivot_interaction'].columns):
                if mask[i, j]:
                    # Search in original DataFrame
                    row_data = self.df[(self.df['sel1'] == sel1) &
                                       (self.df['sel2'] == sel2) &
                                       (self.df[
                                            'interaction_name'] == interaction)]
                    if not row_data.empty:
                        note1 = row_data['note1'].iloc[0] if pd.notna(
                            row_data['note1'].iloc[0]) else ''
                        note2 = row_data['note2'].iloc[0] if pd.notna(
                            row_data['note2'].iloc[0]) else ''
                        row.append([sel1, note1, sel2, note2])
                    else:
                        row.append(['', '', '', ''])
                else:
                    row.append(['', '', '', ''])
            customdata.append(row)

        return np.array(customdata)

    def _create_legend_entry(self, interaction, interaction_to_num):
        """
        Create legend entry for an interaction.

        Args:
            interaction (str): Interaction name
            interaction_to_num (dict): Mapping of interactions to numbers

        Returns:
            go.Scatter: Legend trace
        """
        return go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(
                size=10,
                color=all_interactions_colors[interaction],
                symbol='square',
                line=dict(color='rgba(128, 128, 128, 0.5)', width=1)
            ),
            name=f"{interaction} ({interaction_to_num[interaction] + 1})",
            showlegend=True,
            legendgroup=interaction
        )

    def create_layout(self, axisx, axisy):
        """
        Create the layout configuration for the plot.

        Args:
            axisx (str): X-axis title
            axisy (str): Y-axis title

        Returns:
            dict: Layout configuration
        """
        return dict(
            width=self.width,
            height=self.height,
            margin=dict(l=50, r=150, t=50, b=50),
            showlegend=True,
            legend=dict(
                title=dict(
                    text="<b>Interaction Types (Priority)</b>",
                    font=dict(family="Roboto", size=14,
                              color='rgb(26, 26, 26)')
                ),
                yanchor="top",
                y=0.99,
                xanchor="left",
                x=1.02,
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='rgba(0,0,0,0.2)',
                borderwidth=1,
                font=dict(family="Roboto", size=14, color='rgb(26, 26, 26)')
            ),
            paper_bgcolor='white',
            plot_bgcolor='white',
            xaxis=dict(
                title=axisx,
                title_font=dict(family="Roboto", size=16,
                                color='rgb(26, 26, 26)'),
                showgrid=True,
                gridwidth=1,
                gridcolor='rgba(0, 0, 0, 0.15)',
                linewidth=1,
                linecolor='rgba(0, 0, 0, 0.25)',
                tickangle=45,
                tickfont=dict(family="Roboto", size=14,
                              color='rgb(26, 26, 26)')
            ),
            yaxis=dict(
                title=axisy,
                title_font=dict(family="Roboto", size=16,
                                color='rgb(26, 26, 26)'),
                showgrid=True,
                gridwidth=1,
                gridcolor='rgba(0, 0, 0, 0.15)',
                linewidth=1,
                linecolor='rgba(0, 0, 0, 0.25)',
                tickfont=dict(family="Roboto", size=14,
                              color='rgb(26, 26, 26)')
            )
        )

    def create_heatmap_plot(self, axisx, axisy):
        """
        Create the complete heatmap plot.

        Args:
            axisx (str): X-axis title
            axisy (str): Y-axis title

        Returns:
            go.Figure: Complete plotly figure
        """
        if self.df is None or self.df.empty:
            return None

        fig = go.Figure()

        # Add base heatmap
        base_trace = self.create_base_heatmap()
        fig.add_trace(base_trace)

        # Add interaction traces
        interaction_traces = self.create_interaction_traces()
        for trace in interaction_traces:
            fig.add_trace(trace)

        # Update layout
        layout = self.create_layout(axisx, axisy)
        fig.update_layout(layout)

        return fig

# =============================================================================
#
# =============================================================================
# from intermap.shiny.app.icsv import CSVFilter
#
# pickle = '/home/fajardo01/03_Fajardo_Hub/02_InterMap/visualizations/data/last_version/hmr_pickle/prot-dna_InterMap.pickle'
# cfg = '/home/fajardo01/03_Fajardo_Hub/02_InterMap/visualizations/data/last_version/hmr_pickle/prot-dna_InterMap.cfg'
# csv_obj = CSVFilter(pickle, cfg)
# master_df = csv_obj.master
#
# heatmap = HeatmapPlot(master_df, (800, 600), show_prevalence=True)
# fig = heatmap.create_heatmap_plot("Selection 1", "Selection 2")
