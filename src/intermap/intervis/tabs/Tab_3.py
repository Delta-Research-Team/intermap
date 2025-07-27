"""
Tab 3: Lifetime Plot Implementation
Contains the LifetimePlot class for creating interaction lifetime box plots.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from intermap.intervis.app.css import all_interactions_colors


def process_lifetime_data(df):
    """Process data for lifetime box plot with frame ranges."""
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


class LifetimePlot:
    """
    A class to create interactive lifetime box plots for molecular interactions.
    """

    def __init__(self, df, plot_size=(800, 600), show_prevalence=False):
        """
        Initialize the LifetimePlot with the provided DataFrame.

        Args:
            df (pd.DataFrame): DataFrame containing interaction data
            plot_size (tuple): Size of the plot (width, height)
            show_prevalence (bool): Whether to show prevalence in the plot
        """
        self.df = df
        self.width = plot_size[0]
        self.height = plot_size[1]
        self.show_prevalence = show_prevalence
        self.processed_data = None

    def process_data(self):
        """
        Process the DataFrame for lifetime visualization.

        Returns:
            pd.DataFrame: Processed data ready for plotting
        """
        self.processed_data = process_lifetime_data(self.df)
        return self.processed_data

    def create_box_traces(self):
        """
        Create box traces for each interaction type.

        Returns:
            list: List of plotly box traces
        """
        if self.processed_data is None or self.processed_data.empty:
            return []

        traces = []

        # Convert to numpy arrays for better performance
        pairs = self.processed_data['pair'].values
        lifetimes = self.processed_data['lifetime'].values
        interactions = self.processed_data['interaction_name'].values
        prevalences = self.processed_data['prevalence'].values
        start_frames = self.processed_data['start_frame'].values
        end_frames = self.processed_data['end_frame'].values

        shown_interactions = set()

        unique_interactions = np.unique(interactions)
        for interaction_name in unique_interactions:
            mask = interactions == interaction_name
            show_legend = interaction_name not in shown_interactions
            shown_interactions.add(interaction_name)

            trace = go.Box(
                x=pairs[mask],
                y=lifetimes[mask],
                name=interaction_name,
                marker=dict(
                    color=all_interactions_colors[interaction_name],
                    opacity=0.7,
                    size=4,
                    outliercolor='rgba(50,50,50,0.7)'
                ),
                line=dict(
                    color=all_interactions_colors[interaction_name],
                    width=2
                ),
                boxmean=True,
                notched=False,
                boxpoints='outliers',
                jitter=0.3,
                whiskerwidth=0.7,
                fillcolor='rgba(255,255,255,0.5)',

                customdata=np.column_stack((
                    prevalences[mask],
                    start_frames[mask],
                    end_frames[mask]
                )),
                hovertemplate=(
                    "<b>%{x}</b><br><br>"
                    f"Interaction: {interaction_name}<br>"
                    "Prevalence: %{customdata[0]:.1f}%<br>"
                    "Lifetime: %{y} frames<br>"
                    "Frame range: %{customdata[1]} - %{customdata[2]}<br>"
                    "<extra></extra>"
                ),
                showlegend=show_legend,
                legendgroup=interaction_name,
                hoverlabel=dict(
                    bgcolor=all_interactions_colors[interaction_name],
                    font_size=14,
                    font_family="Roboto"
                ),
            )
            traces.append(trace)

        return traces

    def create_layout(self, axisx, axisy):
        """
        Create the layout configuration for the plot.

        Args:
            axisx (str): X-axis title
            axisy (str): Y-axis title

        Returns:
            dict: Layout configuration
        """
        return {
            'width': self.width,
            'height': self.height,
            'title': {
                'text': "<b>Interaction Lifetimes Distribution</b>",
                'x': 0.5,
                'font': {'family': "Roboto", 'size': 20}
            },
            'showlegend': True,
            'legend': {
                'title': {'text': "<b>Interaction Types</b>",
                          'font': {'family': "Roboto", 'size': 14}},
                'yanchor': "top",
                'y': 0.99,
                'xanchor': "left",
                'x': 1.02,
                'bgcolor': 'rgba(255,255,255,0.9)',
                'bordercolor': 'rgba(0,0,0,0.2)',
                'borderwidth': 1
            },
            'paper_bgcolor': 'white',
            'plot_bgcolor': 'white',
            'boxmode': 'group',
            'margin': {'l': 50, 'r': 150, 't': 50, 'b': 50},
            'xaxis': {
                'title': "<b>Selection Pairs</b>",
                'tickangle': 45,
                'title_font': {'family': "Roboto", 'size': 16},
                'tickfont': {'family': "Roboto", 'size': 14},
                'showgrid': False,
                'zeroline': False,
                'linewidth': 1.5,
                'linecolor': '#d3d3d3',
                'mirror': True
            },
            'yaxis': {
                'title': "<b>Interaction Lifetime (frames)</b>",
                'title_font': {'family': "Roboto", 'size': 16},
                'tickfont': {'family': "Roboto", 'size': 14},
                'showgrid': True,
                'gridcolor': 'rgba(200,200,200,0.2)',
                'zeroline': False,
                'linewidth': 1.5,
                'linecolor': '#d3d3d3',
                'mirror': True,
                'range': [0, None]
            }
        }

    def create_lifetime_plot(self, axisx, axisy, show_prevalence=False):
        """
        Create the complete lifetime box plot.

        Args:
            axisx (str): X-axis title
            axisy (str): Y-axis title
            show_prevalence (bool): Whether to show prevalence (for compatibility)

        Returns:
            go.Figure: Complete plotly figure
        """
        if self.df is None or self.df.empty:
            return None

        # Process data
        self.process_data()

        if self.processed_data is None or self.processed_data.empty:
            return None

        fig = go.Figure()

        # Add box traces
        box_traces = self.create_box_traces()
        for trace in box_traces:
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
# # Para lifetime
# lifetime_plot = LifetimePlot(master_df, (800, 600), show_prevalence=True)
# fig = lifetime_plot.create_lifetime_plot("Selection 1", "Selection 2")
