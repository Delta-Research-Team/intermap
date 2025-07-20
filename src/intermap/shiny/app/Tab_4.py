"""
Tab 3: Time Series Plot Implementation
Contains the TimeSeriesPlot class for creating interactions over time plots.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly_resampler import FigureResampler

from intermap.shiny.app.css import all_interactions_colors


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


class TimeSeriesPlot:
    """
    A class to create interactive time series plots for molecular interactions over time.
    """

    def __init__(self, df, plot_size=(800, 600)):
        """
        Initialize the TimeSeriesPlot with the provided DataFrame.

        Args:
            df (pd.DataFrame): DataFrame containing interaction data
            plot_size (tuple): Size of the plot (width, height)
        """
        self.df = df
        self.width = plot_size[0]
        self.height = plot_size[1]
        self.processed_data = None

    def process_data(self):
        """
        Process the DataFrame for time series visualization.

        Returns:
            dict: Processed data ready for plotting
        """
        self.processed_data = process_time_series_data(self.df)
        return self.processed_data

    def create_subplot_figure(self):
        """
        Create the subplot structure for the time series plot.

        Returns:
            FigureResampler: Plotly figure with subplots configured
        """
        fig = make_subplots(
            rows=2, cols=2,
            column_widths=[0.8, 0.2],
            row_heights=[0.2, 0.8],
            vertical_spacing=0.02,
            horizontal_spacing=0.02,
            specs=[[{"type": "histogram"}, {"type": "histogram"}],
                   [{"type": "scatter"}, {"type": "histogram"}]]
        )

        return FigureResampler(fig, default_n_shown_samples=2000)

    def add_legend_traces(self, fig):
        """
        Add legend traces for interaction types.

        Args:
            fig: Plotly figure to add traces to
        """
        if not self.processed_data:
            self.process_data()

        for interaction in sorted(self.processed_data['scatter_df'][
                                      'interaction_name'].unique()):
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode='markers',
                    marker=dict(
                        symbol='square',
                        size=10,
                        color=all_interactions_colors[interaction]
                    ),
                    name=interaction,
                    showlegend=True
                )
            )

    def add_scatter_traces(self, fig):
        """
        Add scatter traces for interaction data points.

        Args:
            fig: Plotly figure to add traces to
        """
        if not self.processed_data:
            self.process_data()

        for pair in self.processed_data['scatter_df'][
            'selection_pair'].unique():
            pair_data = self.processed_data['scatter_df'][
                self.processed_data['scatter_df']['selection_pair'] == pair]

            fig.add_trace(
                go.Scattergl(
                    x=pair_data['frame'],
                    y=[pair] * len(pair_data),
                    mode='markers',
                    marker=dict(
                        symbol='square',
                        size=8,
                        color=pair_data['interaction_name'].map(
                            all_interactions_colors),
                        opacity=0.7
                    ),
                    name=pair,
                    showlegend=False,
                    customdata=pair_data[
                        ['interaction_name', 'prevalence']].values,
                    hovertemplate=(
                        "Frame: %{x}<br>"
                        "Selection Pair: %{y}<br>"
                        "Interaction: %{customdata[0]}<br>"
                        "Prevalence: %{customdata[1]:.1f}%"
                        "<extra></extra>"
                    )
                ),
                row=2, col=1,
                limit_to_view=True
            )

    def add_histogram_trace(self, fig):
        """
        Add histogram trace for frame distribution.

        Args:
            fig: Plotly figure to add traces to
        """
        if not self.processed_data:
            self.process_data()

        fig.add_trace(
            go.Histogram(
                x=self.processed_data['scatter_df']['frame'],
                marker=dict(
                    color='rgb(64,81,181)',
                    opacity=0.7,
                    line=dict(color='gray', width=1)
                ),
                xbins=dict(
                    size=1,
                    start=self.processed_data['scatter_df'][
                              'frame'].min() - 0.5,
                    end=self.processed_data['scatter_df']['frame'].max() + 0.5
                ),
                name='Interactions per Frame',
                showlegend=False,
                hovertemplate=(
                    "Frame: %{x}<br>"
                    "n-Inters: %{y}<br>"
                    "<extra></extra>"
                )
            ),
            row=1, col=1
        )

    def add_prevalence_bar_trace(self, fig):
        """
        Add prevalence bar trace.

        Args:
            fig: Plotly figure to add traces to
        """
        if not self.processed_data:
            self.process_data()

        fig.add_trace(
            go.Bar(
                y=self.processed_data['prevalence_data'].index,
                x=self.processed_data['prevalence_data'].values,
                orientation='h',
                marker=dict(
                    color=self.processed_data['interaction_colors'],
                    opacity=0.7,
                    line=dict(color='#1a1a1a', width=1)
                ),
                name='Prevalence (%)',
                showlegend=False,
                hovertemplate=(
                    "Selection Pair: %{y}<br>"
                    "Prevalence: %{x:.1f}%<br>"
                    "<extra></extra>"
                )
            ),
            row=2, col=2
        )

    def configure_layout(self, fig, axisx, axisy):
        """
        Configure the layout for the time series plot.

        Args:
            fig: Plotly figure to configure
            axisx (str): X-axis title
            axisy (str): Y-axis title
        """
        fig.update_layout(
            width=self.width,
            height=self.height * 1.2,
            showlegend=True,
            paper_bgcolor='white',
            plot_bgcolor='white',
            margin=dict(l=50, r=150, t=50, b=50),
            dragmode='zoom',
            hovermode='closest',
            legend=dict(
                title="<b>Interaction Types</b>",
                yanchor="top",
                y=0.99,
                xanchor="left",
                x=1.02,
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='rgba(0,0,0,0.2)',
                borderwidth=1,
                font=dict(family="Roboto", size=14, color='rgb(26, 26, 26)')
            )
        )

    def configure_axes(self, fig):
        """
        Configure axes styling for the plot.

        Args:
            fig: Plotly figure to configure
        """
        # Main scatter plot axes
        fig.update_xaxes(
            row=2, col=1,
            title=dict(
                text="<b>Frame Number</b>",
                font=dict(family="Roboto", size=16, color='rgb(26, 26, 26)')
            ),
            gridcolor='rgba(0,0,0,0.15)',
            zeroline=True,
            zerolinecolor='rgba(0,0,0,0.25)',
            domain=[0, 0.8]
        )

        fig.update_yaxes(
            row=2, col=1,
            title=dict(
                text="<b>Selection Pairs</b>",
                font=dict(family="Roboto", size=16, color='rgb(26, 26, 26)')
            ),
            gridcolor='rgba(0,0,0,0.15)',
            zeroline=True,
            zerolinecolor='rgba(0,0,0,0.25)',
            domain=[0, 0.75]
        )

        # Histogram axes
        fig.update_xaxes(
            row=1, col=1,
            showticklabels=False,
            matches='x3',
            domain=[0, 0.8]
        )

        fig.update_yaxes(
            title=dict(
                text="<b>n-Inters</b>",
                font=dict(family="Roboto", size=16, color='rgb(26, 26, 26)')
            ),
            row=1, col=1,
            domain=[0.85, 1]
        )

        # Prevalence bar axes
        fig.update_xaxes(
            title=dict(
                text="<b>Prevalence (%)</b>",
                font=dict(family="Roboto", size=16, color='rgb(26, 26, 26)')
            ),
            row=2, col=2,
            domain=[0.85, 1]
        )

        fig.update_yaxes(
            row=2, col=2,
            matches='y3',
            showticklabels=False,
            domain=[0, 0.75]
        )

        # Update all axis fonts
        for axis in fig.layout:
            if axis.startswith('xaxis') or axis.startswith('yaxis'):
                fig.layout[axis].update(
                    tickfont=dict(family="Roboto", size=16,
                                  color='rgb(26, 26, 26)')
                )

    def create_time_series_plot(self, axisx, axisy):
        """
        Create the complete time series plot.

        Args:
            axisx (str): X-axis title
            axisy (str): Y-axis title

        Returns:
            go.Figure: Complete plotly figure
        """
        if self.df is None or self.df.empty:
            return None

        # Create subplot structure
        fig = self.create_subplot_figure()

        # Process data
        self.process_data()

        # Add all traces
        self.add_legend_traces(fig)
        self.add_scatter_traces(fig)
        self.add_histogram_trace(fig)
        self.add_prevalence_bar_trace(fig)

        # Configure layout and axes
        self.configure_layout(fig, axisx, axisy)
        self.configure_axes(fig)

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
# # time series
# time_series_plot = TimeSeriesPlot(master_df, (800, 600))
# fig = time_series_plot.create_time_series_plot("Selection 1", "Selection 2")
