"""
Plotting components for the InterMap Visualizations app.
"""

import pandas as pd

from intermap.intervis.tabs.heatmap import HeatMap
from intermap.intervis.tabs.Tab_2 import PrevalencePlot
from intermap.intervis.tabs.Tab_3 import LifetimePlot
from intermap.intervis.tabs.Tab_4 import TimeSeriesPlot
from intermap.intervis.tabs.Tab_5 import InterNetwork


###############################################################################
# Tab 1
###############################################################################
def create_plot(df, width, height, axisx, axisy, show_prevalence=False):
    """Create interactive heatmap visualization using HeatmapPlot class."""
    if df.empty:
        return None

    heatmap = HeatMap(
        df=df,
        width=width,
        height=height,
        show_prevalence=show_prevalence)

    return heatmap.create_heatmap_plot(axisx, axisy)


###############################################################################
# Tab 2
###############################################################################
def create_sel1_interactions_plot(df, width, height, axisx, axisy):
    """
    Function to create sel1 interactions plot.

    Args:
        df (pd.DataFrame): DataFrame containing interaction data
        width (int): Plot width
        height (int): Plot height
        axisx (str): X-axis title
        axisy (str): Y-axis title

    Returns:
        go.Figure: Plotly figure
    """
    prevalence_plot = PrevalencePlot(df, (width, height), 'sel1')
    return prevalence_plot.create_prevalence_plot(axisx, axisy)


def create_sel2_interactions_plot(df, width, height, axisx, axisy):
    """
    Function to create sel2 interactions plot.

    Args:
        df (pd.DataFrame): DataFrame containing interaction data
        width (int): Plot width
        height (int): Plot height
        axisx (str): X-axis title
        axisy (str): Y-axis title

    Returns:
        go.Figure: Plotly figure
    """
    prevalence_plot = PrevalencePlot(df, (width, height), 'sel2')
    return prevalence_plot.create_prevalence_plot(axisx, axisy)


###############################################################################
# Tab 3
###############################################################################
def create_lifetime_plot(df, width, height, axisx, axisy,
                         show_prevalence=False):
    """Create lifetime violin plot using LifetimePlot class."""
    if df.empty:
        return None

    lifetime_plot = LifetimePlot(
        df=df,
        plot_size=(width, height),
        show_prevalence=show_prevalence)

    return lifetime_plot.create_lifetime_plot(axisx, axisy, show_prevalence)


###############################################################################
# Tab 4
###############################################################################
def create_interactions_over_time_plot(df, width, height, axisx, axisy):
    """
    Convenience function to create interactions over time plot.

    Args:
        df (pd.DataFrame): DataFrame containing interaction data
        width (int): Plot width
        height (int): Plot height
        axisx (str): X-axis title
        axisy (str): Y-axis title

    Returns:
        go.Figure: Plotly figure
    """
    time_series_plot = TimeSeriesPlot(df, (width, height))
    return time_series_plot.create_time_series_plot(axisx, axisy)


###############################################################################
# Tab 5
###############################################################################
def create_network_plot(df, width, height, axisx, axisy, network_params=None):
    """Create interactive network visualization using InterNetwork class."""
    if df.empty:
        return None

    # Default parameters if none provided
    default_params = {
        'gravity': -200,
        'central_gravity': 0.005,
        'spring_length': 200,
        'spring_constant': 0.5,
        'avoid_overlap': 0.8,
        'stabilization_iterations': 1000,
        'physics_enabled': True
    }

    params = network_params if network_params else default_params

    network = InterNetwork(
        master_df=df,
        plot_size=(width, height),
        node_sizes=(
            params.get('min_node_size', 20), params.get('max_node_size', 50)),
        edge_widths=(
            params.get('min_edge_width', 5), params.get('max_edge_width', 15)),
        network_params=params)

    return network.create_network_plot()
