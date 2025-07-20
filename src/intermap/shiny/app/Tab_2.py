"""
Tab 2: Prevalence Plot
Contains the PrevalencePlot class for creating sel1 and sel2 prevalence plots.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from intermap.shiny.app.css import all_interactions_colors

"""
def process_prevalence_data2(df, selection_column, sort_by='note'):

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

"""

def process_prevalence_data(df, selection_column, batch_size=250):
    """Process data for interaction plots.

    Args:
        df: DataFrame with interaction data
        selection_column: Column to use for selection ('sel1' for sel1, 'sel2' for sel2)
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


class PrevalencePlot:
    """
    A class to create interactive prevalence plots for molecular interactions.
    """

    def __init__(self, df, plot_size=(800, 600), plot_type='sel1'):
        """
        Initialize the PrevalencePlot with the provided DataFrame.

        Args:
            df (pd.DataFrame): DataFrame containing interaction data
            plot_size (tuple): Size of the plot (width, height)
            plot_type (str): Type of plot ('sel1' for sel1, 'sel2' for sel2)
        """
        self.df = df
        self.width = plot_size[0]
        self.height = plot_size[1]
        self.plot_type = plot_type
        self.selection_column = 'sel1' if plot_type == 'sel1' else 'sel2'
        self.processed_data = None

    def process_data(self):
        """
        Process the DataFrame for prevalence visualization.

        Returns:
            list: Processed data ready for plotting
        """
        self.processed_data = process_prevalence_data(self.df, self.selection_column)
        return self.processed_data

    def create_bar_traces(self):
        """
        Create bar traces for each interaction in the data.

        Returns:
            list: List of plotly bar traces
        """
        if not self.processed_data:
            self.process_data()

        traces = []

        for item in self.processed_data:
            x_value = item['sel1'] if self.plot_type == 'sel1' else item['sel2']
            other_sel = item['sel2'] if self.plot_type == 'sel1' else item['sel1']

            trace = go.Bar(
                name=item['interaction_name'],
                x=[x_value],
                y=[item['prevalence']],
                marker=dict(
                    color=item['color'],
                    line=dict(color='#1a1a1a', width=1)
                ),
                showlegend=item['show_legend'],
                legendgroup=item['interaction_name'],
                hovertemplate=(
                        f"<b>Selection_{1 if self.plot_type == 'sel1' else 2}:</b> %{{x}}<br>" +
                        "<b>Interaction:</b> " + item['interaction_name'] + "<br>" +
                        f"<b>Selection_{2 if self.plot_type == 'sel1' else 1}:</b> " + other_sel + "<br>" +
                        "<b>Prevalence:</b> " + f"{item['prevalence']:.1f}%" + "<br>" +
                        "<b>Annotation:</b> " + item['annotation'] +
                        "<extra></extra>"
                )
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
        x_title = axisx if self.plot_type == 'sel1' else axisy

        return dict(
            width=self.width,
            height=self.height,
            barmode='overlay',
            xaxis_title=dict(
                text=x_title,
                font=dict(family="Roboto", size=16, color='rgb(26, 26, 26)')
            ),
            yaxis_title=dict(
                text="<b>Interaction Prevalence (%)</b>",
                font=dict(family="Roboto", size=16, color='rgb(26, 26, 26)')
            ),
            showlegend=True,
            legend=dict(
                title=dict(
                    text="<b>Interaction Types</b>",
                    font=dict(family="Roboto", size=16, color='rgb(26, 26, 26)')
                ),
                yanchor="top",
                y=0.99,
                xanchor="left",
                x=1.02,
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='rgba(0,0,0,0.2)',
                borderwidth=1,
                font=dict(family="Roboto", size=16, color='rgb(26, 26, 26)')
            ),
            margin=dict(l=50, r=150, t=50, b=50),
            paper_bgcolor='white',
            plot_bgcolor='white',
            bargap=0.15,
            bargroupgap=0.1
        )

    def update_axes_style(self, fig):
        """
        Update axes styling for the plot.

        Args:
            fig: Plotly figure to update
        """
        fig.update_xaxes(
            tickangle=45,
            gridcolor='rgba(0,0,0,0.15)',
            zeroline=True,
            zerolinecolor='rgba(0,0,0,0.25)',
            tickfont=dict(family="Roboto", size=16, color='rgb(26, 26, 26)')
        )

        fig.update_yaxes(
            range=[0, 100],
            gridcolor='rgba(0,0,0,0.15)',
            zeroline=True,
            zerolinecolor='rgba(0,0,0,0.25)',
            tickfont=dict(family="Roboto", size=16, color='rgb(26, 26, 26)')
        )

    def create_prevalence_plot(self, axisx, axisy):
        """
        Create the complete prevalence plot.

        Args:
            axisx (str): X-axis title
            axisy (str): Y-axis title

        Returns:
            go.Figure: Complete plotly figure
        """
        if self.df is None or self.df.empty:
            return None

        fig = go.Figure()

        # Add bar traces
        bar_traces = self.create_bar_traces()
        for trace in bar_traces:
            fig.add_trace(trace)

        # Update layout
        layout = self.create_layout(axisx, axisy)
        fig.update_layout(layout)

        # Update axes styling
        self.update_axes_style(fig)

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
# # sel1
# sel1_plot = PrevalencePlot(master_df, (800, 600), 'sel1')
# fig1 = sel1_plot.create_prevalence_plot("Selection 1", "Selection 2")
#
# # sel2
# sel2_plot = PrevalencePlot(master_df, (800, 600), 'sel2')
# fig2 = sel2_plot.create_prevalence_plot("Selection 1", "Selection 2")
