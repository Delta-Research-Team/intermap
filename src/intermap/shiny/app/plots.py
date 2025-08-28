"""
Plotting components for the InterMap Visualizations app.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly_resampler import FigureResampler

from intermap.shiny.app.css import all_interactions_colors


def process_heatmap_data(df):
    """Process data for heatmap plot."""
    interaction_priority = {
        'Anionic': 1, 'Cationic': 2, 'HBDonor': 3, 'HBAcceptor': 4,
        'MetalDonor': 5, 'MetalAcceptor': 6, 'PiCation': 7, 'CationPi': 8,
        'PiStacking': 9, 'FaceToFace': 10, 'EdgeToFace': 11, 'XBDonor': 12,
        'XBAcceptor': 13, 'Hydrophobic': 14, 'VdWContact': 15,
        'CloseContact': 16
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

    present_interactions = sorted(priority_df['interaction_name'].unique(),
                                  key=lambda x: interaction_priority[x])

    return {
        'pivot_interaction': pivot_interaction,
        'pivot_prevalence': pivot_prevalence,
        'pivot_prevalence_rounded': pivot_prevalence_rounded,
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
                'color': all_interactions_colors[row['interaction_name']]
            })

    return batched_data


def process_time_series_data(df):
    """Process data for time series plot."""
    interaction_abbreviations = {
        'HBDonor': 'HBD', 'HBAcceptor': 'HBA', 'Cationic': 'Cat',
        'Anionic': 'Ani', 'Water Bridge': 'WB', 'PiStacking': 'πS',
        'PiCation': 'πC', 'CationPi': 'Cπ', 'FaceToFace': 'F2F',
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


def create_plot(df, width, height, show_prevalence=False):
    """Create heatmap plot."""
    data = process_heatmap_data(df)

    fig = go.Figure()

    fig.add_trace(go.Heatmap(
        z=[[0] * len(data['pivot_interaction'].columns)] * len(
            data['pivot_interaction'].index),
        x=data['pivot_interaction'].columns,
        y=data['pivot_interaction'].index,
        showscale=False,
        colorscale=[[0, '#FEFBF6'], [1, '#FEFBF6']],
        hoverongaps=False,
        hoverinfo='skip',
        showlegend=False
    ))

    interaction_to_num = {inter: i for i, inter in
                          enumerate(data['present_interactions'])}

    for interaction in data['present_interactions']:
        mask = data['pivot_interaction'].values == interaction
        z_values = np.where(mask, interaction_to_num[interaction], None)
        text_matrix = np.where(mask, data['pivot_prevalence_rounded'].values,
                               '')

        fig.add_trace(go.Heatmap(
            z=z_values,
            x=data['pivot_interaction'].columns,
            y=data['pivot_interaction'].index,
            text=text_matrix,
            texttemplate="%{text}" if show_prevalence else "",
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
            hovertemplate=(
                    "<b>Selection_2:</b> %{y}<br>" +
                    f"<b>Interaction:</b> {interaction}<br>" +
                    "<b>Selection_1:</b> %{x}<br>" +
                    "<b>Prevalence:</b> %{text}%<extra></extra>"
            ),
            legendgroup=interaction,
            showlegend=False,
            visible=True,
            xgap=1,
            ygap=1
        ))

        fig.add_trace(go.Scatter(
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
        ))

    fig.update_layout(
        width=width,
        height=height,
        margin=dict(l=50, r=150, t=50, b=50),
        showlegend=True,
        legend=dict(
            title=dict(
                text="<b>Interaction Types (Priority)</b>",
                font=dict(family="Roboto", size=14, color='rgb(26, 26, 26)')
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
            title="<b>Selection 1</b>",
            title_font=dict(family="Roboto", size=16, color='rgb(26, 26, 26)'),
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(0, 0, 0, 0.15)',
            linewidth=1,
            linecolor='rgba(0, 0, 0, 0.25)',
            tickangle=45,
            tickfont=dict(family="Roboto", size=14, color='rgb(26, 26, 26)')
        ),
        yaxis=dict(
            title="<b>Selection 2</b>",
            title_font=dict(family="Roboto", size=16, color='rgb(26, 26, 26)'),
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(0, 0, 0, 0.15)',
            linewidth=1,
            linecolor='rgba(0, 0, 0, 0.25)',
            tickfont=dict(family="Roboto", size=14, color='rgb(26, 26, 26)')
        )
    )

    return fig


def create_ligand_interactions_plot(df, width, height):
    """Create ligand interactions plot."""
    batched_data = process_prevalence_data(df, 'sel1')

    fig = go.Figure()

    for item in batched_data:
        fig.add_trace(go.Bar(
            name=item['interaction_name'],
            x=[item['sel1']],
            y=[item['prevalence']],
            marker=dict(
                color=item['color'],
                line=dict(color='#1a1a1a', width=1)
            ),
            showlegend=item['show_legend'],
            legendgroup=item['interaction_name'],
            hovertemplate=(
                    "<b>Selection_1:</b> %{x}<br>"
                    "<b>Interaction:</b> " + item[
                        'interaction_name'] + "<br>" +
                    "<b>Selection_2:</b> " + item['sel2'] + "<br>" +
                    "<b>Prevalence:</b> " + f"{item['prevalence']:.1f}%" +
                    "<extra></extra>"
            )
        ))

    fig.update_layout(
        width=width,
        height=height,
        barmode='overlay',
        xaxis_title=dict(
            text="<b>Selection 1</b>",
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

    return fig


def create_receptor_interactions_plot(df, width, height):
    """Create receptor interactions plot."""
    batched_data = process_prevalence_data(df, 'sel2')

    fig = go.Figure()

    for item in batched_data:
        fig.add_trace(go.Bar(
            name=item['interaction_name'],
            x=[item['sel2']],
            y=[item['prevalence']],
            marker=dict(
                color=item['color'],
                line=dict(color='#1a1a1a', width=1)
            ),
            showlegend=item['show_legend'],
            legendgroup=item['interaction_name'],
            hovertemplate=(
                    "<b>Selection_2:</b> %{x}<br>" +
                    "<b>Interaction:</b> " + item[
                        'interaction_name'] + "<br>" +
                    "<b>Selection_1:</b> " + item['sel1'] + "<br>" +
                    "<b>Prevalence:</b> " + f"{item['prevalence']:.1f}%" +
                    "<extra></extra>"
            )
        ))

    fig.update_layout(
        width=width,
        height=height,
        barmode='overlay',
        xaxis_title=dict(
            text="<b>Selection 2</b>",
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

    return fig


def create_interactions_over_time_plot(df, width, height):
    """Create interactions over time plot."""
    data = process_time_series_data(df)

    fig = make_subplots(
        rows=2, cols=2,
        column_widths=[0.8, 0.2],
        row_heights=[0.2, 0.8],
        vertical_spacing=0.02,
        horizontal_spacing=0.02,
        specs=[[{"type": "histogram"}, {"type": "histogram"}],
               [{"type": "scatter"}, {"type": "histogram"}]]
    )

    fig = FigureResampler(fig, default_n_shown_samples=2000)

    for pair in data['scatter_df']['selection_pair'].unique():
        pair_data = data['scatter_df'][
            data['scatter_df']['selection_pair'] == pair]

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

    fig.add_trace(
        go.Histogram(
            x=data['scatter_df']['frame'],
            marker=dict(
                color='rgb(64,81,181)',
                opacity=0.7,
                line=dict(color='gray', width=1)
            ),
            xbins=dict(
                size=1,
                start=data['scatter_df']['frame'].min() - 0.5,
                end=data['scatter_df']['frame'].max() + 0.5
            ),
            name='Interactions per Frame',
            hovertemplate=(
                "Frame: %{x}<br>"
                "n-Inters: %{y}<br>"
                "<extra></extra>"
            )
        ),
        row=1, col=1
    )

    fig.add_trace(
        go.Bar(
            y=data['prevalence_data'].index,
            x=data['prevalence_data'].values,
            orientation='h',
            marker=dict(
                color=data['interaction_colors'],
                opacity=0.7,
                line=dict(color='#1a1a1a', width=1)
            ),
            name='Prevalence (%)',
            hovertemplate=(
                "Selection Pair: %{y}<br>"
                "Prevalence: %{x:.1f}%<br>"
                "<extra></extra>"
            )
        ),
        row=2, col=2
    )

    fig.update_layout(
        width=width,
        height=height * 1.2,
        showlegend=False,
        paper_bgcolor='white',
        plot_bgcolor='white',
        margin=dict(l=50, r=50, t=50, b=50),
        dragmode='zoom',
        hovermode='closest'
    )

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

    fig.update_xaxes(
        row=1, col=1,
        showticklabels=False,
        matches='x3',
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

    fig.update_yaxes(
        title=dict(
            text="<b>n-Inters</b>",
            font=dict(family="Roboto", size=16, color='rgb(26, 26, 26)')
        ),
        row=1, col=1,
        domain=[0.85, 1]
    )

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

    for axis in fig.layout:
        if axis.startswith('xaxis') or axis.startswith('yaxis'):
            fig.layout[axis].update(
                tickfont=dict(family="Roboto", size=16,
                              color='rgb(26, 26, 26)')
            )

    return fig
