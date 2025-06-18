"""
Plotting components for the InterMap Visualizations app.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly_resampler import FigureResampler

from intermap.shiny.app.css import all_interactions_colors

from intermap.shiny.app.icsv import process_heatmap_data, process_prevalence_data, process_time_series_data
from rdkit.sping.colors import white


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

        # Crear matrices de customdata
        customdata = []
        for i, sel2 in enumerate(data['pivot_interaction'].index):
            row = []
            for j, sel1 in enumerate(data['pivot_interaction'].columns):
                if mask[i, j]:
                    # Buscar en el DataFrame original
                    row_data = df[(df['sel1'] == sel1) &
                                  (df['sel2'] == sel2) &
                                  (df['interaction_name'] == interaction)]
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

        customdata = np.array(customdata)

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

    # Modificamos el layout para usar axisx y axisy
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
            title=f"",
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
            title=f"",
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
                    "<b>Selection_1:</b> %{x}<br>" +
                    "<b>Interaction:</b> " + item[
                        'interaction_name'] + "<br>" +
                    "<b>Selection_2:</b> " + item['sel2'] + "<br>" +
                    "<b>Prevalence:</b> " + f"{item['prevalence']:.1f}%" + "<br>" +
                    "<b>Annotation:</b> " + item['annotation'] +
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
                    "<b>Prevalence:</b> " + f"{item['prevalence']:.1f}%" + "<br>" +
                    "<b>Annotation:</b> " + item['annotation'] +
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


    for interaction in sorted(data['scatter_df']['interaction_name'].unique()):
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
            showlegend=False,
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
            showlegend=False,
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


def process_lifetime_data(df):
    """Process data for lifetime violin plot with frame ranges."""
    from bitarray import bitarray as ba, util as bu
    import numpy as np

    # Dictionary for interaction name abbreviations
    interaction_abbreviations = {
        'HBDonor': 'HBD', 'HBAcceptor': 'HBA', 'Cationic': 'Cat',
        'Anionic': 'Ani', 'Water Bridge': 'WB', 'PiStacking': 'πS',
        'PiCation': 'πC', 'CationPi': 'Cπ', 'FaceToFace': 'F2F',
        'EdgeToFace': 'E2F', 'MetalDonor': 'MD', 'MetalAcceptor': 'MA',
        'VdWContact': 'VdW', 'CloseContact': 'CC', 'Hydrophobic': 'Hyd',
        'XBAcceptor': 'XBA', 'XBDonor': 'XBD'
    }

    # Convert DataFrame columns to numpy arrays for faster access
    sel1_array = df['sel1'].values
    sel2_array = df['sel2'].values
    interaction_array = df['interaction_name'].values
    prevalence_array = df['prevalence'].values
    timeseries_array = df['timeseries'].values

    data_list = []

    for idx in range(len(df)):
        array = ba(timeseries_array[idx])
        intervals = bu.intervals(array)
        # Filter intervals where contact exists (flag == 1)
        contact_intervals = [x for x in intervals if x[0] == 1]

        if contact_intervals:
            abbrev = interaction_abbreviations.get(interaction_array[idx],
                                                 interaction_array[idx])
            pair_name = f"{sel1_array[idx]} - {sel2_array[idx]} ({abbrev})"

            # Convert contact_intervals to numpy array for faster operations
            intervals_array = np.array(contact_intervals)

            # Extract all starts and ends at once
            starts = intervals_array[:, 1]
            ends = intervals_array[:, 2]
            lifetimes = ends - starts

            # Create arrays for repeated values using numpy
            n_intervals = len(intervals_array)
            pairs = np.repeat(pair_name, n_intervals)
            interactions = np.repeat(interaction_array[idx], n_intervals)
            prevalences = np.repeat(prevalence_array[idx], n_intervals)

            # Add all intervals at once
            for i in range(n_intervals):
                data_list.append({
                    'pair': pairs[i],
                    'lifetime': lifetimes[i],
                    'interaction_name': interactions[i],
                    'prevalence': prevalences[i],
                    'start_frame': starts[i],
                    'end_frame': ends[i]
                })

    return pd.DataFrame(data_list)


def create_lifetime_plot(df, width, height, show_prevalence=False):
    """Create fully vectorized violin plot for interaction lifetimes."""
    # Process data
    processed_df = process_lifetime_data(df)

    if processed_df.empty:
        return None

    fig = go.Figure()

    # Convert to numpy arrays once for better performance
    pairs = processed_df['pair'].values
    lifetimes = processed_df['lifetime'].values
    interactions = processed_df['interaction_name'].values
    prevalences = processed_df['prevalence'].values
    start_frames = processed_df['start_frame'].values
    end_frames = processed_df['end_frame'].values

    # Keep track of shown interactions for legend
    shown_interactions = set()

    # Create violin plot for each interaction type while maintaining original order
    unique_interactions = np.unique(interactions)
    for interaction_name in unique_interactions:
        # Create mask using numpy for better performance
        mask = interactions == interaction_name
        show_legend = interaction_name not in shown_interactions
        shown_interactions.add(interaction_name)

        fig.add_trace(go.Violin(
            x=pairs[mask],
            y=lifetimes[mask],
            name=interaction_name,
            fillcolor=all_interactions_colors[interaction_name],
            line=dict(width=2,
                     color=all_interactions_colors[interaction_name]),
            meanline=dict(visible=True, color='rgba(0,0,0,0.5)', width=2),
            points='outliers',
            pointpos=0,
            marker=dict(
                color='rgba(50,50,50,0.7)',
                size=4,
                symbol='circle'
            ),
            box=dict(visible=True, width=0.1),
            opacity=0.7,
            side='both',
            width=0.7,
            spanmode='soft',
            scalemode='width',
            hoverlabel=dict(
                bgcolor=all_interactions_colors[interaction_name],
                font_size=14,
                font_family="Roboto"
            ),
            # Include start and end frames in customdata
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
            legendgroup=interaction_name
        ))

    # Layout configuration
    layout_config = {
        'width': width,
        'height': height,
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
        'violinmode': 'overlay',
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
            'showgrid': False,
            'zeroline': False,
            'linewidth': 1.5,
            'linecolor': '#d3d3d3',
            'mirror': True
        }
    }

    fig.update_layout(layout_config)

    return fig
