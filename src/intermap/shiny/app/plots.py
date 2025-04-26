"""
Plotting components for the InterMap Visualizations app.

"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from intermap.shiny.app.css import all_interactions_colors
from plotly_resampler import FigureResampler

from shiny import ui

def create_plot(df, width, height,
                show_prevalence=False):

    interaction_priority = {
        # Strong interactions
        'Anionic': 1,
        'Cationic': 2,
        'HBDonor': 3,
        'HBAcceptor': 4,

        # Medium-strength interactions
        'MetalDonor': 5,
        'MetalAcceptor': 6,

        # π interactions
        'PiCation': 7,
        'CationPi': 8,
        'PiStacking': 9,
        'FaceToFace': 10,
        'EdgeToFace': 11,

        # Halogen bonding
        'XBDonor': 12,
        'XBAcceptor': 13,

        # Weak interactions
        'Hydrophobic': 14,
        'VdWContact': 15,
        'CloseContact': 16
    }

    # Optimizar el procesamiento de grupos
    df['priority'] = df['interaction_name'].map(
        interaction_priority)
    priority_df = (df.sort_values(
        ['sel1', 'sel2', 'priority', 'prevalence'],
        ascending=[True, True, True, False])
                   .groupby(['sel1', 'sel2'])
                   .first()
                   .reset_index())

    # Pivot de interacciones y permanencias - optimizado
    pivot_interaction = pd.pivot_table(priority_df,
                                       values='interaction_name',
                                       index='sel1',
                                       columns='sel2',
                                       aggfunc='first',
                                       fill_value='')

    pivot_prevalence = pd.pivot_table(priority_df,
                                      values='prevalence',
                                      index='sel1',
                                      columns='sel2',
                                      aggfunc='first',
                                      fill_value="")

    # Pre-calcular valores redondeados
    pivot_prevalence_rounded = pivot_prevalence.round(1).astype(str)

    # Generar la escala de colores para las interacciones activas
    present_interactions = sorted(priority_df['interaction_name'].unique(),
                                  key=lambda x: interaction_priority[x])

    interaction_to_num = {inter: i for i, inter in
                          enumerate(present_interactions)}

    # Crear la figura
    fig = go.Figure()

    # Añadir las interacciones de fondo
    fig.add_trace(go.Heatmap(
        z=[[0] * len(pivot_interaction.columns)] * len(
            pivot_interaction.index),
        x=pivot_interaction.columns,
        y=pivot_interaction.index,
        showscale=False,
        colorscale=[[0, '#FEFBF6'], [1, '#FEFBF6']],
        hoverongaps=False,
        hoverinfo='skip',
        showlegend=False
    ))

    # Crear un heatmap por cada tipo de interacción
    for interaction in present_interactions:
        # Crear máscara para esta interacción usando NumPy para mayor eficiencia
        mask = pivot_interaction.values == interaction
        z_values = np.where(mask, interaction_to_num[interaction], None)

        # Crear matriz de texto eficientemente
        text_matrix = np.where(mask, pivot_prevalence_rounded.values, '')

        fig.add_trace(go.Heatmap(
            z=z_values,
            x=pivot_interaction.columns,
            y=pivot_interaction.index,
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
                font=dict(
                    family="Roboto",
                    size=15,
                    color='rgb(26, 26, 26)'
                )
            ),
            hovertemplate=(
                    "<b>Selection_1:</b> %{y}<br>" +
                    f"<b>Interaction:</b> {interaction}<br>" +
                    "<b>Selection_2:</b> %{x}<br>" +
                    "<b>Prevalence:</b> %{text}%<extra></extra>"
            ),
            legendgroup=interaction,
            showlegend=False,
            visible=True,
            xgap=1,
            ygap=1,
        ))

        # Añadir el elemento de leyenda personalizado
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(
                size=10,
                color=all_interactions_colors[interaction],
                symbol='square',
                line=dict(
                    color='rgba(128, 128, 128, 0.5)',
                    width=1
                )
            ),
            name=f"{interaction} ({interaction_priority[interaction]})",
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
                font=dict(
                    family="Roboto",
                    size=14,
                    color='rgb(26, 26, 26)'
                )
            ),
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.02,
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='rgba(0,0,0,0.2)',
            borderwidth=1,
            font=dict(
                family="Roboto",
                size=14,
                color='rgb(26, 26, 26)'
            )
        ),
        paper_bgcolor='white',
        plot_bgcolor='white',
        xaxis=dict(
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(0, 0, 0, 0.15)',
            linewidth=1,
            linecolor='rgba(0, 0, 0, 0.25)',
            tickangle=45,
            tickfont=dict(
                family="Roboto",
                size=14,
                color='rgb(26, 26, 26)'
            )
        ),
        yaxis=dict(
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(0, 0, 0, 0.15)',
            linewidth=1,
            linecolor='rgba(0, 0, 0, 0.25)',
            tickfont=dict(
                family="Roboto",
                size=14,
                color='rgb(26, 26, 26)'
            )
        ),
        hoverlabel=dict(
            bgcolor='white',
            font_size=14,
            font_family="Roboto",
            align='left'
        )
    )

    return fig


def create_ligand_interactions_plot(df, width, height):
    # 773

    # Ordenar datos usando NumPy
    sort_idx = np.lexsort((
        -df['prevalence'].to_numpy(),
        df['sel1'].to_numpy()
    ))
    receptor_interactions = df.iloc[sort_idx]

    # Crear figura
    fig = go.Figure()
    legend_entries = set()

    # Procesar datos en lotes
    batch_size = 250
    unique_proteins = receptor_interactions['sel1'].unique()

    for i in range(0, len(unique_proteins), batch_size):
        batch_proteins = unique_proteins[i:i + batch_size]
        batch_data = receptor_interactions[
            receptor_interactions['sel1'].isin(batch_proteins)]

        # Procesar cada proteína en el lote
        for _, row in batch_data.iterrows():
            show_legend = row['interaction_name'] not in legend_entries
            if show_legend:
                legend_entries.add(row['interaction_name'])

            solid_color = all_interactions_colors[row['interaction_name']]

            # Añadir barra
            fig.add_trace(go.Bar(
                name=row['interaction_name'],
                x=[row['sel1']],
                y=[row['prevalence']],
                marker=dict(
                    color=solid_color,
                    line=dict(
                        color='#1a1a1a',
                        width=1),
                ),
                showlegend=show_legend,
                legendgroup=row['interaction_name'],
                hovertemplate=(
                        "<b>Selection_1:</b> %{x}<br>"
                        "<b>Interaction:</b> " + row[
                            'interaction_name'] + "<br>" +
                        "<b>Selection_2:</b> " + row['sel2'] + "<br>" +
                        "<b>Prevalence:</b> " + f"{row['prevalence']:.1f}%" +
                        "<extra></extra>"
                )
            ))

    fig.update_layout(
        width=width,
        height=height,
        barmode='overlay',
        xaxis_title="Selection 1 Atoms",
        yaxis_title="Interaction Prevalence (%)",
        showlegend=True,
        legend=dict(
            title=dict(
                text="<b>Interaction Types</b>",
                font=dict(
                    family="Roboto",
                    size=15,
                    color='rgb(26, 26, 26)'
                )
            ),
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.02,
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='rgba(0,0,0,0.2)',
            borderwidth=1,
            font=dict(
                family="Roboto",
                size=15,
                color='rgb(26, 26, 26)'
            ),

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
        tickfont=dict(
            family="Roboto",
            size=14,
            color='rgb(26, 26, 26)'
        )
    )

    fig.update_yaxes(
        range=[0, 100],
        gridcolor='rgba(0,0,0,0.15)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.25)',
        tickfont=dict(
            family="Roboto",
            size=14,
            color='rgb(26, 26, 26)'
        )
    )

    return fig


def create_receptor_interactions_plot(df, width, height):
    # Ordenar datos usando NumPy
    sort_idx = np.lexsort((
        -df['prevalence'].to_numpy(),
        df['sel2'].to_numpy()
    ))
    receptor_interactions = df.iloc[sort_idx]

    # Crear figura
    fig = go.Figure()
    legend_entries = set()

    # Procesar datos en lotes
    batch_size = 250
    unique_proteins = receptor_interactions['sel2'].unique()

    for i in range(0, len(unique_proteins), batch_size):
        batch_proteins = unique_proteins[i:i + batch_size]
        batch_data = receptor_interactions[
            receptor_interactions['sel2'].isin(batch_proteins)]

        # Procesar cada proteína en el lote
        for _, row in batch_data.iterrows():
            show_legend = row['interaction_name'] not in legend_entries
            if show_legend:
                legend_entries.add(row['interaction_name'])

            solid_color = all_interactions_colors[row['interaction_name']]

            # Añadir barra
            fig.add_trace(go.Bar(
                name=row['interaction_name'],
                x=[row['sel2']],
                y=[row['prevalence']],
                marker=dict(
                    color=solid_color,
                    line=dict(
                        color='#1a1a1a',
                        width=1),
                ),
                showlegend=show_legend,
                legendgroup=row['interaction_name'],
                hovertemplate=(
                        "<b>Selection_2:</b> %{x}<br>" +
                        "<b>Interaction:</b> " + row[
                            'interaction_name'] + "<br>" +
                        "<b>Selection_1:</b> " + row['sel1'] + "<br>" +
                        "<b>Prevalence:</b> " + f"{row['prevalence']:.1f}%" +
                        "<extra></extra>"
                )
            ))

    fig.update_layout(
        width=width,
        height=height,
        barmode='overlay',
        xaxis_title="Receptor Atoms",
        yaxis_title="Interaction Prevalence (%)",
        showlegend=True,
        legend=dict(
            title=dict(
                text="<b>Interaction Types</b>",
                font=dict(
                    family="Roboto",
                    size=15,
                    color='rgb(26, 26, 26)'
                )
            ),
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.02,
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='rgba(0,0,0,0.2)',
            borderwidth=1,
            font=dict(
                family="Roboto",
                size=15,
                color='rgb(26, 26, 26)'
            )
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
        tickfont=dict(
            family="Roboto",
            size=14,
            color='rgb(26, 26, 26)'
        )
    )

    fig.update_yaxes(
        range=[0, 100],
        gridcolor='rgba(0,0,0,0.15)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.25)',
        tickfont=dict(
            family="Roboto",
            size=14,
            color='rgb(26, 26, 26)'
        )
    )

    return fig


# Primero, definimos las abreviaturas para cada tipo de interacción
interaction_abbreviations = {
    'HBDonor': 'HBD',
    'HBAcceptor': 'HBA',
    'Cationic': 'Cat',
    'Anionic': 'Ani',
    'Water Bridge': 'WB',
    'PiStacking': 'πS',
    'PiCation': 'πC',
    'CationPi': 'Cπ',
    'FaceToFace': 'F2F',
    'EdgeToFace': 'E2F',
    'MetalDonor': 'MD',
    'MetalAcceptor': 'MA',
    'VdWContact': 'VdW',
    'CloseContact': 'CC',
    'Hydrophobic': 'Hyd',
    'XBAcceptor': 'XBA',
    'XBDonor': 'XBD'
}

from plotly_resampler import FigureResampler


def create_interactions_over_time_plot(df, width, height):
    """
    Create an interactive scatter plot with marginal histograms using plotly-resampler.
    """

    # Crear pares de selecciones con abreviaturas de interacción
    df['selection_pair'] = df['sel1'] + ' - ' + \
                           df['sel2'] + ' (' + \
                           df['interaction_name'].map(
                               interaction_abbreviations) + ')'

    ui.notification_show(
        "We are preparing the plot with your current selections...\n"
        "A lot of data huh  ;)"
        , type="message", duration=20)

    # Procesar series de tiempo
    frame_interactions = []
    for idx, row in df.iterrows():
        timeseries = np.array(list(row['timeseries']), dtype=int)
        frames_with_interaction = np.where(timeseries == 1)[0]
        for frame in frames_with_interaction:
            frame_interactions.append({
                'selection_pair': row['selection_pair'],
                'frame': frame,
                'prevalence': row['prevalence'],
                'interaction_name': row['interaction_name']
            })


    # Crear DataFrame para el scatter plot
    scatter_df = pd.DataFrame(frame_interactions)

    # Crear la figura con subplots
    fig = make_subplots(
        rows=2, cols=2,
        column_widths=[0.8, 0.2],
        row_heights=[0.2, 0.8],
        vertical_spacing=0.02,
        horizontal_spacing=0.02,
        specs=[[{"type": "histogram"}, {"type": "histogram"}],
               [{"type": "scatter"}, {"type": "histogram"}]]
    )

    # Inicializar FigureResampler
    fig = FigureResampler(fig, default_n_shown_samples=2000)

    # Scatter plot principal
    for pair in scatter_df['selection_pair'].unique():
        mask = scatter_df['selection_pair'] == pair
        pair_data = scatter_df[mask]

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
                hovertemplate=(
                        "Frame: %{x}<br>" +
                        "Selection Pair: %{y}<br>" +
                        "Interaction: %{customdata[0]}<br>" +
                        "Prevalence: %{customdata[1]:.1f}%<br>" +
                        "<extra></extra>"
                ),
                customdata=pair_data[['interaction_name', 'prevalence']].values
            ),
            row=2, col=1,
            limit_to_view=True
            # Importante: esto permite el resampling sin reordenar
        )

    ui.notification_show("Hope you like it!", type="message",
                         duration=3)

    # Histograma superior (interacciones por frame)
    fig.add_trace(
        go.Histogram(
            x=scatter_df['frame'],
            marker=dict(
                color='rgb(64,81,181)',
                opacity=0.7,
                line=dict(
                    color='gray',
                    width=1
                )
            ),
            xbins=dict(
                size=1,
                start=scatter_df['frame'].min() - 0.5,
                end=scatter_df['frame'].max() + 0.5
            ),
            name='Interactions per Frame',
            hovertemplate=(
                    "Frame: %{x}<br>" +
                    "n-Inters: %{y}<br>" +
                    "<extra></extra>"
            )
        ),
        row=1, col=1
    )

    # Histograma derecho (prevalencia por par)
    pairs_with_interactions = df['selection_pair'].unique()
    prevalence_data = df.groupby('selection_pair')[
        'prevalence'].mean()

    interaction_colors = []
    for pair in prevalence_data.index:
        interaction_abbrev = pair.split('(')[1].rstrip(')')
        for full_name, abbrev in interaction_abbreviations.items():
            if abbrev == interaction_abbrev:
                interaction_colors.append(all_interactions_colors[full_name])
                break

    fig.add_trace(
        go.Bar(
            y=prevalence_data.index,
            x=prevalence_data.values,
            orientation='h',
            marker=dict(
                color=interaction_colors,
                opacity=0.7,
                line=dict(
                    color='#1a1a1a',
                    width=1
                )
            ),
            name='Prevalence (%)',
            hovertemplate=(
                    "Selection Pair: %{y}<br>" +
                    "Prevalence: %{x:.1f}%<br>" +
                    "<extra></extra>"
            )
        ),
        row=2, col=2
    )


    # Actualizar layout
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

    # Configurar ejes vinculados
    fig.update_xaxes(
        row=2, col=1,
        title_text="Frame Number",
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
        title_text="Selection Pairs",
        gridcolor='rgba(0,0,0,0.15)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.25)',
        domain=[0, 0.75]
    )

    fig.update_yaxes(
        title_text="n-Inters",
        row=1, col=1,
        domain=[0.85, 1]
    )

    fig.update_xaxes(
        title_text="Prevalence (%)",
        row=2, col=2,
        domain=[0.85, 1]
    )

    fig.update_yaxes(
        row=2, col=2,
        matches='y3',
        showticklabels=False,
        domain=[0, 0.75]
    )

    # Actualizar fuentes y estilos
    for axis in fig.layout:
        if axis.startswith('xaxis') or axis.startswith('yaxis'):
            fig.layout[axis].update(
                tickfont=dict(family="Roboto", size=14,
                              color='rgb(26, 26, 26)')
            )

    return fig


# version without resampler
"""
def create_interactions_over_time_plot(df, width, height,
                                       prevalence_threshold):
    if df is None or df.empty:
        return None

    # Filtrado inicial
    mask = (df['interaction_name'].isin(selected_interactions) &
            (df['prevalence'] >= prevalence_threshold))
    df_filtered = df[mask]

    if df_filtered.empty:
        return None

    # Procesamiento de búsqueda
    search_status = get_search_state()
    if search_status['active'] and search_status['search_term']:
        search_term = search_status['search_term'].upper()
        search_mask = (df_filtered['sel1'].str.contains(search_term,
                                                             case=False) |
                       df_filtered['sel2'].str.contains(search_term,
                                                             case=False))
        df_filtered = df_filtered[search_mask]

    if df_filtered.empty:
        return None

    # Crear pares de selecciones con abreviaturas de interacción
    df_filtered['selection_pair'] = df_filtered['sel1'] + ' - ' + \
                                    df_filtered['sel2'] + ' (' + \
                                    df_filtered['interaction_name'].map(interaction_abbreviations) + ')'

    # Procesar series de tiempo
    frame_interactions = []
    for idx, row in df_filtered.iterrows():
        timeseries = np.array(list(row['timeseries']), dtype=int)
        frames_with_interaction = np.where(timeseries == 1)[0]
        for frame in frames_with_interaction:
            frame_interactions.append({
                'selection_pair': row['selection_pair'],
                'frame': frame,
                'prevalence': row['prevalence'],
                'interaction_name': row['interaction_name']
            })

    # Crear DataFrame para el scatter plot
    scatter_df = pd.DataFrame(frame_interactions)

    # Crear la figura con subplots
    fig = make_subplots(
        rows=2, cols=2,
        column_widths=[0.8, 0.2],
        row_heights=[0.2, 0.8],
        vertical_spacing=0.02,
        horizontal_spacing=0.02,
        specs=[[{"type": "histogram"}, {"type": "histogram"}],
               [{"type": "scatter"}, {"type": "histogram"}]]
    )

    # Scatter plot principal
    fig.add_trace(
        go.Scatter(
            x=scatter_df['frame'],
            y=scatter_df['selection_pair'],
            mode='markers',
            marker=dict(
                symbol='square',
                size=8,
                color=scatter_df['interaction_name'].map(all_interactions_colors),
                opacity=0.7
            ),
            name='Interactions',
            hovertemplate=(
                "Frame: %{x}<br>" +
                "Selection Pair: %{y}<br>" +
                "Interaction: %{customdata[0]}<br>" +
                "Prevalence: %{customdata[1]:.1f}%<br>" +
                "<extra></extra>"
            ),
            customdata=scatter_df[['interaction_name', 'prevalence']].values
        ),
        row=2, col=1
    )

    # Histograma superior (interacciones por frame)
    fig.add_trace(
        go.Histogram(
            x=scatter_df['frame'],
            marker=dict(
                color='blue',
                opacity=0.7,
                line=dict(
                    color='#1a1a1a',
                    # Mismo color de contorno que los otros plots
                    width=1
                )
            ),
            xbins=dict(
                size=1,  # Una barra por frame
                start=scatter_df['frame'].min() - 0.5,
                end=scatter_df['frame'].max() + 0.5
            ),
            name='Interactions per Frame',
            hovertemplate=(
                    "Frame: %{x}<br>" +
                    "n-Inters: %{y}<br>" +
                    "<extra></extra>"
            )
        ),
        row=1, col=1
    )

    # Histograma derecho (prevalencia por par)
    # Extraer el tipo de interacción de los paréntesis para cada par
    pairs_with_interactions = df_filtered['selection_pair'].unique()
    prevalence_data = df_filtered.groupby('selection_pair')[
        'prevalence'].mean()

    # Extraer el tipo de interacción de los paréntesis para cada par
    interaction_colors = []
    for pair in prevalence_data.index:
        # Extraer la abreviatura entre paréntesis
        interaction_abbrev = pair.split('(')[1].rstrip(')')
        # Encontrar la interacción completa que corresponde a esta abreviatura
        for full_name, abbrev in interaction_abbreviations.items():
            if abbrev == interaction_abbrev:
                interaction_colors.append(all_interactions_colors[full_name])
                break

    fig.add_trace(
        go.Bar(
            y=prevalence_data.index,
            x=prevalence_data.values,  # Convertir a porcentaje
            orientation='h',
            marker=dict(
                color=interaction_colors,
                # Usar los colores específicos de cada interacción
                opacity=0.7,
                line=dict(
                    color='#1a1a1a',
                    # Mismo color de contorno que los otros plots
                    width=1
                )
            ),
            name='Prevalence (%)',
            hovertemplate=(
                    "Selection Pair: %{y}<br>" +
                    "Prevalence: %{x:.1f}%<br>" +
                    "<extra></extra>"
            )
        ),
        row=2, col=2
    )

    # Actualizar layout
    fig.update_layout(
        width=width,
        height=height,
        showlegend=False,
        paper_bgcolor='white',
        plot_bgcolor='white',
        margin=dict(l=50, r=50, t=50, b=50),
        dragmode='zoom',
        hovermode='closest'
    )

    # Configurar ejes vinculados
    fig.update_xaxes(
        row=2, col=1,
        title_text="Frame Number",
        gridcolor='rgba(0,0,0,0.15)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.25)',
        domain=[0, 0.8]
    )

    fig.update_xaxes(
        row=1, col=1,
        showticklabels=False,
        matches='x3',  # Vincula con el eje X del scatter
        domain=[0, 0.8]
    )

    fig.update_yaxes(
        row=2, col=1,
        title_text="Selection Pairs",
        gridcolor='rgba(0,0,0,0.15)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.25)',
        domain=[0, 0.75]
    )

    # Actualizar histograma superior
    fig.update_yaxes(
        title_text="Count",
        row=1, col=1,
        domain=[0.85, 1]
    )

    # Actualizar histograma derecho
    fig.update_xaxes(
        title_text="Prevalence (%)",
        row=2, col=2,
        domain=[0.85, 1]
    )

    fig.update_yaxes(
        row=2, col=2,
        matches='y3',  # Vincula con el eje Y del scatter
        showticklabels=False,
        domain=[0, 0.75]
    )

    # Actualizar fuentes y estilos
    for axis in fig.layout:
        if axis.startswith('xaxis') or axis.startswith('yaxis'):
            fig.layout[axis].update(
                tickfont=dict(family="Roboto", size=14,
                              color='rgb(26, 26, 26)')
            )

    return fig
"""
