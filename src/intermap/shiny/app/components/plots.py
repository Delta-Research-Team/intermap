"""
Plotting components for the InterMap Visualizations app.

"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..css import all_interactions_colors
search_state = None

def get_search_state():
    """
    Safely get the current search state.

    Returns:
        dict: Current search state with 'active' and 'search_term' keys
    """
    if search_state is None:
        return {'active': False, 'search_term': ''}
    return search_state.get()


def initialize_search_state(state):
    """
    Initialize the global search state for the plots.

    Args:
        state: A reactive value object to store search state
    """
    global search_state
    search_state = state



def create_plot(df, width, height, selected_interactions, prevalence_threshold,
                show_prevalence=False):
    if df is None or df.empty:
        return None

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

    # Optimizar el filtrado inicial usando NumPy
    df_filtered = df[df['interaction_name'].isin(selected_interactions)]
    df_filtered = df_filtered[
        df_filtered['prevalence'] >= prevalence_threshold]

    search_status = get_search_state()
    if search_status['active'] and search_status['search_term']:
        search_term = search_status['search_term'].upper()
        mask = (df_filtered['sel1'].str.contains(search_term,
                                                      case=False) |
                df_filtered['sel2'].str.contains(search_term, case=False))
        df_filtered = df_filtered[mask]

    if df_filtered.empty:
        return None

    # Optimizar el procesamiento de grupos
    df_filtered['priority'] = df_filtered['interaction_name'].map(
        interaction_priority)
    priority_df = (df_filtered.sort_values(
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


def create_ligand_interactions_plot(df, width, height,
                                    selected_interactions,
                                    prevalence_threshold):
    if df is None or df.empty:
        return None

    # Convertir a arrays de NumPy para operaciones más rápidas
    df_filtered = df[df['interaction_name'].isin(selected_interactions)]

    # Optimizar búsqueda
    search_status = get_search_state()
    if search_status['active'] and search_status['search_term']:
        search_term = search_status['search_term'].upper()
        mask = df_filtered['sel1'].str.contains(search_term, case=False)
        if mask.any():
            df_filtered = df_filtered[mask]
        else:
            return None

    # Filtrado por umbral usando NumPy
    df_filtered = df_filtered[
        df_filtered['prevalence'] >= prevalence_threshold]

    if df_filtered.empty:
        return None

    # Ordenar datos usando NumPy
    sort_idx = np.lexsort((
        -df_filtered['prevalence'].to_numpy(),
        df_filtered['sel1'].to_numpy()
    ))
    receptor_interactions = df_filtered.iloc[sort_idx]

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


def create_receptor_interactions_plot(df, width, height,
                                      selected_interactions,
                                      prevalence_threshold):
    if df is None or df.empty:
        return None

    # Convertir a arrays de NumPy para operaciones más rápidas
    df_filtered = df[df['interaction_name'].isin(selected_interactions)]

    # Optimizar búsqueda
    search_status = get_search_state()
    if search_status['active'] and search_status['search_term']:
        search_term = search_status['search_term'].upper()
        mask = df_filtered['sel2'].str.contains(search_term, case=False)
        if mask.any():
            df_filtered = df_filtered[mask]
        else:
            return None

    # Filtrado por umbral usando NumPy
    df_filtered = df_filtered[
        df_filtered['prevalence'] >= prevalence_threshold]

    if df_filtered.empty:
        return None

    # Ordenar datos usando NumPy
    sort_idx = np.lexsort((
        -df_filtered['prevalence'].to_numpy(),
        df_filtered['sel2'].to_numpy()
    ))
    receptor_interactions = df_filtered.iloc[sort_idx]

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


def create_interactions_over_time_plot(df, width, height,
                                       selected_interactions,
                                       prevalence_threshold):
    """
    Create an interactive scatter plot with marginal histograms using plotly.
    """
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

    # Crear pares de selecciones
    df_filtered['selection_pair'] = df_filtered['sel1'] + ' - ' + \
                                    df_filtered['sel2']

    # Procesar series de tiempo
    frame_interactions = []
    for idx, row in df_filtered.iterrows():
        timeseries = np.array(list(row['timeseries']), dtype=int)
        frames_with_interaction = np.where(timeseries == 1)[0]
        for frame in frames_with_interaction:
            frame_interactions.append({
                'selection_pair': row['selection_pair'],
                'frame': frame,
                'prevalence': row['prevalence']
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
                color='black'
            ),
            name='Interactions'
        ),
        row=2, col=1
    )

    # Histograma superior (interacciones por frame)
    fig.add_trace(
        go.Histogram(
            x=scatter_df['frame'],
            nbinsx=50,
            name='Interactions per Frame'
        ),
        row=1, col=1
    )

    # Histograma derecho (prevalencia por par)
    prevalence_data = df_filtered.groupby('selection_pair')[
        'prevalence'].mean()
    fig.add_trace(
        go.Bar(
            y=prevalence_data.index,
            x=prevalence_data.values * 100,  # Convertir a porcentaje
            orientation='h',
            name='Prevalence (%)'
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
        margin=dict(l=50, r=50, t=50, b=50)
    )

    # Actualizar ejes
    fig.update_xaxes(
        title_text="Frame Number",
        gridcolor='rgba(0,0,0,0.15)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.25)',
        row=2, col=1
    )

    fig.update_yaxes(
        title_text="Selection Pairs",
        gridcolor='rgba(0,0,0,0.15)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.25)',
        row=2, col=1
    )

    # Actualizar histograma superior
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=1)

    # Actualizar histograma derecho
    fig.update_xaxes(title_text="Prevalence (%)", row=2, col=2)
    fig.update_yaxes(showticklabels=False, row=2, col=2)

    return fig


"""

def create_interactions_over_time_plot(df, width, height,
                                       selected_interactions,
                                       prevalence_threshold):

    if df is None or df.empty:
        return None

    # Filtrado vectorizado usando una sola máscara
    mask = (df['interaction_name'].isin(selected_interactions) &
            (df['prevalence'] >= prevalence_threshold))
    df_filtered = df[mask]

    if df_filtered.empty:
        return None

    # Filtrado de búsqueda vectorizado
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

    # Convertir las series de tiempo de string a arrays numpy
    timeseries_arrays = [np.array(list(ts), dtype=int) for ts in
                         df_filtered['timeseries']]
    n_frames = len(timeseries_arrays[0])  # Número de frames

    # Crear diccionario para almacenar información por frame
    frame_data = {}
    interactions_array = df_filtered['interaction_name'].values

    # Procesar cada frame de manera vectorizada
    for frame_idx in range(n_frames):
        # Obtener máscara de interacciones activas en este frame
        frame_mask = np.array([ts[frame_idx] == 1 for ts in timeseries_arrays])

        if np.any(frame_mask):
            # Obtener interacciones activas en este frame
            frame_interactions = interactions_array[frame_mask]

            # Contar interacciones de manera vectorizada
            unique, counts = np.unique(frame_interactions, return_counts=True)
            counts_dict = dict(zip(unique, counts))

            if counts_dict:
                frame_data[frame_idx] = {
                    'counts': pd.Series(counts_dict),
                    'total': len(frame_interactions)
                }

    # Crear figura
    fig = go.Figure()

    # Procesar cada tipo de interacción
    unique_interactions = df_filtered['interaction_name'].unique()

    for interaction_type in sorted(unique_interactions):
        x_values = []
        y_values = []
        hover_texts = []

        for frame_idx, data in frame_data.items():
            if interaction_type in data['counts']:
                count = data['counts'][interaction_type]
                x_values.append(frame_idx)
                y_values.append(count)

                # Construir texto del hover
                interactions_text = "<br>".join(
                    f"• <span style='color: {all_interactions_colors[int_type]}'>{int_type}</span>: {data['counts'].get(int_type, 0)}"
                    for int_type in sorted(data['counts'].index)
                )

                hover_text = (
                    f"<b>Frame {frame_idx}</b><br><br>"
                    f"<b>Interactions Present:</b><br>"
                    f"{interactions_text}<br><br>"
                    f"<b>Total Interactions:</b> {data['total']}"
                )
                hover_texts.append(hover_text)

        if x_values:  # Solo añadir si hay datos para este tipo de interacción
            fig.add_trace(go.Scatter(
                x=x_values,
                y=y_values,
                mode='markers',
                name=interaction_type,
                marker=dict(
                    size=10,
                    color=all_interactions_colors[interaction_type],
                    symbol='circle',
                    line=dict(color='white', width=1)
                ),
                text=hover_texts,
                hoverinfo='text',
                hoverlabel=dict(
                    bgcolor='rgba(255,255,255,0.95)',
                    bordercolor=all_interactions_colors[interaction_type],
                    font=dict(family='Arial', size=14, color='black')
                )
            ))

    # El resto del código de actualización del diseño permanece igual
    fig.update_layout(
        width=width,
        height=height,
        xaxis_title="Frame Number",
        yaxis_title="Number of Interactions",
        showlegend=True,
        legend=dict(
            title=dict(
                text="<b>Interaction Types</b>",
                font=dict(family="RobotoRoboto", size=14,
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
        margin=dict(l=50, r=150, t=50, b=50),
        paper_bgcolor='white',
        plot_bgcolor='white',
        hovermode='closest',
        hoverdistance=100,
        hoverlabel=dict(
            font_family="Roboto",
            font_size=14,
            align='left'
        )
    )

    fig.update_xaxes(
        gridcolor='rgba(0,0,0,0.15)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.25)',
        tickfont=dict(family="Roboto", size=14, color='rgb(26, 26, 26)'),
        showgrid=True,
        nticks=10,
        tickmode='auto',
        minor_griddash='dot',
        minor_gridwidth=1,
        minor_gridcolor='rgba(0,0,0,0.1)'
    )

    fig.update_yaxes(
        gridcolor='rgba(0,0,0,0.15)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.25)',
        tickfont=dict(family="Roboto", size=14, color='rgb(26, 26, 26)'),
        showgrid=True,
        nticks=8,
        tickmode='auto',
        minor_griddash='dot',
        minor_gridwidth=1,
        minor_gridcolor='rgba(0,0,0,0.1)'
    )

    return fig

"""
