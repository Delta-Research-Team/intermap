"""
Plotting components for the InterMap Visualizations app.

"""

import pandas as pd
import plotly.graph_objects as go
import numpy as np
from ..config import all_interactions_colors
from ..utils.helpers import calculate_prevalence

# Variable global para el estado de búsqueda
search_state = None


# Al inicio del archivo, después de las importaciones, agregar:
COLUMN_MAPPING = {
    'sel1_atom': 'ligand',
    'sel2_atom': 'protein',
    'interaction_name': 'interaction'
}


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

    # Clonar el DataFrame y renombrar columnas
    df = df.copy()
    if 'ligand' not in df.columns:  # Si aún no están renombradas
        df = df.rename(columns=COLUMN_MAPPING)


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
        'CloseContacts': 16
    }

    # Filtrado de interacciones seleccionadas
    df_filtered = df[df['interaction'].isin(selected_interactions)]
    df_filtered['prevalence'] = df_filtered.apply(calculate_prevalence, axis=1)

    search_status = get_search_state()
    if search_status['active'] and search_status['search_term']:
        search_term = search_status['search_term'].upper()
        mask = (df_filtered['ligand'].str.contains(search_term, case=False) |
                df_filtered['protein'].str.contains(search_term, case=False))
        df_filtered = df_filtered[mask]

    if df_filtered.empty:
        return None

    df_filtered = df_filtered[
        df_filtered['prevalence'] >= prevalence_threshold]

    priority_interactions = []

    for (lig, prot), group in df_filtered.groupby(['ligand', 'protein']):
        if not group.empty:
            group['priority'] = group['interaction'].map(interaction_priority)
            group = group.sort_values(['priority', 'prevalence'],
                                      ascending=[True, False])

            priority_int = group.iloc[0]
            priority_interactions.append({
                'ligand': lig,
                'protein': prot,
                'interaction': priority_int['interaction'],
                'prevalence': priority_int['prevalence']
            })

    priority_df = pd.DataFrame(priority_interactions)

    # Pivot de interacciones y permanencias
    pivot_interaction = pd.pivot_table(priority_df,
                                       values='interaction',
                                       index='ligand',
                                       columns='protein',
                                       aggfunc='first',
                                       fill_value='')

    pivot_prevalence = pd.pivot_table(priority_df,
                                      values='prevalence',
                                      index='ligand',
                                      columns='protein',
                                      aggfunc='first',
                                      fill_value=0)

    fig = go.Figure()

    # Añadir las interacciones de fondo, incluso las desactivadas
    fig.add_trace(go.Heatmap(
        z=[[0] * len(pivot_interaction.columns)] * len(
            pivot_interaction.index),
        x=pivot_interaction.columns,
        y=pivot_interaction.index,
        showscale=False,
        colorscale=[[0, '#FEFBF6'], [1, '#FEFBF6']],
        hoverongaps=False,
        hoverinfo='skip',
    ))

    # Generar la escala de colores para las interacciones activas
    present_interactions = sorted(priority_df['interaction'].unique(),
                                  key=lambda x: interaction_priority[x])

    color_scale = []
    for i, interaction in enumerate(present_interactions):
        pos = i / (len(present_interactions) - 1) if len(
            present_interactions) > 1 else 0.5
        color_scale.extend([[pos, all_interactions_colors[interaction]],
                            [pos, all_interactions_colors[interaction]]])

    if len(present_interactions) == 1:
        color = all_interactions_colors[present_interactions[0]]
        color_scale = [[0, color], [1, color]]

    interaction_to_num = {inter: i for i, inter in
                          enumerate(present_interactions)}
    pivot_numerical = pivot_interaction.replace(interaction_to_num)

    # Mostrar la permanencia si es necesario
    text_values = pivot_prevalence.round(1).astype(str)
    if show_prevalence:
        text_values = text_values.mask(pivot_prevalence == 0, '')
    else:
        text_values = pivot_interaction.applymap(lambda x: '')

    fig.add_trace(go.Heatmap(
        z=pivot_numerical,
        x=pivot_interaction.columns,
        y=pivot_interaction.index,
        text=text_values,
        texttemplate="%{text}",
        textfont={"size": 9, "family": "Arial", "color": "black"},
        showscale=False,
        colorscale=color_scale,
        hoverongaps=False,
        hovertemplate=(
                "<b>Selection_1:</b> %{y}<br>" +
                "<b>Selection_2:</b> %{x}<br>" +
                "<b>Interaction:</b> %{customdata}<br>" +
                "<b>Prevalence:</b> %{text}%<extra></extra>"
        ),
        customdata=pivot_interaction.values,
        xgap=1,
        ygap=1,
    ))

    # Añadir las leyendas para las interacciones
    for interaction in present_interactions:
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(
                size=10,
                color=all_interactions_colors[interaction],
                symbol='square',
                line=dict(  # Añadir borde a los marcadores de la leyenda
                    color='rgba(128, 128, 128, 0.5)',
                    width=1
                )
            ),
            name=f"{interaction} ({interaction_priority[interaction]})",
            showlegend=True
        ))

    fig.update_layout(
        width=width,
        height=height,
        margin=dict(l=50, r=50, t=50, b=50),
        showlegend=True,
        legend=dict(
            title=dict(
                text="<b>Interaction Types (Priority)</b>",
                font=dict(
                    family="Ubuntu Mono",
                    size=12,
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
                family="Ubuntu Mono",
                size=12,
                color='rgb(26, 26, 26)'
            )
        ),
        xaxis=dict(
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(0, 0, 0, 0.15)',
            linewidth=1,
            linecolor='rgba(0, 0, 0, 0.25)',
            tickfont=dict(
                family="Ubuntu Mono",
                size=12,
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
                family="Ubuntu Mono",
                size=12,
                color='rgb(26, 26, 26)'
            )
        )
    )

    return fig


def create_ligand_interactions_plot(df, width, height,
                                    selected_interactions,
                                    prevalence_threshold):
    if df is None or df.empty:
        return None

    # Clonar el DataFrame y renombrar columnas
    df = df.copy()
    if 'ligand' not in df.columns:  # Si aún no están renombradas
        df = df.rename(columns=COLUMN_MAPPING)


    # Convertir a arrays de NumPy para operaciones más rápidas
    df_filtered = df[df['interaction'].isin(selected_interactions)].copy()

    # Vectorizar el cálculo de permanencia
    if 'frames_present' in df_filtered.columns and 'total_frames' in df_filtered.columns:
        df_filtered['prevalence'] = (df_filtered['frames_present'].to_numpy() /
                                     df_filtered[
                                         'total_frames'].to_numpy()) * 100
    else:
        df_filtered['prevalence'] = df_filtered.apply(calculate_prevalence,
                                                      axis=1)

    # Optimizar búsqueda
    search_status = get_search_state()
    if search_status['active'] and search_status['search_term']:
        search_term = search_status['search_term'].upper()
        mask = df_filtered['ligand'].str.contains(search_term, case=False)
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
        df_filtered['ligand'].to_numpy()
    ))
    receptor_interactions = df_filtered.iloc[sort_idx]

    # Crear figura
    fig = go.Figure()
    legend_entries = set()

    # Procesar datos en lotes
    batch_size = 50
    unique_proteins = receptor_interactions['ligand'].unique()

    for i in range(0, len(unique_proteins), batch_size):
        batch_proteins = unique_proteins[i:i + batch_size]
        batch_data = receptor_interactions[
            receptor_interactions['ligand'].isin(batch_proteins)]

        # Procesar cada proteína en el lote
        for _, row in batch_data.iterrows():
            show_legend = row['interaction'] not in legend_entries
            if show_legend:
                legend_entries.add(row['interaction'])

            solid_color = all_interactions_colors[row['interaction']]

            # Añadir barra
            fig.add_trace(go.Bar(
                name=row['interaction'],
                x=[row['ligand']],
                y=[row['prevalence']],
                marker=dict(
                    color=solid_color,
                    line=dict(
                        color='#1a1a1a',
                        width=1),
                ),
                showlegend=show_legend,
                legendgroup=row['interaction'],
                hovertemplate=(
                        "<b>Selection_1:</b> %{x}<br>"
                        "<b>Interaction:</b> " + row['interaction'] + "<br>" +
                        "<b>Selection_2:</b> " + row['protein'] + "<br>" +
                        "<b>Prevalence:</b> " + f"{row['prevalence']:.1f}%" +
                        "<extra></extra>"
                )
            ))

    fig.update_layout(
        width=width,
        height=height,
        barmode='overlay',
        xaxis_title="Ligand Atoms",
        yaxis_title="Interaction Prevalence (%)",
        showlegend=True,
        legend=dict(
            title=dict(
                text="<b>Interaction Types</b>",
                font=dict(
                    family="Ubuntu Mono",
                    size=12,
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
                family="Ubuntu Mono",
                size=12,
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
            family="Ubuntu Mono",
            size=12,
            color='rgb(26, 26, 26)'
        )
    )

    fig.update_yaxes(
        range=[0, 100],
        gridcolor='rgba(0,0,0,0.15)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.25)',
        tickfont=dict(
            family="Ubuntu Mono",
            size=12,
            color='rgb(26, 26, 26)'
        )
    )

    return fig


def create_receptor_interactions_plot(df, width, height,
                                      selected_interactions,
                                      prevalence_threshold):
    if df is None or df.empty:
        return None

    # Clonar el DataFrame y renombrar columnas
    df = df.copy()
    if 'ligand' not in df.columns:  # Si aún no están renombradas
        df = df.rename(columns=COLUMN_MAPPING)

    # Convertir a arrays de NumPy para operaciones más rápidas
    df_filtered = df[df['interaction'].isin(selected_interactions)].copy()

    # Vectorizar el cálculo de permanencia
    if 'frames_present' in df_filtered.columns and 'total_frames' in df_filtered.columns:
        df_filtered['prevalence'] = (df_filtered['frames_present'].to_numpy() /
                                     df_filtered[
                                         'total_frames'].to_numpy()) * 100
    else:
        df_filtered['prevalence'] = df_filtered.apply(calculate_prevalence,
                                                      axis=1)

    # Optimizar búsqueda
    search_status = get_search_state()
    if search_status['active'] and search_status['search_term']:
        search_term = search_status['search_term'].upper()
        mask = df_filtered['protein'].str.contains(search_term, case=False)
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
        df_filtered['protein'].to_numpy()
    ))
    receptor_interactions = df_filtered.iloc[sort_idx]

    # Crear figura
    fig = go.Figure()
    legend_entries = set()

    # Procesar datos en lotes
    batch_size = 50
    unique_proteins = receptor_interactions['protein'].unique()

    for i in range(0, len(unique_proteins), batch_size):
        batch_proteins = unique_proteins[i:i + batch_size]
        batch_data = receptor_interactions[
            receptor_interactions['protein'].isin(batch_proteins)]

        # Procesar cada proteína en el lote
        for _, row in batch_data.iterrows():
            show_legend = row['interaction'] not in legend_entries
            if show_legend:
                legend_entries.add(row['interaction'])

            solid_color = all_interactions_colors[row['interaction']]

            # Añadir barra
            fig.add_trace(go.Bar(
                name=row['interaction'],
                x=[row['protein']],
                y=[row['prevalence']],
                marker=dict(
                    color=solid_color,
                    line=dict(
                        color='#1a1a1a',
                        width=1),
                ),
                showlegend=show_legend,
                legendgroup=row['interaction'],
                hovertemplate=(
                        "<b>Selection_2:</b> %{x}<br>" +
                        "<b>Interaction:</b> " + row['interaction'] + "<br>" +
                        "<b>Selection_1:</b> " + row['ligand'] + "<br>" +
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
                    family="Ubuntu Mono",
                    size=12,
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
                family="Ubuntu Mono",
                size=12,
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
            family="Ubuntu Mono",
            size=12,
            color='rgb(26, 26, 26)'
        )
    )

    fig.update_yaxes(
        range=[0, 100],
        gridcolor='rgba(0,0,0,0.15)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.25)',
        tickfont=dict(
            family="Ubuntu Mono",
            size=12,
            color='rgb(26, 26, 26)'
        )
    )

    return fig


def create_interactions_over_time_plot(df, width, height, selected_interactions,
                                     prevalence_threshold):
    """
    Create a plot showing interactions over time frames.
    Shows all interactions present in each frame with a unified hover display.
    """
    if df is None or df.empty:
        return None

    # Filtrar por interacciones seleccionadas usando el nombre correcto de la columna
    df_filtered = df[df['interaction_name'].isin(selected_interactions)].copy()

    # Aplicar filtro de prevalencia
    df_filtered['prevalence'] = df_filtered.apply(calculate_prevalence, axis=1)
    df_filtered = df_filtered[df_filtered['prevalence'] >= prevalence_threshold]

    if df_filtered.empty:
        return None

    # Aplicar filtro de búsqueda si está activo
    search_status = get_search_state()
    if search_status['active'] and search_status['search_term']:
        search_term = search_status['search_term'].upper()
        mask = (df_filtered['sel1_atom'].str.contains(search_term, case=False) |
                df_filtered['sel2_atom'].str.contains(search_term, case=False))
        df_filtered = df_filtered[mask]

    if df_filtered.empty:
        return None

    # Obtener columnas de frames (todas las columnas numéricas)
    frame_columns = [col for col in df_filtered.columns if str(col).isdigit()]

    if not frame_columns:
        return None

    # Crear diccionario para almacenar información por frame y por tipo
    frame_data = {}
    all_interaction_types = set()

    # Analizar cada frame
    for frame_idx, frame_col in enumerate(frame_columns):
        # Obtener datos donde la interacción está presente (==1)
        frame_interactions = df_filtered[df_filtered[frame_col] == 1]

        # Contar interacciones por tipo
        interaction_counts = frame_interactions['interaction_name'].value_counts()

        # Solo almacenar frames donde hay interacciones
        if not interaction_counts.empty:
            frame_data[frame_idx] = {
                'counts': interaction_counts,
                'total': interaction_counts.sum()
            }
            all_interaction_types.update(interaction_counts.index)

    # Crear figura
    fig = go.Figure()

    # Añadir un trace por cada tipo de interacción
    for interaction_type in sorted(all_interaction_types):
        x_values = []
        y_values = []
        hover_texts = []

        for frame_idx, data in frame_data.items():
            if interaction_type in data['counts']:
                count = data['counts'][interaction_type]
                x_values.append(frame_idx)
                y_values.append(count)

                # Construir texto del hover para este punto
                other_interactions = "<br>".join([
                    f"• <span style='color: {all_interactions_colors[int_type]}'>{int_type}</span>: {count}"
                    for int_type, count in data['counts'].items()
                ])

                hover_text = (
                    f"<b>Frame {frame_idx}</b><br><br>"
                    f"<b>Interactions Present:</b><br>"
                    f"{other_interactions}<br><br>"
                    f"<b>Total Interactions:</b> {data['total']}"
                )
                hover_texts.append(hover_text)

        # Añadir trace para este tipo de interacción
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
                    line=dict(
                        color='white',
                        width=1
                    )
                ),
                text=hover_texts,
                hoverinfo='text',
                hoverlabel=dict(
                    bgcolor='rgba(255,255,255,0.95)',
                    bordercolor=all_interactions_colors[interaction_type],
                    font=dict(
                        family='Arial',
                        size=12,
                        color='black'
                    )
                )
            ))

    # Actualizar el diseño
    fig.update_layout(
        width=width,
        height=height,
        xaxis_title="Frame Number",
        yaxis_title="Number of Interactions",
        showlegend=True,
        legend=dict(
            title=dict(
                text="<b>Interaction Types</b>",
                font=dict(
                    family="Ubuntu Mono",
                    size=12,
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
                family="Ubuntu Mono",
                size=12,
                color='rgb(26, 26, 26)'
            )
        ),
        margin=dict(l=50, r=150, t=50, b=50),
        paper_bgcolor='white',
        plot_bgcolor='white',
        hovermode='closest',
        hoverdistance=100,
        hoverlabel=dict(
            font_family="Ubuntu Mono",
            font_size=12,
            align='left'
        )
    )

    fig.update_xaxes(
        gridcolor='rgba(0,0,0,0.15)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.25)',
        tickfont=dict(
            family="Ubuntu Mono",
            size=12,
            color='rgb(26, 26, 26)'
        ),
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
        tickfont=dict(
            family="Ubuntu Mono",
            size=12,
            color='rgb(26, 26, 26)'
        ),
        showgrid=True,
        nticks=8,
        tickmode='auto',
        minor_griddash='dot',
        minor_gridwidth=1,
        minor_gridcolor='rgba(0,0,0,0.1)'
    )

    return fig
