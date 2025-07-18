"""
Plotting components for the InterMap Visualizations app.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly_resampler import FigureResampler

from intermap.shiny.app.css import all_interactions_colors

from intermap.shiny.app.icsv import process_heatmap_data, process_prevalence_data, process_time_series_data, process_lifetime_data
from rdkit.sping.colors import white
import networkx as nx
from pyvis.network import Network

from intermap.shiny.app.graph import InterNetwork


def create_plot(df, width, height, axisx, axisy, show_prevalence=False):
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
            title=axisx,
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
            title=axisy,
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

def create_ligand_interactions_plot(df, width, height, axisx, axisy):
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
            text=axisx,
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

def create_receptor_interactions_plot(df, width, height, axisx, axisy):
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
            text=axisy,
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

def create_interactions_over_time_plot(df, width, height, axisx, axisy):
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

def create_lifetime_plot(df, width, height, axisx, axisy, show_prevalence=False):
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


def create_network_plot(df, width, height, axisx, axisy):
    """Create interactive network visualization using InterNetwork class."""
    if df.empty:
        return None

    network = InterNetwork(
        master_df=df,
        plot_size=(width, height),
        node_sizes=(20, 50),
        edge_widths=(5, 15))

    return network.create_network_plot()



###############################################################################
# 3D AREA
###############################################################################

import MDAnalysis as mda
from MDAnalysis.coordinates.PDB import PDBWriter
import tempfile
import os

class InteractiveVisualization:
    def __init__(self, df, universe):
        self.df = df
        self.universe = universe
        self.current_frame = 0
        self.viewer = None
        self.max_frames = self.get_max_frames()

    def get_max_frames(self):
        if self.df.empty:
            return 0

        max_len = 0
        for _, row in self.df.iterrows():
            timeseries = row['timeseries']
            if isinstance(timeseries, str):
                max_len = max(max_len, len(timeseries))

        return max_len - 1 if max_len > 0 else 0

    def update_frame(self, frame_index):
        if frame_index < 0 or frame_index > self.max_frames:
            print(f"Frame {frame_index} fuera de rango (0-{self.max_frames})")
            return

        self.current_frame = frame_index

        df_filtered = filter_interactions_by_frame(self.df, frame_index)

        self.viewer = create_3d_view(
            df_filtered, self.universe,
            width=900, height=700,
            show_protein=True,
            show_interactions=True,
            molecule_style="cartoon",
            sphere_size=0.8,
            line_width=0.15,
            frame_index=frame_index
        )

        return self.viewer


def create_aligned_pdb_from_trajectory(universe, frame_index=0):
    """
    Crear PDB alineado con las coordenadas actuales del trajectory

    Args:
        universe: Universo MDAnalysis
        frame_index: Frame de la trayectoria a usar (0 = primer frame)

    Returns:
        str: Contenido PDB alineado
    """
    try:
        universe.trajectory[frame_index]

        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb',
                                                delete=False)
        temp_filename = temp_file.name
        temp_file.close()
        with PDBWriter(temp_filename) as writer:
            writer.write(universe.atoms)
        with open(temp_filename, 'r') as f:
            pdb_content = f.read()
        os.unlink(temp_filename)
        return pdb_content

    except Exception as e:
        return None


def verify_coordinates_alignment(universe, pdb_content):
    """
    Verificar que las coordenadas estén alineadas

    Args:
        universe: Universo MDAnalysis
        pdb_content: Contenido PDB generado

    Returns:
        bool: True si están alineadas
    """
    try:
        lines = pdb_content.split('\n')
        atom_lines = [line for line in lines if line.startswith('ATOM')]

        if len(atom_lines) == 0:
            return False

        first_atom_line = atom_lines[0]
        pdb_x = float(first_atom_line[30:38])
        pdb_y = float(first_atom_line[38:46])
        pdb_z = float(first_atom_line[46:54])

        mda_pos = universe.atoms[0].position

        diff = np.linalg.norm(
            [pdb_x - mda_pos[0], pdb_y - mda_pos[1], pdb_z - mda_pos[2]])

        print(f"Diferencia de coordenadas: {diff:.4f} Å")

        if diff < 0.01:
            print("Coordenadas perfectamente alineadas")
            return True
        else:
            print("Coordenadas no alineadas")
            return False

    except Exception as e:
        print(f"Error verificando alineación: {e}")
        return False


def create_3d_view(df, universe, width=900, height=760, show_protein=True,
                   show_interactions=True, molecule_style="cartoon",
                   sphere_size=0.5, line_width=0.15, frame_index=0,
                   apply_transparency=False):
    """
    Crear visualización 3D con control de transparencia opcional
    """
    try:
        import py3Dmol

        viewer = py3Dmol.view(width=width, height=height)
        viewer.setBackgroundColor('white')

        pdb_content = create_aligned_pdb_from_trajectory(universe, frame_index)
        if not pdb_content:
            return None

        viewer.addModel(pdb_content, "pdb")
        viewer.setStyle({}, {})

        if show_protein:
            style = get_protein_style(molecule_style)
            viewer.setStyle({'model': -1}, style)
            print(f"✅ Estilo aplicado: {molecule_style}")

        interaction_atoms = set()
        if show_interactions and df is not None and not df.empty:
            interaction_atoms = add_interactions_with_dashed_lines(
                viewer, df, universe, sphere_size, line_width
            )

        if apply_transparency and show_protein and interaction_atoms:
            apply_atom_transparency(viewer, interaction_atoms, molecule_style)

        add_hover_labels(viewer, universe)

        viewer.zoomTo()
        viewer.render()

        print("Visualización creada")
        return viewer

    except Exception as e:
        print(f"Error en create_3d_view: {e}")
        return None

def get_protein_style(molecule_style):
    """
    Obtener estilo de proteína - verificar que funcione correctamente
    """
    styles = {
        "cartoon": {"cartoon": {"color": "spectrum"}},
        "stick": {"stick": {"radius": 0.2, "color": "spectrum"}},
        "sphere": {"sphere": {"radius": 0.5, "color": "spectrum"}},
        "ball_and_stick": {
            "stick": {"radius": 0.2, "color": "spectrum"},
            "sphere": {"radius": 0.3, "color": "spectrum"}
        },
        "cartoon_and_stick": {
            "cartoon": {"color": "spectrum"},
            "stick": {"radius": 0.2, "color": "spectrum"}
        }
    }

    selected_style = styles.get(molecule_style,
                                {"cartoon": {"color": "spectrum"}})
    print(f"Estilo seleccionado: {molecule_style} -> {selected_style}")

    return selected_style

def add_interactions_with_dashed_lines(viewer, df, universe, sphere_size,
                                       line_width):
    """
    Añadir interacciones y recopilar átomos participantes
    """
    print(f"🔗 Procesando {len(df)} interacciones...")

    added_count = 0
    interaction_atoms = set()  # Para rastrear átomos participantes

    for idx, row in df.iterrows():
        try:
            sel1_info = parse_selection_info(row['sel1'])
            sel2_info = parse_selection_info(row['sel2'])

            if not sel1_info or not sel2_info:
                continue

            atoms1 = get_atoms_from_selection(universe, sel1_info)
            atoms2 = get_atoms_from_selection(universe, sel2_info)

            if len(atoms1) == 0 or len(atoms2) == 0:
                continue

            # Encontrar átomos más cercanos
            min_distance = float('inf')
            closest_pair = None

            for atom1 in atoms1:
                for atom2 in atoms2:
                    distance = np.linalg.norm(atom1.position - atom2.position)
                    if distance < min_distance:
                        min_distance = distance
                        closest_pair = (atom1, atom2)

            if closest_pair and min_distance < 15.0:
                atom1, atom2 = closest_pair
                color = all_interactions_colors.get(row['interaction_name'],
                                                    '#888888')

                # Agregar átomos a la lista de participantes
                interaction_atoms.add(atom1.index)
                interaction_atoms.add(atom2.index)

                # Añadir línea discontinua
                add_dashed_line_segments(
                    viewer, atom1.position, atom2.position,
                    color, line_width, segments=8,
                    interaction_name=row['interaction_name'],
                    prevalence=row.get('prevalence', 0.0)
                )

                # Añadir esferas
                add_interaction_spheres(
                    viewer, atom1, atom2, color, sphere_size,
                    row['interaction_name'], row.get('prevalence', 0.0)
                )

                added_count += 1

        except Exception as e:
            print(f"⚠️ Error procesando interacción {idx}: {e}")
            continue

    print(f"✅ {added_count} interacciones procesadas")
    return interaction_atoms

def add_dashed_line_segments(viewer, pos1, pos2, color, line_width,
                             segments=8, interaction_name=None, prevalence=0.0):
    """
    Crear línea discontinua con datos de interacción para hover
    """
    try:
        start = np.array(pos1)
        end = np.array(pos2)
        direction = end - start
        total_length = np.linalg.norm(direction)

        if total_length == 0:
            return

        direction = direction / total_length
        segment_length = total_length / segments
        dash_length = segment_length * 0.6

        interaction_id = f"interaction_{abs(hash(f'{interaction_name}_{pos1[0]}_{pos2[0]}'))}"

        for i in range(segments):
            if i % 2 == 0:
                segment_start = start + direction * (i * segment_length)
                segment_end = start + direction * (i * segment_length + dash_length)

                viewer.addCylinder({
                    'start': {
                        'x': float(segment_start[0]),
                        'y': float(segment_start[1]),
                        'z': float(segment_start[2])
                    },
                    'end': {
                        'x': float(segment_end[0]),
                        'y': float(segment_end[1]),
                        'z': float(segment_end[2])
                    },
                    'radius': line_width,
                    'color': color,
                    'alpha': 0.8,
                    'interaction_name': interaction_name,
                    'prevalence': prevalence,
                    'interaction_id': interaction_id,
                    'is_interaction': True
                })

    except Exception as e:
        print(f"Error creando línea discontinua: {e}")

def add_interaction_spheres(viewer, atom1, atom2, color, sphere_size, interaction_name, prevalence=0.0):
    """
    Añadir esferas de interacción simplificadas
    """
    try:
        interaction_id = f"interaction_{abs(hash(f'{atom1.index}_{atom2.index}'))}"

        viewer.addSphere({
            'center': {
                'x': float(atom1.position[0]),
                'y': float(atom1.position[1]),
                'z': float(atom1.position[2])
            },
            'radius': sphere_size,
            'color': color,
            'alpha': 0.8,
            'interaction_name': interaction_name,
            'prevalence': prevalence,
            'interaction_id': interaction_id,
            'is_interaction': True
        })

        viewer.addSphere({
            'center': {
                'x': float(atom2.position[0]),
                'y': float(atom2.position[1]),
                'z': float(atom2.position[2])
            },
            'radius': sphere_size,
            'color': color,
            'alpha': 0.8,
            'interaction_name': interaction_name,
            'prevalence': prevalence,
            'interaction_id': interaction_id,
            'is_interaction': True
        })

    except Exception as e:
        print(f"Error añadiendo esferas: {e}")

def add_hover_labels(viewer, universe):
    """
    Hover labels simplificados y eficientes
    """
    try:
        viewer.setHoverable({}, True, """
        function(atom, viewer, event, container) {
            if (!atom) return;

            // Limpiar tooltips existentes
            const existing = container.querySelectorAll('.hover-tooltip');
            existing.forEach(tip => tip.remove());

            let info = '';
            let backgroundColor = 'rgba(0,0,0,0.85)';

            // Verificar si es una interacción
            if (atom.is_interaction && atom.interaction_name) {
                info = `🔗 ${atom.interaction_name} (${atom.prevalence}%)`;
                backgroundColor = 'rgba(64, 81, 181, 0.9)';
            } else {
                info = `${atom.atom || 'N/A'} - ${atom.resn || 'N/A'}${atom.resi || ''}`;
                backgroundColor = 'rgba(0,0,0,0.85)';
            }

            // Crear tooltip simple
            const tooltip = document.createElement('div');
            tooltip.className = 'hover-tooltip';
            tooltip.innerHTML = info;
            tooltip.style.cssText = `
                position: absolute;
                background: ${backgroundColor};
                color: white;
                padding: 6px 10px;
                border-radius: 4px;
                font-size: 11px;
                font-family: 'Roboto', sans-serif;
                z-index: 10000;
                pointer-events: none;
                box-shadow: 0 2px 6px rgba(0,0,0,0.3);
                white-space: nowrap;
            `;

            // Posicionar tooltip
            const rect = container.getBoundingClientRect();
            tooltip.style.left = (event.clientX - rect.left + 10) + 'px';
            tooltip.style.top = (event.clientY - rect.top - 5) + 'px';

            container.appendChild(tooltip);

            // Remover después de 2 segundos
            setTimeout(() => {
                if (tooltip.parentNode) {
                    tooltip.parentNode.removeChild(tooltip);
                }
            }, 2000);
        }
        """)

        print("Hover labels simplificados configurados")

    except Exception as e:
        print(f"Error configurando hover labels: {e}")

def parse_selection_info(selection_string):
    """
    Parsear información de selección desde el string

    Args:
        selection_string: String como 'TYR_92_234_CG_4473'

    Returns:
        dict: Información parseada o None si hay error
    """
    try:
        parts = selection_string.split('_')
        if len(parts) >= 3:
            return {
                'resname': parts[0],
                'resid': int(parts[1]),
                'atom_name': parts[3] if len(parts) > 3 else None,
                'atom_id': int(parts[4]) if len(parts) > 4 else None
            }
        return None
    except Exception as e:
        print(f"Error parseando selección {selection_string}: {e}")
        return None

def get_atoms_from_selection(universe, sel_info):
    """
    Obtener átomos de MDAnalysis basado en información de selección

    Args:
        universe: Universo MDAnalysis
        sel_info: Información de selección parseada

    Returns:
        list: Lista de átomos
    """
    try:
        selection = f"resname {sel_info['resname']} and resid {sel_info['resid']}"

        if sel_info['atom_name']:
            selection += f" and name {sel_info['atom_name']}"

        atoms = universe.select_atoms(selection)
        return atoms

    except Exception as e:
        print(f"Error obteniendo átomos: {e}")
        return []

def filter_interactions_by_frame(df, frame_index):
    """
    Filtrar interacciones que ocurren en un frame específico
    """
    filtered_rows = []

    for idx, row in df.iterrows():
        timeseries = row['timeseries']
        if isinstance(timeseries, str) and len(timeseries) > frame_index:
            if timeseries[frame_index] == '1':
                filtered_rows.append(row)

    return pd.DataFrame(filtered_rows) if filtered_rows else pd.DataFrame()


def apply_atom_transparency(viewer, interaction_atoms,
                            molecule_style="cartoon"):
    """
    Aplicar transparencia respetando el estilo de visualización seleccionado
    """
    try:
        if molecule_style == "cartoon":
            transparent_style = {'cartoon': {'opacity': 0.4}}
            opaque_style = {'cartoon': {'opacity': 1.0}}
        elif molecule_style == "stick":
            transparent_style = {'stick': {'opacity': 0.4}}
            opaque_style = {'stick': {'opacity': 1.0}}
        elif molecule_style == "sphere":
            transparent_style = {'sphere': {'opacity': 0.4}}
            opaque_style = {'sphere': {'opacity': 1.0}}
        elif molecule_style == "ball_and_stick":
            transparent_style = {
                'stick': {'opacity': 0.4},
                'sphere': {'opacity': 0.4}
            }
            opaque_style = {
                'stick': {'opacity': 1.0},
                'sphere': {'opacity': 1.0}
            }
        elif molecule_style == "cartoon_and_stick":
            transparent_style = {
                'cartoon': {'opacity': 0.4},
                'stick': {'opacity': 0.4}
            }
            opaque_style = {
                'cartoon': {'opacity': 1.0},
                'stick': {'opacity': 1.0}
            }
        else:
            transparent_style = {'cartoon': {'opacity': 0.4}}
            opaque_style = {'cartoon': {'opacity': 1.0}}

        viewer.setStyle({}, transparent_style)

        for atom_index in interaction_atoms:
            viewer.setStyle({'index': atom_index}, opaque_style)


    except Exception as e:
        print(f"Error aplicando transparencia: {e}")
