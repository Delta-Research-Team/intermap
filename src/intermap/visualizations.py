import pandas as pd
from shiny import App, render, ui, reactive
import threading
import webbrowser

import pandas as pd
import plotly.graph_objects as go
from shiny import App, reactive, render, ui

all_interactions_colors = {
    'HBDonor': '#E41A1C',
    'HBAcceptor': '#377EB8',
    'Hydrophobic': '#4DAF4A',
    'VdWContact': '#FF7F00',
    'CloseContact': '#984EA3',
    'PiStacking': '#F781BF',
    'PiCation': '#A65628',
    'CationPi': '#999999',
    'Anionic': '#B2182B',
    'Cationic': '#2166AC',
    'MetalDonor': '#A8A8A8',
    'MetalAcceptor': '#8C510A',
    'FaceToFace': '#762A83',
    'EdgeToFace': '#1B7837',
    'XBAcceptor': '#5E3C99',
    'XBDonor': '#2B5C8C'
}


def get_image_base64(image_path):
    """Convierte una imagen a base64."""
    try:
        with open(image_path, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode()
            return f"data:image/png;base64,{encoded_string}"
    except Exception as e:
        print(f"Error al codificar la imagen: {e}")
        return None



def generate_interaction_choices(df):
    if df is None:
        return {k: f"{k}" for k in all_interactions_colors.keys()}

    counts = df['interaction'].value_counts()
    return {
        k: ui.span(
            k,
            ui.span(
                f" ({counts.get(k, 0)})",
                style="color: #666; font-size: 0.9em;"
            )
        ) for k in all_interactions_colors.keys()
    }


app_ui = ui.page_fluid(
    ui.tags.head(
        ui.tags.script(src="https://cdn.plot.ly/plotly-latest.min.js"),
        ui.tags.style("""
            .welcome-section {
                display: flex;
                justify-content: space-between;
                align-items: center;
                margin-bottom: 20px;
                padding: 20px;
                background: white;
                width: 100%;
            }
            .welcome-text {
                text-align: left;
                padding-left: 20px;
            }
            .welcome-title {
                font-size: 36px;
                font-weight: bold;
                margin: 0;
                color: #333;
            }
            .welcome-subtitle {
                font-size: 18px;
                color: #666;
                margin-top: 5px;
            }
            .welcome-image {
                max-width: 200px;
                height: auto;
                margin-right: 20px;
            }
            .interaction-filter {
                background: #f8f9fa;
                padding: 15px;
                border-radius: 8px;
                margin-bottom: 15px;
            }
            .checkbox-group {
                max-height: 300px;
                overflow-y: auto;
                padding: 10px;
                border: 1px solid #ddd;
                border-radius: 4px;
            }
            .interaction-count {
                color: #666;
                font-size: 0.9em;
                margin-left: 4px;
            }
            .plot-container {
                margin-bottom: 20px;
                border: 1px solid #ddd;
                border-radius: 4px;
                padding: 10px;
            }
            .search-container {
                display: flex;
                align-items: center;
                gap: 10px;
                margin-bottom: 10px;
            }
            .search-button {
                background-color: #007bff;
                color: white;
                border: none;
                padding: 6px 12px;
                border-radius: 4px;
                cursor: pointer;
            }
            .search-button:hover {
                background-color: #0056b3;
            }
            .not-found {
                color: red;
                margin-top: 5px;
                font-size: 0.9em;
            }
        """)
    ),
    ui.div(
        {"class": "welcome-section"},
        ui.div(
            {"class": "welcome-text"},
            ui.h1("Welcome to InterMap Visualizations!", {"class": "welcome-title"}),
            ui.p("Revolutionizing Molecular Interaction Analysis with High-Speed Computation",
                 {"class": "welcome-subtitle", "style": "font-style: italic;"})
        ),
        ui.img(
            {"src": get_image_base64(
                "/home/fajardo/03_Fajardo_Hub/02_InterMap/Last_version/intermap/docs/figs/imap.png"),
                "class": "welcome-image",
                "alt": "InterMap Logo"}
        )
    ),
    ui.row(
        ui.column(3,
                  ui.div(
                      {"class": "interaction-filter"},
                      ui.h4("Data Input"),
                      ui.input_file("file", "Upload CSV file:",
                                    accept=[".csv"]),
                      ui.hr(),

                      ui.h4("Permanence Filter"),
                      ui.input_slider(
                          "permanence_threshold",
                          "Filter by Permanence (%)",
                          min=0,
                          max=100,
                          value=0,
                          step=1
                      ),
                      ui.hr(),
                      ui.input_switch(
                          "show_permanence",
                          "Show Permanence Values",
                          value=False
                      ),
                      ui.hr(),
                      ui.h4("Residue Filter"),
                      ui.div(
                          {"class": "search-container"},
                          ui.input_text(
                              "residue_filter",
                              "",
                              placeholder="Enter residue name...",
                              width="100%"
                          ),
                          ui.input_action_button(
                              "search_button",
                              "Find",
                              class_="search-button"
                          )
                      ),
                      ui.output_ui("residue_not_found"),
                      ui.hr(),
                      ui.h4("Interaction Filter"),
                      ui.div(
                          {"class": "checkbox-group"},
                          ui.output_ui("interaction_checkboxes")
                      ),
                      ui.hr(),
                      ui.h4("Plot Settings"),
                      ui.input_numeric("plot_width", "Plot Width:", value=1000,
                                       min=400, max=2000),
                      ui.input_numeric("plot_height", "Plot Height:",
                                       value=800, min=300, max=1500),
                      ui.hr(),

                  )
                  ),
        ui.column(9,
                  ui.div(
                      {"class": "plot-container"},
                      ui.h4("Interaction Heatmap"),
                      ui.output_ui("interaction_plot")
                  ),
                  ui.div(
                      {"class": "plot-container"},
                      ui.h4("Ligand Atoms Interaction Analysis"),
                      ui.output_ui("ligand_interactions_plot")
                  ),
                  ui.div(
                      {"class": "plot-container"},
                      ui.h4("Receptor Atoms Interaction Analysis"),
                      ui.output_ui("receptor_interactions_plot")
                  ),
                  ui.div(
                      {"class": "plot-container"},
                      ui.h4("Interaction Types Over Time"),
                      ui.output_ui("interactions_over_time_plot")
                  ),
                  )
    )
)


def server(input, output, session):
    data_store = reactive.Value({'df': None})
    search_state = reactive.Value({'search_term': '', 'active': False})

    @reactive.Effect
    @reactive.event(input.search_button)
    def _():
        search_state.set({
            'search_term': input.residue_filter(),
            'active': True
        })

    def calculate_permanence(row):
        frame_data = row[3:]
        frame_data = frame_data[frame_data.notna()]
        total_frames = len(frame_data)
        frames_with_interaction = sum(frame_data == 1)
        return round((frames_with_interaction / total_frames) * 100,
                     1) if total_frames > 0 else 0

    def _adjust_color_transparency(color, opacity):
        if color.startswith('#'):
            r = int(color[1:3], 16)
            g = int(color[3:5], 16)
            b = int(color[5:7], 16)
            return f'rgba({r},{g},{b},{opacity})'
        elif color.startswith('rgb'):
            return color.replace('rgb', 'rgba').replace(')', f',{opacity})')
        return color

    def create_plot(df, width, height, selected_interactions, permanence_threshold, show_permanence=False):
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

        df_filtered = df[df['interaction'].isin(selected_interactions)]
        df_filtered['permanence'] = df_filtered.apply(calculate_permanence, axis=1)

        if search_state.get()['active'] and search_state.get()['search_term']:
            search_term = search_state.get()['search_term'].upper()
            mask = (df_filtered['ligand'].str.contains(search_term, case=False) |
                    df_filtered['protein'].str.contains(search_term, case=False))
            df_filtered = df_filtered[mask]

        if df_filtered.empty:
            return None

        df_filtered = df_filtered[df_filtered['permanence'] >= permanence_threshold]

        priority_interactions = []

        for (lig, prot), group in df_filtered.groupby(['ligand', 'protein']):
            if not group.empty:
                group['priority'] = group['interaction'].map(interaction_priority)
                group = group.sort_values(['priority', 'permanence'],
                                          ascending=[True, False])

                priority_int = group.iloc[0]
                priority_interactions.append({
                    'ligand': lig,
                    'protein': prot,
                    'interaction': priority_int['interaction'],
                    'permanence': priority_int['permanence']
                })

        priority_df = pd.DataFrame(priority_interactions)

        pivot_interaction = pd.pivot_table(priority_df,
                                           values='interaction',
                                           index='ligand',
                                           columns='protein',
                                           aggfunc='first',
                                           fill_value='')

        pivot_permanence = pd.pivot_table(priority_df,
                                          values='permanence',
                                          index='ligand',
                                          columns='protein',
                                          aggfunc='first',
                                          fill_value=0)

        fig = go.Figure()

        fig.add_trace(go.Heatmap(
            z=[[0] * len(pivot_interaction.columns)] * len(pivot_interaction.index),
            x=pivot_interaction.columns,
            y=pivot_interaction.index,
            showscale=False,
            colorscale=[[0, '#FEFBF6'], [1, '#FEFBF6']],
            hoverongaps=False,
            hoverinfo='skip',
        ))

        present_interactions = sorted(priority_df['interaction'].unique(),
                                      key=lambda x: interaction_priority[x])

        color_scale = []
        for i, interaction in enumerate(present_interactions):
            pos = i / (len(present_interactions) - 1) if len(present_interactions) > 1 else 0.5
            color_scale.extend([[pos, all_interactions_colors[interaction]],
                                [pos, all_interactions_colors[interaction]]])

        if len(present_interactions) == 1:
            color = all_interactions_colors[present_interactions[0]]
            color_scale = [[0, color], [1, color]]

        interaction_to_num = {inter: i for i, inter in enumerate(present_interactions)}
        pivot_numerical = pivot_interaction.replace(interaction_to_num)

        text_values = pivot_permanence.round(1).astype(str)
        if show_permanence:
            text_values = text_values.mask(pivot_permanence == 0, '')
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
            hovertemplate="Ligand: %{y}<br>Protein: %{x}<br>" +
                          "Interaction: %{customdata}<br>" +
                          "Permanence: %{text}%<extra></extra>",
            customdata=pivot_interaction.values
        ))

        for interaction in present_interactions:
            fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                marker=dict(
                    size=10,
                    color=all_interactions_colors[interaction],
                    symbol='square',
                ),
                name=f"{interaction} ({interaction_priority[interaction]})",
                showlegend=True
            ))

        fig.update_layout(
            width=width,
            height=height,
            margin=dict(l=50, r=150, t=80, b=50),
            showlegend=True,
            legend=dict(
                title="Interaction Types (Priority)",
                yanchor="top",
                y=0.99,
                xanchor="left",
                x=1.02,
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='rgba(0,0,0,0.2)',
                borderwidth=1
            ),
            xaxis_title="Protein Residues",
            yaxis_title="Ligand Atoms",
            annotations=[
                dict(
                    text="Generated by: LQC3",
                    xref="paper",
                    yref="paper",
                    x=1.02,
                    y=-0.1,
                    showarrow=False,
                    font=dict(size=8, color="gray")
                )
            ],
            autosize=False,
            paper_bgcolor='white',
            plot_bgcolor='white'
        )

        return fig

    def create_ligand_interactions_plot(df, width, height, selected_interactions, permanence_threshold):
        if df is None or df.empty:
            return None

        df_filtered = df[df['interaction'].isin(selected_interactions)]
        df_filtered['permanence'] = df_filtered.apply(calculate_permanence, axis=1)

        if search_state.get()['active'] and search_state.get()['search_term']:
            search_term = search_state.get()['search_term'].upper()
            if df_filtered['ligand'].str.contains(search_term, case=False).any():
                df_filtered = df_filtered[
                    df_filtered['ligand'].str.contains(search_term, case=False)]
            else:
                return None


        df_filtered = df_filtered[df_filtered['permanence'] >= permanence_threshold]

        if df_filtered.empty:
            return None


        ligand_interactions = df_filtered.groupby(['ligand', 'interaction'])[
            'permanence'].mean().reset_index()

        fig = go.Figure()
        legend_entries = set()


        for ligand in sorted(ligand_interactions['ligand'].unique()):
            ligand_data = ligand_interactions[ligand_interactions['ligand'] == ligand]

            sorted_interactions = ligand_data.sort_values('permanence', ascending=True)

            prev_height = 0

            for _, row in sorted_interactions.iterrows():
                if row['interaction'] in selected_interactions:
                    show_legend = row['interaction'] not in legend_entries
                    legend_entries.add(row['interaction'])

                    current_height = row['permanence']
                    bar_height = current_height - prev_height

                    fig.add_trace(go.Bar(
                        name=row['interaction'],
                        x=[ligand],
                        y=[bar_height],
                        marker_color=all_interactions_colors[row['interaction']],
                        marker=dict(
                            color=all_interactions_colors[row['interaction']],
                            line=dict(
                                color=all_interactions_colors[row['interaction']],
                                width=1
                            )
                        ),
                        opacity=1.0,
                        showlegend=show_legend,
                        base=prev_height,
                        hovertemplate=(
                                "Ligand: %{x}<br>" +
                                "Interaction: " + row['interaction'] + "<br>" +
                                "permanence: " + f"{current_height:.1f}%" +
                                "<extra></extra>"
                        )
                    ))

                    prev_height = current_height

        fig.update_layout(
            width=width,
            height=height,
            barmode='overlay',
            xaxis_title="Ligand Atoms",
            yaxis_title="Interaction Permanence (%)",
            showlegend=True,
            legend=dict(
                title="Interaction Types",
                yanchor="top",
                y=0.99,
                xanchor="left",
                x=1.02,
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='rgba(0,0,0,0.2)',
                borderwidth=1
            ),
            margin=dict(l=50, r=150, t=50, b=50),
            paper_bgcolor='white',
            plot_bgcolor='white',
            legend_traceorder='normal',
            bargap=0.15,
            bargroupgap=0.1
        )

        fig.update_xaxes(tickangle=45)
        fig.update_yaxes(range=[0, 100])

        return fig

    def create_receptor_interactions_plot(df, width, height,
                                          selected_interactions,
                                          permanence_threshold):
        if df is None or df.empty:
            return None

        df_filtered = df[df['interaction'].isin(selected_interactions)]
        df_filtered['permanence'] = df_filtered.apply(calculate_permanence,
                                                      axis=1)

        if search_state.get()['active'] and search_state.get()['search_term']:
            search_term = search_state.get()['search_term'].upper()
            if df_filtered['protein'].str.contains(search_term,
                                                   case=False).any():
                df_filtered = df_filtered[
                    df_filtered['protein'].str.contains(search_term,
                                                        case=False)]
            else:
                return None

        df_filtered = df_filtered[
            df_filtered['permanence'] >= permanence_threshold]

        if df_filtered.empty:
            return None

        receptor_interactions = \
        df_filtered.groupby(['protein', 'interaction'])[
            'permanence'].mean().reset_index()

        fig = go.Figure()
        legend_entries = set()

        for residue in receptor_interactions['protein'].unique():
            residue_data = receptor_interactions[
                receptor_interactions['protein'] == residue]
            sorted_interactions = residue_data.sort_values('permanence',
                                                           ascending=False)

            for _, row in sorted_interactions.iterrows():
                if row['interaction'] in selected_interactions:
                    show_legend = row['interaction'] not in legend_entries
                    legend_entries.add(row['interaction'])

                    fig.add_trace(go.Bar(
                        name=row['interaction'],
                        x=[residue],
                        y=[row['permanence']],
                        marker_color=all_interactions_colors[
                            row['interaction']],
                        marker=dict(
                            color=all_interactions_colors[row['interaction']],
                            line=dict(
                                color=all_interactions_colors[
                                    row['interaction']],
                                width=1
                            )
                        ),
                        opacity=1.0,
                        showlegend=show_legend
                    ))

        fig.update_layout(
            width=width,
            height=height,
            barmode='overlay',
            xaxis_title="Receptor Residues",
            yaxis_title="Interaction Permanence (%)",
            showlegend=True,
            legend=dict(
                title="Interaction Types",
                yanchor="top",
                y=0.99,
                xanchor="left",
                x=1.02,
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='rgba(0,0,0,0.2)',
                borderwidth=1
            ),
            margin=dict(l=50, r=150, t=50, b=50),
            paper_bgcolor='white',
            plot_bgcolor='white',
            legend_traceorder='normal',
            bargap=0.15,
            bargroupgap=0.1
        )

        fig.update_xaxes(tickangle=45)
        fig.update_yaxes(range=[0, 100])

        return fig

    def create_interactions_over_time_plot(df, width, height,
                                           selected_interactions,
                                           permanence_threshold):
        if df is None or df.empty:
            return None

        df_filtered = df[df['interaction'].isin(selected_interactions)]

        if search_state.get()['active'] and search_state.get()['search_term']:
            search_term = search_state.get()['search_term'].upper()
            mask = (df_filtered['ligand'].str.contains(search_term,
                                                       case=False) |
                    df_filtered['protein'].str.contains(search_term,
                                                        case=False))
            df_filtered = df_filtered[mask]

        if df_filtered.empty:
            return None

        frame_columns = df_filtered.columns[3:]
        interaction_counts = {interaction: [] for interaction in
                              selected_interactions}
        frame_numbers = list(range(1, len(frame_columns) + 1))
        present_interactions = set()

        for frame in frame_columns:
            frame_data = df_filtered[df_filtered[frame] == 1]
            counts = frame_data['interaction'].value_counts()

            for interaction in selected_interactions:
                count = counts.get(interaction, 0)
                interaction_counts[interaction].append(count)
                if count > 0:
                    present_interactions.add(interaction)

        fig = go.Figure()

        for interaction in present_interactions:
            fig.add_trace(go.Scatter(
                x=frame_numbers,
                y=interaction_counts[interaction],
                mode='markers',
                name=interaction,
                marker=dict(
                    color=all_interactions_colors[interaction],
                    size=8,
                    symbol='circle',
                    line=dict(
                        color='rgba(255,255,255,0.5)',
                        width=1
                    )
                ),
                hovertemplate="Frame: %{x}<br>" +
                              "Count: %{y}<br>" +
                              f"Type: {interaction}<extra></extra>"
            ))

        fig.update_layout(
            width=width,
            height=height,
            xaxis_title="Frame Number",
            yaxis_title="Number of Interactions",
            showlegend=True,
            legend=dict(
                title="Interaction Types",
                yanchor="top",
                y=0.99,
                xanchor="left",
                x=1.02,
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='rgba(0,0,0,0.2)',
                borderwidth=1
            ),
            margin=dict(l=50, r=150, t=50, b=50),
            paper_bgcolor='white',
            plot_bgcolor='white',
            hovermode='x unified'
        )

        fig.update_xaxes(
            gridcolor='rgba(0,0,0,0.1)',
            zeroline=True,
            zerolinecolor='rgba(0,0,0,0.2)'
        )
        fig.update_yaxes(
            gridcolor='rgba(0,0,0,0.1)',
            zeroline=True,
            zerolinecolor='rgba(0,0,0,0.2)'
        )

        return fig

    @output
    @render.ui
    def residue_not_found():
        if not search_state.get()['active'] or data_store.get()['df'] is None:
            return ui.div()

        df = data_store.get()['df']
        search_term = search_state.get()['search_term'].upper()

        if not search_term:
            return ui.div()

        ligand_matches = df['ligand'].str.contains(search_term,
                                                   case=False).any()
        protein_matches = df['protein'].str.contains(search_term,
                                                     case=False).any()

        if not ligand_matches and not protein_matches:
            return ui.p("Not Found", style="color: red; margin-top: 5px;")
        return ui.div()

    @output
    @render.ui
    @reactive.event(input.file, input.plot_width, input.plot_height,
                    input.selected_interactions, input.permanence_threshold,
                    input.show_permanence, input.search_button)
    def interaction_plot():
        file_infos = input.file()
        if not file_infos:
            return ui.p("Please upload a CSV file.")

        file_info = file_infos[0]
        try:
            if data_store.get()['df'] is None:
                df = pd.read_csv(file_info['datapath'])
                data_store.set({'df': df})
            else:
                df = data_store.get()['df']

            fig = create_plot(
                df,
                input.plot_width(),
                input.plot_height(),
                input.selected_interactions(),
                input.permanence_threshold(),
                input.show_permanence()
            )

            if fig is None:
                return ui.p("No interactions found with the current filters.")

            plot_html = fig.to_html(
                include_plotlyjs='cdn',
                full_html=False,
                config={'responsive': True}
            )

            return ui.tags.div(
                ui.HTML(plot_html),
                style="width:100%; height:800px; border:1px solid #ddd; margin: 10px 0;"
            )

        except Exception as e:
            import traceback
            print(f"Error in plot generation:")
            print(traceback.format_exc())
            return ui.div(
                ui.p(f"Error: {str(e)}", style="color:red;"),
                ui.pre(traceback.format_exc())
            )



    @output
    @render.ui
    @reactive.event(input.file, input.plot_width, input.plot_height,
                    input.selected_interactions, input.permanence_threshold,
                    input.search_button)
    def ligand_interactions_plot():
        data = data_store.get()
        if data['df'] is None:
            return ui.p("Please upload a CSV file.")

        try:
            fig = create_ligand_interactions_plot(
                data['df'],
                input.plot_width(),
                input.plot_height() // 2,
                input.selected_interactions(),
                input.permanence_threshold()
            )

            if fig is None:
                return ui.p(
                    "No ligand interactions found with the current filters.")

            plot_html = fig.to_html(
                include_plotlyjs='cdn',
                full_html=False,
                config={'responsive': True}
            )

            return ui.tags.div(
                ui.HTML(plot_html),
                style=f"width:100%; height:{input.plot_height() // 2}px;"
            )

        except Exception as e:
            import traceback
            print(f"Error in plot generation:")
            print(traceback.format_exc())
            return ui.div(
                ui.p(f"Error: {str(e)}", style="color:red;"),
                ui.pre(traceback.format_exc())
            )

    @output
    @render.ui
    @reactive.event(input.file, input.plot_width, input.plot_height,
                    input.selected_interactions, input.permanence_threshold,
                    input.search_button)
    def receptor_interactions_plot():
        data = data_store.get()
        if data['df'] is None:
            return ui.p("Please upload a CSV file.")

        try:
            fig = create_receptor_interactions_plot(
                data['df'],
                input.plot_width(),
                input.plot_height() // 2,
                input.selected_interactions(),
                input.permanence_threshold()
            )

            if fig is None:
                return ui.p(
                    "No receptor interactions found with the current filters.")

            plot_html = fig.to_html(
                include_plotlyjs='cdn',
                full_html=False,
                config={'responsive': True}
            )

            return ui.tags.div(
                ui.HTML(plot_html),
                style=f"width:100%; height:{input.plot_height() // 2}px;"
            )

        except Exception as e:
            import traceback
            print(f"Error in plot generation:")
            print(traceback.format_exc())
            return ui.div(
                ui.p(f"Error: {str(e)}", style="color:red;"),
                ui.pre(traceback.format_exc())
            )

    @output
    @render.ui
    @reactive.event(input.file, input.plot_width, input.plot_height,
                    input.selected_interactions, input.permanence_threshold,
                    input.search_button)
    def interactions_over_time_plot():
        data = data_store.get()
        if data['df'] is None:
            return ui.p("Please upload a CSV file.")

        try:
            fig = create_interactions_over_time_plot(
                data['df'],
                input.plot_width(),
                input.plot_height() // 2,
                input.selected_interactions(),
                input.permanence_threshold()
            )

            if fig is None:
                return ui.p(
                    "No interactions over time found with the current filters.")

            plot_html = fig.to_html(
                include_plotlyjs='cdn',
                full_html=False,
                config={'responsive': True}
            )

            return ui.tags.div(
                ui.HTML(plot_html),
                style=f"width:100%; height:{input.plot_height() // 2}px;"
            )

        except Exception as e:
            import traceback
            print(f"Error in plot generation:")
            print(traceback.format_exc())
            return ui.div(
                ui.p(f"Error: {str(e)}", style="color:red;"),
                ui.pre(traceback.format_exc())
            )

    @output
    @render.ui
    def interaction_checkboxes():
        df = data_store.get()['df']
        choices = generate_interaction_choices(df)

        return ui.input_checkbox_group(
            "selected_interactions",
            "",
            choices=choices,
            selected=list(all_interactions_colors.keys())
        )

    @reactive.Effect
    @reactive.event(input.file)
    def _():
        file_infos = input.file()
        if file_infos and file_infos[0] is not None:
            df = pd.read_csv(file_infos[0]['datapath'])
            data_store.set({'df': df})

    @session.download(filename="interaction_map.png")
    async def download_png():
        data = data_store.get()
        if data['df'] is None:
            return None

        fig = create_plot(
            data['df'],
            input.plot_width(),
            input.plot_height(),
            input.selected_interactions(),
            input.permanence_threshold()
        )
        if fig is not None:
            return fig.to_image(format="png", scale=4)

    @session.download(filename="interaction_map.html")
    async def download_html():
        data = data_store.get()
        if data['df'] is None:
            return None

        fig = create_plot(
            data['df'],
            input.plot_width(),
            input.plot_height(),
            input.selected_interactions(),
            input.permanence_threshold()
        )
        if fig is not None:
            return fig.to_html()





app = App(app_ui, server)


def find_free_port(start_port=8001, max_port=8999):
    import socket
    for port in range(start_port, max_port + 1):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(('127.0.0.1', port))
                return port
        except OSError:
            continue
    raise OSError("No available port found")


def open_browser(port):
    url = f"http://127.0.0.1:{port}"
    webbrowser.open_new(url)


if __name__ == "__main__":
    import warnings

    warnings.filterwarnings('ignore')

    try:
        port = find_free_port()
        print(f"Starting server on port {port}")

        threading.Timer(1.25, open_browser, args=[port]).start()

        app.run(host="127.0.0.1", port=port)
    except Exception as e:
        print(f"Error starting the server: {e}")
