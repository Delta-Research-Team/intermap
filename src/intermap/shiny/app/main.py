"""
Main application file for InterMap Visualizations.
Integrates UI components with server-side logic.
"""

import os
import tempfile

from shiny import App, reactive, render, ui

from .css import ERROR_MESSAGES
from .icsv import CSVFilter, traspose, sortby
from .plots import (create_interactions_over_time_plot,
                    create_ligand_interactions_plot, create_plot,
                    create_receptor_interactions_plot)
from .ui import app_ui


def server(input, output, session):
    """Define server logic for the app."""

    # =========================================================================
    # Reactive States
    # =========================================================================
    csv = reactive.Value(None)
    filtered_idx = reactive.Value(None)
    csv_filtered = reactive.Value(None)

    @reactive.Effect
    @reactive.event(input.csv_file, input.top_file)
    def initialize_filter():
        """Initialize CSVFilter when both files are uploaded."""
        try:
            # Verificar que ambos archivos estén presentes
            csv_infos = input.csv_file()
            top_infos = input.top_file()

            if not csv_infos or not top_infos:
                csv.set(None)
                filtered_idx.set(None)
                csv_filtered.set(None)
                return

            # Crear directorio temporal
            temp_dir = tempfile.mkdtemp()

            # Guardar ambos archivos en el directorio temporal
            csv_path = os.path.join(temp_dir, csv_infos[0]["name"])
            top_path = os.path.join(temp_dir, top_infos[0]["name"])

            # Copiar los archivos usando shutil
            import shutil
            shutil.copy2(csv_infos[0]["datapath"], csv_path)
            shutil.copy2(top_infos[0]["datapath"], top_path)

            # Inicializar CSVFilter con ambos archivos
            master_instance = CSVFilter(
                csv=csv_path,
                topo=top_path
            )

            csv.set(master_instance)

            ui.notification_show(
                "Files loaded successfully!",
                duration=3000,
                type="message"
            )

            # Seleccionar automáticamente la primera pestaña
            ui.update_navs("plot_tabs", selected="Sele1 vs Sele2")

        except Exception as e:
            print(f"Error in initialize_filter: {str(e)}")
            csv.set(None)
            filtered_idx.set(None)
            csv_filtered.set(None)
            ui.notification_show(
                f"Error loading files: {str(e)}",
                duration=5000,
                type="error"
            )

    @reactive.Effect
    @reactive.event(input.plot_button,
                    input.selected_interactions,
                    input.selected_annotations,
                    input.prevalence_threshold,
                    input.sort_by)
    def update_filtered_idx():
        """Update filtered indices based on current UI selections."""
        if csv.get() is None:
            filtered_idx.set(None)
            csv_filtered.set(None)
            return

        master_instance = csv.get()

        # Obtain idx's using CSVFilter
        mda_idx, mda_status = master_instance.by_mda(
            input.mda_selection() if input.mda_selection() else 'all')
        prev_idx, prev_status = master_instance.by_prevalence(
            input.prevalence_threshold())
        inter_idx, inter_status = master_instance.by_inters(tuple(
            input.selected_interactions()) if input.selected_interactions() else input.selected_interactions())
        annot_idx, notes_status = master_instance.by_notes(tuple(
            input.selected_annotations()) if input.selected_annotations() else input.selected_annotations())

        # Status section
        if mda_status == -1:
            ui.notification_show(
                "There is not information for this MDAnalysis selection"
                , type="message", duration=5)
        elif prev_status == -1:
            ui.notification_show(
                "There are no interactions with this prevalence or higher."
                , type="message", duration=5)
        elif inter_status == -1:
            ui.notification_show(
                "There is not information for this Interaction selection"
                , type="message", duration=5)
        elif notes_status == -1:
            ui.notification_show(
                "There is not information for this Annotation selection"
                , type="message", duration=5)

        # Update reactive states
        df_idx = set.intersection(mda_idx, prev_idx, inter_idx, annot_idx)
        if df_idx:
            filtered_df = csv.get().master.iloc[list(df_idx)]

            # Apply transpose if checkbox is checked
            if input.transpose_data():
                filtered_df = traspose(filtered_df)

            # Apply sorting based on radio button selection
            filtered_df = sortby(filtered_df, input.sort_by())

            filtered_idx.set(df_idx)
            csv_filtered.set(filtered_df)
        else:
            filtered_idx.set(None)
            csv_filtered.set(None)

    # =========================================================================
    # Rendering helper functions
    # =========================================================================
    def check_data():
        """Helper function to check if filtered data is available."""
        if csv_filtered.get() is None:
            return ui.p(ERROR_MESSAGES["nothing_selected"])
        return None

    def render_plot(create_plot_func, is_main=False):
        """Helper function to render plots with common logic."""
        error = check_data()
        if error is not None:
            return error

        plot_args = [
            csv_filtered.get(),
            input.plot_width(),
            input.plot_height() if is_main else input.plot_height() // 1.5,
        ]

        if is_main:
            plot_args.append(input.show_prevalence())

        fig = create_plot_func(*plot_args)

        if fig is None:
            return ui.p(ERROR_MESSAGES["nothing_selected"])

        plot_html = fig.to_html(
            include_plotlyjs="cdn",
            full_html=False,
            config={"responsive": True}
        )

        style = (
            "width:100%; height:800px; border:1px solid white; margin: 10px 0;"
            if is_main else
            f"width:100%; height:{input.plot_height() // 2}px;")

        return ui.tags.div(ui.HTML(plot_html), style=style)

    def create_checkbox_group(input_id, data_keys_getter):
        """Helper function to create checkbox groups."""
        if csv.get() is None:
            return ui.input_checkbox_group(input_id, "", [])

        choices = list(data_keys_getter(csv.get()))
        return ui.input_checkbox_group(input_id, "", choices, selected=choices)

    # =========================================================================
    # Plot Rendering Functions
    # =========================================================================
    @output
    @render.ui
    @reactive.event(input.plot_button)
    def interaction_plot():
        """Render the main interaction heatmap plot."""
        return render_plot(create_plot, is_main=True)

    @output
    @render.ui
    @reactive.event(input.plot_button)
    def ligand_interactions_plot():
        """Render the ligand interactions plot."""
        return render_plot(create_ligand_interactions_plot)

    @output
    @render.ui
    @reactive.event(input.plot_button)
    def receptor_interactions_plot():
        """Render the receptor interactions plot."""
        return render_plot(create_receptor_interactions_plot)

    @output
    @render.ui
    @reactive.event(input.plot_button)
    def interactions_over_time_plot():
        """Render the interactions over time plot."""
        return render_plot(create_interactions_over_time_plot)

    # =========================================================================
    # Checkbox Rendering
    # =========================================================================
    @output
    @render.ui
    def interaction_checkboxes():
        """Render the interaction checkboxes."""
        return create_checkbox_group("selected_interactions",
                                     lambda x: x.inters2df.keys())

    @output
    @render.ui
    def annotation_checkboxes():
        """Render the annotation checkboxes."""
        return create_checkbox_group("selected_annotations",
                                     lambda x: x.notes2df.keys())

    @reactive.Effect
    @reactive.event(input.select_all_interactions)
    def update_interaction_selections():
        """Handle select all/deselect all for interactions."""
        if csv.get() is None:
            return

        choices = list(csv.get().inters2df.keys())
        if input.select_all_interactions():
            ui.update_checkbox_group("selected_interactions", selected=choices)
        else:
            ui.update_checkbox_group("selected_interactions", selected=[])

    @reactive.Effect
    @reactive.event(input.select_all_annotations)
    def update_annotation_selections():
        """Handle select all/deselect all for annotations."""
        if csv.get() is None:
            return

        choices = list(csv.get().notes2df.keys())
        if input.select_all_annotations():
            ui.update_checkbox_group("selected_annotations", selected=choices)
        else:
            ui.update_checkbox_group("selected_annotations", selected=[])

app = App(app_ui, server)
