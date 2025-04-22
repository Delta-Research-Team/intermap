"""
Main application file for InterMap Visualizations.
Integrates UI components with server-side logic.
"""

from shiny import App, reactive, render, ui

from .css import ERROR_MESSAGES
from .icsv import CSVFilter
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
    @reactive.event(input.file)
    def initialize_filter():
        """Initialize CSVFilter when a file is uploaded."""
        try:
            file_infos = input.file()
            if not file_infos:
                csv.set(None)
                filtered_idx.set(None)
                csv_filtered.set(None)
                return

            csv_path = file_infos[0]["datapath"]
            master_instance = CSVFilter(csv_path)

            print(f"Topology path found: {master_instance.topo_path}")

            csv.set(master_instance)
            ui.notification_show("File loaded successfully!", type="message",
                                 duration=5)
        except Exception as e:
            print(f"Error in initialize_filter: {str(e)}")
            csv.set(None)
            filtered_idx.set(None)
            csv_filtered.set(None)

    @reactive.Effect
    @reactive.event(input.mda_selection_submit,
                    input.selected_interactions,
                    input.selected_annotations,
                    input.prevalence_threshold)
    def update_filtered_idx():
        """Update filtered indices based on current UI selections."""
        if csv.get() is None:
            filtered_idx.set(None)
            csv_filtered.set(None)
            return

        master_instance = csv.get()

        # Obtain idx's using CSVFilter
        mda_idx, _ = master_instance.by_mda(
            'all' if not input.mda_selection() else input.mda_selection())
        prev_idx, _ = master_instance.by_prevalence(
            input.prevalence_threshold())
        inter_idx, _ = master_instance.by_inters(tuple(
            input.selected_interactions()) if input.selected_interactions() else 'all')
        annot_idx, _ = master_instance.by_notes(tuple(
            input.selected_annotations()) if input.selected_annotations() else 'all')

        # Update reactive states
        df_idx = set.intersection(mda_idx, prev_idx, inter_idx, annot_idx)
        if df_idx:
            filtered_idx.set(df_idx)
            csv_filtered.set(csv.get().master.iloc[list(df_idx)])
        else:
            filtered_idx.set(None)
            csv_filtered.set(None)

    # =========================================================================
    # Rendering helper functions
    # =========================================================================
    def check_data():
        """Helper function to check if filtered data is available."""
        if csv_filtered.get() is None:
            return ui.p(ERROR_MESSAGES["no_file"])
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
            input.selected_interactions(),
            input.prevalence_threshold(),
        ]

        if is_main:
            plot_args.append(input.show_prevalence())

        fig = create_plot_func(*plot_args)

        if fig is None:
            return ui.p(ERROR_MESSAGES["no_interactions"])

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
    def interaction_plot():
        """Render the main interaction heatmap plot."""
        return render_plot(create_plot, is_main=True)

    @output
    @render.ui
    def ligand_interactions_plot():
        """Render the ligand interactions plot."""
        return render_plot(create_ligand_interactions_plot)

    @output
    @render.ui
    def receptor_interactions_plot():
        """Render the receptor interactions plot."""
        return render_plot(create_receptor_interactions_plot)

    @output
    @render.ui
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

app = App(app_ui, server)
