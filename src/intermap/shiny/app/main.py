"""
Main application file for InterMap Visualizations.
Integrates UI components with server-side logic.
"""
import os

import numpy as np
import pandas as pd
from shiny import App, reactive, render, ui

from intermap.shiny.app.plots import (create_interactions_over_time_plot,
                                      create_ligand_interactions_plot, create_plot,
                                      create_receptor_interactions_plot)
from intermap.shiny.app.sele_shiny import CSVFilter
from intermap.shiny.app.ui import app_ui
from .css import all_interactions_colors, ERROR_MESSAGES
from intermap.shiny.app.helpers import (find_topology_file, generate_interaction_choices,
                                        validate_mda_selection)


# =============================================================================
# Data Frame Processing
# =============================================================================


def process_timeseries_numpy(df):
    """
    Process timeseries data from the df into frame-by-frame columns.

    Args:
        df (pd.DataFrame): df containing a 'timeseries' column with binary
            strings representing the presence (1) or absence (0) of
            interactions per frame.

    Returns:
        pd.DataFrame: Original df with additional numbered columns (0, 1, ...)
            representing each frame's binary value
    """
    if df.empty:
        return df

    timeseries = df["timeseries"].str.strip()
    first_series = timeseries.iloc[0]
    num_frames = len(first_series)
    matrix = np.array([list(x) for x in timeseries.values], dtype=int)
    time_df = pd.DataFrame(
        matrix, columns=[str(i) for i in range(num_frames)], index=df.index
    )
    return pd.concat([df, time_df], axis=1)


def generate_annotation_choices(df):
    """
    Generate choices for annotation checkboxes from note1 and note2 columns.
    Excludes empty or NaN values.

    Args:
        df (pd.DataFrame): DataFrame containing note1 and note2 columns

    Returns:
        list: Sorted list of unique annotations from both columns, excluding empty values
    """
    if df is None:
        return []

    # Filtrar valores vacíos y NaN de ambas columnas antes de obtener valores únicos
    note1_values = set(df['note1'][df['note1'].notna() & (df['note1'] != '')].unique())
    note2_values = set(df['note2'][df['note2'].notna() & (df['note2'] != '')].unique())

    # Unir los conjuntos y ordenar
    all_notes = sorted(note1_values.union(note2_values))

    return all_notes


def filter_by_annotations(df, selected_annotations):
    """
    Filter DataFrame based on selected annotations.

    Args:
        df (pd.DataFrame): DataFrame to filter
        selected_annotations (list): List of selected annotation values

    Returns:
        pd.DataFrame: Filtered DataFrame
    """
    if not selected_annotations:
        return df

    mask = (
        df['note1'].isin(selected_annotations) |
        df['note2'].isin(selected_annotations)
    )
    return df[mask]


# =============================================================================
# Shiny Application Server
# =============================================================================


def server(input, output, session):
    """Define server logic for the app."""

    # =========================================================================
    # Reactive States
    # =========================================================================

    data_store = reactive.Value({"df": None})
    topology_state = reactive.Value({"has_topology": False, "path": None})

    @reactive.Effect
    def _():
        """Global error handler"""
        try:
            yield
        except Exception as e:
            import traceback
            print("Error in application:")
            print(traceback.format_exc())
            ui.notification_show(f"Error: {str(e)}",
                                 type="error",
                                 duration=None)

    @reactive.Effect
    @reactive.event(input.file)
    async def check_topology():
        """Check for topology file when CSV is uploaded."""

        file_infos = input.file()
        if not file_infos:
            topology_state.set({"has_topology": False, "path": None})
            await session.send_custom_message(
                "updateTopologyIndicator", {"hasTopology": False})
            return

        # Obtain file paths
        temp_csv_path = file_infos[0]["datapath"]
        original_filename = file_infos[0]["name"]
        current_dir = os.getcwd()
        original_path = os.path.join(current_dir, original_filename)

        # Obtain file paths
        has_topo, topo_path = find_topology_file(temp_csv_path, original_path)
        topology_state.set({"has_topology": has_topo, "path": topo_path})

        # Update visual indicator
        await session.send_custom_message(
            "updateTopologyIndicator", {"hasTopology": has_topo})
        if has_topo:
            ui.notification_show("Topology file found!", type="message",
                                 duration=5)
        else:
            ui.notification_show(
                ERROR_MESSAGES["no_topology"], type="warning", duration=5)

    @reactive.Effect
    @reactive.event(input.file)
    def validate_file():
        """Validate file format."""
        file_infos = input.file()
        if file_infos:
            file_name = file_infos[0]["name"]
            if not file_name.endswith(".csv"):
                ui.notification_show(
                    ERROR_MESSAGES["invalid_file"],
                    type="error",
                    style="font-family: Roboto; font-style: italic;",
                )
                return

    # =========================================================================
    # Render the Interaction Plot (first plot)
    # =========================================================================
    @output
    @render.ui
    @reactive.event(
        input.file,
        input.plot_width,
        input.plot_height,
        input.selected_interactions,
        input.selected_annotations,
        input.prevalence_threshold,
        input.show_prevalence,
        input.mda_selection,
    )
    def interaction_plot():
        """
        Render the main interaction heatmap plot.
        """
        file_infos = input.file()
        if not file_infos:
            return ui.p(
                ERROR_MESSAGES["no_file"],
                style="font-family: Roboto; font-style: italic;",
            )

        try:
            if data_store.get()["df"] is None:
                csv_path = file_infos[0]["datapath"]
                selection = input.mda_selection()
                topo_state = topology_state.get()

                if topo_state["has_topology"] and selection:
                    # Validate MDAnalysis selection
                    is_valid, error_msg = validate_mda_selection(selection)
                    if not is_valid:
                        ui.notification_show(
                            ERROR_MESSAGES["invalid_sele"].format(error_msg),
                            type="error", )
                        return ui.p(
                            ERROR_MESSAGES["invalid_sele"].format(error_msg))

                    try:
                        # Use CSVFilter if topology and selection are present
                        filter_obj = CSVFilter(topo_state["path"],
                                               selection, csv_path)
                        df = filter_obj.df
                    except Exception as e:
                        ui.notification_show(str(e), type="error")
                        return ui.p(f"Error: {str(e)}")
                else:
                    # Normal CSV reading
                    try:
                        df = pd.read_csv(csv_path, skiprows=1, sep=None, engine="python")
                    except:
                        try:
                            df = pd.read_csv(csv_path, skiprows=1, encoding="utf-8")
                        except UnicodeDecodeError:
                            df = pd.read_csv(csv_path, skiprows=1, encoding="latin-1")
                        except Exception as e:
                            print(f"Error reading file: {str(e)}")
                            raise

                # Process timeseries
                df = process_timeseries_numpy(df)
                data_store.set({"df": df})
            else:
                df = data_store.get()["df"]

            # Apply annotations filter
            df = filter_by_annotations(df, input.selected_annotations())
            if df.empty:
                return ui.p(
                    "No data available for selected annotations.",
                    style="color: red; font-family: Roboto;",
                )

            fig = create_plot(
                df,
                input.plot_width(),
                input.plot_height(),
                input.selected_interactions(),
                input.prevalence_threshold(),
                input.show_prevalence(),
            )

            if fig is None:
                return ui.p(ERROR_MESSAGES["no_interactions"])

            # Convert the graph to HTML.
            plot_html = fig.to_html(
                include_plotlyjs="cdn", full_html=False,
                config={"responsive": True})

            # Return the container of the graph.
            return ui.tags.div(
                ui.HTML(plot_html),
                style="width:100%; height:800px; border:1px solid #ddd; "
                      "margin: 10px 0;", )

        except Exception as e:
            import traceback

            print(f"Error in plot generation:")
            print(traceback.format_exc())
            error_message = ERROR_MESSAGES["plot_error"].format(str(e))
            return ui.div(
                ui.p(error_message, style="color:red; "
                                          "font-family: Roboto;"),
                ui.pre(
                    traceback.format_exc(),
                    style="font-family: Roboto; font-size: 12px; "
                          "background-color: #f8f9fa; padding: 10px; "
                          "border-radius: 4px;", ), )

    # =========================================================================
    # Render the Ligand Interaction Plot (second plot)
    # =========================================================================
    @output
    @render.ui
    @reactive.event(
        input.file,
        input.plot_width,
        input.plot_height,
        input.selected_interactions,
        input.selected_annotations,
        input.prevalence_threshold,
    )
    def ligand_interactions_plot():
        """Render the ligand interactions plot.
        """
        # Obtain stored data
        data = data_store.get()
        if data["df"] is None:
            return ui.p(ERROR_MESSAGES["no_file"])

        try:
            df = data["df"]

            # Apply annotations filter
            df = filter_by_annotations(df, input.selected_annotations())
            if df.empty:
                return ui.p("No data available for selected annotations.")

            # Create ligand interaction graph
            fig = create_ligand_interactions_plot(
                df,
                input.plot_width(),
                input.plot_height() // 2,
                input.selected_interactions(),
                input.prevalence_threshold(),
            )

            if fig is None:
                return ui.p("Loading ...")

            # Convert the graph to HTML.
            plot_html = fig.to_html(
                include_plotlyjs="cdn", full_html=False,
                config={"responsive": True})
            return ui.tags.div(
                ui.HTML(plot_html),
                style=f"width:100%; height:{input.plot_height() // 2}px;", )

        except Exception as e:
            import traceback
            print(f"Error in plot generation:")
            print(traceback.format_exc())
            error_message = ERROR_MESSAGES["plot_error"].format(str(e))
            return ui.div(
                ui.p(error_message, style="color:red; "
                                          "font-family: Roboto;"),
                ui.pre(
                    traceback.format_exc(),
                    style="font-family: Roboto; font-size: 12px;", ), )

    # =========================================================================
    # Render the Receptor Interaction Plot (third plot)
    # =========================================================================
    @output
    @render.ui
    @reactive.event(
        input.file,
        input.plot_width,
        input.plot_height,
        input.selected_interactions,
        input.selected_annotations,
        input.prevalence_threshold,
    )
    def receptor_interactions_plot():
        """Render the receptor interactions plot.
        """
        # Obtain stored data
        data = data_store.get()
        if data["df"] is None:
            return ui.p(ERROR_MESSAGES["no_file"])

        try:
            df = data["df"]

            # Apply annotations filter
            df = filter_by_annotations(df, input.selected_annotations())
            if df.empty:
                return ui.p("No data available for selected annotations.")

            # Create receptor interaction graph
            fig = create_receptor_interactions_plot(
                df,
                input.plot_width(),
                input.plot_height() // 2,
                input.selected_interactions(),
                input.prevalence_threshold(),
            )

            if fig is None:
                return ui.p("Loading ...")

            # Convert the graph to HTML.
            plot_html = fig.to_html(
                include_plotlyjs="cdn", full_html=False,
                config={"responsive": True})

            # Return the container of the graph.
            return ui.tags.div(
                ui.HTML(plot_html),
                style=f"width:100%; height:{input.plot_height() // 2}px;", )

        except Exception as e:
            import traceback
            print(f"Error in plot generation:")
            print(traceback.format_exc())
            error_message = ERROR_MESSAGES["plot_error"].format(str(e))
            return ui.div(
                ui.p(error_message, style="color:red; "
                                          "font-family: Roboto;"),
                ui.pre(
                    traceback.format_exc(),
                    style="font-family: Roboto; font-size: 12px;",
                ),
            )

    # =========================================================================
    # Render the Interaction Over Time Plot (fourth plot)
    # =========================================================================
    @output
    @render.ui
    @reactive.event(
        input.file,
        input.plot_width,
        input.plot_height,
        input.selected_interactions,
        input.selected_annotations,
        input.prevalence_threshold,
    )
    def interactions_over_time_plot():
        """Render the interactions over time plot."""
        # Obtain stored data
        data = data_store.get()
        if data["df"] is None:
            return ui.p(ERROR_MESSAGES["no_file"])

        try:
            df = data["df"]

            # Apply annotations filter
            df = filter_by_annotations(df, input.selected_annotations())
            if df.empty:
                return ui.p("No data available for selected annotations.")

            fig = create_interactions_over_time_plot(
                df,
                input.plot_width(),
                input.plot_height() // 2,
                input.selected_interactions(),
                input.prevalence_threshold(),
            )

            if fig is None:
                return ui.p("Loading ...")
            # Convert the graph to HTML.
            plot_html = fig.to_html(
                include_plotlyjs="cdn", full_html=False,
                config={"responsive": True})
            return ui.tags.div(
                ui.HTML(plot_html),
                style=f"width:100%; height:{input.plot_height() // 2}px;", )

        except Exception as e:
            import traceback
            print(f"Error in plot generation:")
            print(traceback.format_exc())
            error_message = ERROR_MESSAGES["plot_error"].format(str(e))
            return ui.div(
                ui.p(error_message, style="color:red; "
                                          "font-family: Roboto;"),
                ui.pre(
                    traceback.format_exc(),
                    style="font-family: Roboto; font-size: 12px;", ), )

    # =========================================================================
    # Render the Interaction Checkboxes
    # =========================================================================
    @output
    @render.ui
    def interaction_checkboxes():
        """Render the interaction checkboxes."""
        df = data_store.get()["df"]
        choices = generate_interaction_choices(df)
        return ui.input_checkbox_group(
            "selected_interactions",
            "",
            choices=choices,
            selected=list(all_interactions_colors.keys()), )

    # =========================================================================
    # Render the Annotation Checkboxes
    # =========================================================================
    @output
    @render.ui
    def annotation_checkboxes():
        """Render the annotation checkboxes."""
        df = data_store.get()["df"]
        choices = generate_annotation_choices(df)
        return ui.input_checkbox_group(
            "selected_annotations",
            "",
            choices=choices,
            selected=choices,  # Por defecto, todas las anotaciones seleccionadas
        )


# Create the Shiny app instance
app = App(app_ui, server)
