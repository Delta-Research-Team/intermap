"""
Main application file for InterMap Visualizations.
Integrates UI components with server-side logic.
"""

import pandas as pd
import numpy as np
from shiny import App, render, reactive, ui, experimental
import os

from .components.ui import app_ui
from .components.plots import (
    create_plot,
    create_ligand_interactions_plot,
    create_receptor_interactions_plot,
    create_interactions_over_time_plot,
    initialize_search_state,
)
from .utils.helpers import (
    generate_interaction_choices,
    find_topology_file,
    validate_mda_selection,
)
from .components.sele_shiny import CSVFilter
from .config import all_interactions_colors, ERROR_MESSAGES

# =============================================================================
# Data Frame Processing
# =============================================================================


def filter_dataframe(df, search_term):
    """
    Filter DataFrame based on search term across specified columns.

    Args:
        df (pd.DataFrame): DataFrame to filter
        search_term (str): Term to search for

    Returns:
        pd.DataFrame: Filtered DataFrame
    """
    if not search_term:
        return df

    search_term = search_term.upper()

    mask = (
        df["sel1_atom"].str.contains(search_term,
                                     case=False, na=False)
        | df["sel2_atom"].str.contains(search_term,
                                       case=False, na=False)
        | df["interaction_name"].str.contains(search_term,
                                              case=False, na=False)
    )

    return df[mask]


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


# =============================================================================
# Shiny Application Server
# =============================================================================


def server(input, output, session):
    """Define server logic for the app."""

    # =========================================================================
    # Reactive States
    # =========================================================================

    data_store = reactive.Value({"df": None})
    search_state = reactive.Value({"search_term": "", "active": False})
    show_state = reactive.Value(False)
    topology_state = reactive.Value({"has_topology": False, "path": None})

    # Initialize search state
    initialize_search_state(search_state)

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
                "updateTopologyIndicator", {"hasTopology": False}
            )
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
        """Validate file and enable/disable show button."""

        file_infos = input.file()
        if file_infos:
            file_name = file_infos[0]["name"]
            if not file_name.endswith(".csv"):
                ui.notification_show(
                    ERROR_MESSAGES["invalid_file"],
                    type="error",
                    style="font-family: Ubuntu Mono; font-style: italic;",)
                ui.update_action_button("show_plots", disabled=True)
                return
            ui.update_action_button("show_plots", disabled=False)
        else:
            ui.update_action_button("show_plots", disabled=True)

        # Reset the plot state
        show_state.set(False)

    @reactive.Effect
    @reactive.event(input.show_plots)
    def handle_show():
        """Handle show button clicks and validate
        MDAnalysis selection if present.
        """
        # Check for topology file and MDAnalysis selection
        if topology_state.get()["has_topology"] and input.mda_selection():
            is_valid, error_msg = validate_mda_selection(input.mda_selection())
            if not is_valid:
                ui.notification_show(
                    ERROR_MESSAGES["invalid_sele"].format(error_msg),
                    type="error",
                    style="font-family: Ubuntu Mono; font-style: italic;",
                )

                return

        # Reset the display state without MDAnalysis selection.
        show_state.set(True)

    @reactive.Effect
    @reactive.event(input.search_button)
    def handle_search():
        """Update search state when search button is clicked."""
        search_state.set({"search_term": input.residue_filter(),
                          "active": True})

    # =========================================================================
    # UI Outputs
    # =========================================================================

    @output
    @render.ui
    def residue_not_found():
        """Show error message when residue is not found."""

        # Check if there is an active search and if data is available.
        if not search_state.get()["active"] or data_store.get()["df"] is None:
            return ui.div()

        # Obtain DataFrame and search term.
        df = data_store.get()["df"]
        search_term = search_state.get()["search_term"].upper()
        if not search_term:
            return ui.div()

        # Filter DataFrame based on the search term.
        filtered_df = filter_dataframe(df, search_term)

        # If there are no results, display an error message.
        if filtered_df.empty:
            return ui.p("Not Found", style="color: red; margin-top: 5px;")
        return ui.div()

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
        input.prevalence_threshold,
        input.show_prevalence,
        input.search_button,
        input.show_plots,
        input.mda_selection,
    )
    def interaction_plot():
        """
        Render the main interaction heatmap plot.
        """

        if not show_state():
            return None

        file_infos = input.file()
        if not file_infos:
            return ui.p(
                ERROR_MESSAGES["no_file"],
                style="font-family: Ubuntu Mono; font-style: italic;",
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
                            type="error",)
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
                        df = pd.read_csv(csv_path, sep=None, engine="python")
                    except:
                        try:
                            df = pd.read_csv(csv_path, encoding="utf-8")
                        except UnicodeDecodeError:
                            df = pd.read_csv(csv_path, encoding="latin-1")
                        except Exception as e:
                            print(f"Error reading file: {str(e)}")
                            raise

                # Clean names and verify required column
                df.columns = df.columns.str.strip()
                required_columns = [
                    "sel1_atom",
                    "sel2_atom",
                    "interaction_name",
                    "prevalence",
                    "timeseries",
                ]
                missing_columns = [
                    col for col in required_columns if col not in df.columns
                ]
                if missing_columns:
                    raise ValueError(
                        f"Missing required columns: " 
                        f"{', '.join(missing_columns)}"
                    )

                # Process timeseries
                df = process_timeseries_numpy(df)
                data_store.set({"df": df})
            else:
                df = data_store.get()["df"]

            # Apply search filter if active
            if search_state.get()["active"]:
                search_term = search_state.get()["search_term"]
                df = filter_dataframe(df, search_term)
                if df.empty:
                    return ui.p(
                        "No matches found for the search term.",
                        style="color: red; font-family: Ubuntu Mono;",
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
                "margin: 10px 0;",)

        except Exception as e:
            import traceback

            print(f"Error in plot generation:")
            print(traceback.format_exc())
            error_message = ERROR_MESSAGES["plot_error"].format(str(e))
            return ui.div(
                ui.p(error_message, style="color:red; "
                                          "font-family: Ubuntu Mono;"),
                ui.pre(
                    traceback.format_exc(),
                    style="font-family: Ubuntu Mono; font-size: 12px; "
                    "background-color: #f8f9fa; padding: 10px; "
                    "border-radius: 4px;",),)

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
        input.prevalence_threshold,
        input.search_button,
        input.show_plots,
    )
    def ligand_interactions_plot():
        """Render the ligand interactions plot.
        """
        if not show_state():
            return None

        # Obtain stored data
        data = data_store.get()
        if data["df"] is None:
            return ui.p(ERROR_MESSAGES["no_file"])

        try:
            df = data["df"]
            # Apply search filter if active
            if search_state.get()["active"]:
                search_term = search_state.get()["search_term"]
                df = filter_dataframe(df, search_term)
                if df.empty:
                    return ui.p("No matches found for the search term.")

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
                style=f"width:100%; height:{input.plot_height() // 2}px;",)

        except Exception as e:
            import traceback
            print(f"Error in plot generation:")
            print(traceback.format_exc())
            error_message = ERROR_MESSAGES["plot_error"].format(str(e))
            return ui.div(
                ui.p(error_message, style="color:red; "
                                          "font-family: Ubuntu Mono;"),
                ui.pre(
                    traceback.format_exc(),
                    style="font-family: Ubuntu Mono; font-size: 12px;",),)

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
        input.prevalence_threshold,
        input.search_button,
        input.show_plots,
    )
    def receptor_interactions_plot():
        """Render the receptor interactions plot.
        """
        if not show_state():
            return None

        # Obtain stored data
        data = data_store.get()
        if data["df"] is None:
            return ui.p(ERROR_MESSAGES["no_file"])

        try:
            df = data["df"]
            # Apply search filter if active
            if search_state.get()["active"]:
                search_term = search_state.get()["search_term"]
                df = filter_dataframe(df, search_term)
                if df.empty:
                    return ui.p("No matches found for the search term.")

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
                style=f"width:100%; height:{input.plot_height() // 2}px;",)

        except Exception as e:
            import traceback
            print(f"Error in plot generation:")
            print(traceback.format_exc())
            error_message = ERROR_MESSAGES["plot_error"].format(str(e))
            return ui.div(
                ui.p(error_message, style="color:red; "
                                          "font-family: Ubuntu Mono;"),
                ui.pre(
                    traceback.format_exc(),
                    style="font-family: Ubuntu Mono; font-size: 12px;",
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
        input.prevalence_threshold,
        input.search_button,
        input.show_plots,
    )
    def interactions_over_time_plot():
        """Render the interactions over time plot."""
        if not show_state():
            return None

        # Obtain stored data
        data = data_store.get()
        if data["df"] is None:
            return ui.p(ERROR_MESSAGES["no_file"])

        try:
            df = data["df"]
            # Apply search filter if active
            if search_state.get()["active"]:
                search_term = search_state.get()["search_term"]
                df = filter_dataframe(df, search_term)
                if df.empty:
                    return ui.p("No matches found for the search term.")

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
                style=f"width:100%; height:{input.plot_height() // 2}px;",)

        except Exception as e:
            import traceback
            print(f"Error in plot generation:")
            print(traceback.format_exc())
            error_message = ERROR_MESSAGES["plot_error"].format(str(e))
            return ui.div(
                ui.p(error_message, style="color:red; "
                                          "font-family: Ubuntu Mono;"),
                ui.pre(
                    traceback.format_exc(),
                    style="font-family: Ubuntu Mono; font-size: 12px;",),)

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
            selected=list(all_interactions_colors.keys()),)


# Create the Shiny app instance
app = App(app_ui, server)
