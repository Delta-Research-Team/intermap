"""
Main application file for InterMap Visualizations.
Integrates UI components with server-side logic.
"""

import os
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from tkinter import filedialog
import tkinter as tk
from shiny import App, reactive, render, ui

from .css import ERROR_MESSAGES
from .icsv import CSVFilter, sortby, transpose
from .plots import (create_interactions_over_time_plot,
                    create_ligand_interactions_plot, create_plot,
                    create_receptor_interactions_plot,
                    create_lifetime_plot, create_network_plot, create_3d_view,
                    InteractiveVisualization, filter_interactions_by_frame)
from .ui import app_ui


def server(input, output, session):
    """Define server logic for the app."""

    # =========================================================================
    # Reactive States
    # =========================================================================
    csv = reactive.Value(None)
    filtered_idx = reactive.Value(None)
    csv_filtered = reactive.Value(None)
    transpose_state = reactive.Value(False)
    frame_visualization = reactive.Value(None)

    @reactive.Effect
    @reactive.event(input.pickle_file, input.config_file)
    def initialize_filter():
        """Initialize CSVFilter when both files are uploaded."""
        try:
            pickle_infos = input.pickle_file()
            config_infos = input.config_file()

            if not pickle_infos or not config_infos:
                csv.set(None)
                filtered_idx.set(None)
                csv_filtered.set(None)
                return

            temp_dir = tempfile.mkdtemp()

            pickle_path = os.path.join(temp_dir, pickle_infos[0]["name"])
            config_path = os.path.join(temp_dir, config_infos[0]["name"])

            shutil.copy2(pickle_infos[0]["datapath"], pickle_path)
            shutil.copy2(config_infos[0]["datapath"], config_path)

            try:
                master_instance = CSVFilter(
                    pickle=pickle_path,
                    cfg=config_path
                )

                csv.set(master_instance)

                ui.notification_show(
                    "Files loaded successfully!",
                    duration=3000,
                    type="message"
                )

                ui.update_navs("plot_tabs", selected="Sele1 vs Sele2")

            except Exception as e:
                print(f"Error initializing CSVFilter: {str(e)}")
                raise e

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
    @reactive.event(input.sort_by)
    async def handle_sort():
        """Handle sorting when sort option changes."""
        if filtered_idx.get() is not None and csv.get() is not None:
            master_instance = csv.get()
            filtered_df = master_instance.master.iloc[
                list(filtered_idx.get())].copy()

            # Apply new sorting directly
            sorted_data = sortby(filtered_df, input.sort_by())

            # Update visualization
            await session.send_custom_message("refresh-plots", {})

    @reactive.Effect
    @reactive.event(input.transpose_button)
    def handle_transpose():
        """Handle the transpose button click."""
        transpose_state.set(not transpose_state.get())

        session.send_custom_message(
            "update-transpose-button",
            {"active": transpose_state.get()}
        )

        if filtered_idx.get() is not None and csv.get() is not None:
            master_instance = csv.get()
            filtered_df = csv.get().master.iloc[
                list(filtered_idx.get())].copy()

            if transpose_state.get():
                filtered_df = transpose(filtered_df)
                master_instance.axisx, master_instance.axisy = master_instance.axisy, master_instance.axisx

            csv_filtered.set(filtered_df)
            session.send_custom_message("refresh-plots", {})

    @reactive.Effect
    @reactive.event(input.plot_button,
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

        try:
            # Obtener idx's usando CSVFilter
            mda_idx, mda_status = master_instance.by_mda(
                input.mda_selection() if input.mda_selection() else 'all')
            prev_idx, prev_status = master_instance.by_prevalence(
                input.prevalence_threshold())
            inter_idx, inter_status = master_instance.by_inters(tuple(
                input.selected_interactions()) if input.selected_interactions() else input.selected_interactions())
            annot_idx, annot_status = master_instance.by_notes(tuple(
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
            elif annot_status == -1:
                ui.notification_show(
                    "There is not information for this Annotation selection"
                    , type="message", duration=5)

            # Update reactive states
            df_idx = set.intersection(mda_idx, prev_idx, inter_idx, annot_idx)

            if df_idx:
                filtered_df = master_instance.master.iloc[list(df_idx)].copy()
                if transpose_state.get():
                    filtered_df = transpose(filtered_df)

                filtered_idx.set(df_idx)
                csv_filtered.set(filtered_df)

                session.send_custom_message("refresh-plots", {})
            else:
                filtered_idx.set(None)
                csv_filtered.set(None)

        except Exception as e:
            print(f"Error in update_filtered_idx: {str(e)}")
            filtered_idx.set(None)
            csv_filtered.set(None)
            ui.notification_show(
                f"Error updating filters: {str(e)}",
                duration=5000,
                type="error"
            )

    @reactive.Effect
    @reactive.event(input.save_plot_trigger)
    def handle_download():
        """Handle saving the current plot."""
        if csv_filtered.get() is None:
            ui.notification_show(
                "No plot to save. Please generate a plot first.",
                duration=3000,
                type="warning"
            )
            return

        try:
            root = tk.Tk()
            root.withdraw()
            save_dir = filedialog.askdirectory(
                title="Select Directory to Save Plot"
            )

            if not save_dir:
                return

            save_dir = Path(save_dir)

            active_tab = input.plot_tabs()
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            master_instance = csv.get()

            if active_tab == "Life Time":
                fig = create_lifetime_plot(csv_filtered.get(),
                                           input.plot_width(),
                                           input.plot_height(),
                                           master_instance.axisx,
                                           master_instance.axisy,
                                           input.show_prevalence())
                filename = f"lifetime_plot_{timestamp}.html"
                full_path = save_dir / filename
                fig.write_html(str(full_path))

            elif active_tab == "Sele1 vs Sele2":
                fig = create_plot(csv_filtered.get(),
                                  input.plot_width(),
                                  input.plot_height(),
                                  master_instance.axisx,
                                  master_instance.axisy,
                                  input.show_prevalence())
                filename = f"sele1_vs_sele2_plot_{timestamp}.html"
                full_path = save_dir / filename
                fig.write_html(str(full_path))

            elif active_tab == "Prevalence":
                fig1 = create_ligand_interactions_plot(csv_filtered.get(),
                                                       input.plot_width(),
                                                       input.plot_height() // 1.5,
                                                       master_instance.axisx,
                                                       master_instance.axisy)
                filename1 = f"selection1_prevalence_{timestamp}.html"
                full_path1 = save_dir / filename1
                fig1.write_html(str(full_path1))

                fig2 = create_receptor_interactions_plot(csv_filtered.get(),
                                                         input.plot_width(),
                                                         input.plot_height() // 1.5,
                                                         master_instance.axisx,
                                                         master_instance.axisy)
                filename2 = f"selection2_prevalence_{timestamp}.html"
                full_path2 = save_dir / filename2
                fig2.write_html(str(full_path2))

                ui.notification_show(
                    f"Plots saved as:\n{filename1}\n{filename2}",
                    duration=5000,
                    type="message"
                )
                return

            elif active_tab == "Time Series":
                fig = create_interactions_over_time_plot(csv_filtered.get(),
                                                         input.plot_width(),
                                                         input.plot_height(),
                                                         master_instance.axisx,
                                                         master_instance.axisy)
                filename = f"time_series_plot_{timestamp}.html"
                full_path = save_dir / filename
                fig.write_html(str(full_path))

            elif active_tab == "Network":
                net = create_network_plot(
                    csv_filtered.get(),
                    input.plot_width(),
                    input.plot_height(),
                    master_instance.axisx,
                    master_instance.axisy
                )
                filename = f"network_plot_{timestamp}.html"
                full_path = save_dir / filename
                net.save_graph(str(full_path))

            if active_tab not in ["Prevalence", "Network"]:
                ui.notification_show(
                    f"Plot saved as:\n{filename}",
                    duration=5000,
                    type="message"
                )
            elif active_tab == "Network":
                ui.notification_show(
                    f"Network saved as:\n{filename}",
                    duration=5000,
                    type="message"
                )

        except Exception as e:
            print(f"Error saving plot: {str(e)}")
            ui.notification_show(
                f"Error saving plot: {str(e)}",
                duration=5000,
                type="error"
            )
        finally:
            if 'root' in locals():
                root.destroy()

    @reactive.Effect
    @reactive.event(input.update_frame_button)
    def update_3d_frame():
        """Update 3D visualization for specific frame."""
        if csv_filtered.get() is None or csv.get() is None:
            ui.notification_show(
                "No data loaded for 3D visualization",
                duration=3000,
                type="warning"
            )
            return

        try:
            master_instance = csv.get()
            frame_index = input.frame_number()

            # Create or update visualization
            if frame_visualization.get() is None:
                viz = InteractiveVisualization(csv_filtered.get(),
                                               master_instance.universe)
                frame_visualization.set(viz)
            else:
                viz = frame_visualization.get()
                viz.df = csv_filtered.get()

            # Update max frame in UI
            max_frames = viz.get_max_frames()
            ui.update_numeric("frame_number", max=max_frames)

            # Update frame info
            session.send_custom_message(
                "update-frame-info",
                {"current": frame_index, "max": max_frames}
            )

            # Update visualization
            viewer = viz.update_frame(frame_index)

            if viewer:
                ui.notification_show(
                    f"Frame {frame_index} loaded successfully",
                    duration=2000,
                    type="message"
                )
            else:
                ui.notification_show(
                    f"Error loading frame {frame_index}",
                    duration=3000,
                    type="error"
                )

        except Exception as e:
            print(f"Error updating 3D frame: {e}")
            ui.notification_show(
                f"Error updating frame: {str(e)}",
                duration=5000,
                type="error"
            )

    # =========================================================================
    # Rendering helper functions
    # =========================================================================
    def check_data():
        """Helper function to check if filtered data is available."""
        if csv_filtered.get() is None:
            return ui.p(ERROR_MESSAGES["nothing_selected"])
        return None

    def render_plot(create_plot_func, is_main=False, is_network=False):
        """Helper function to render plots with common logic."""
        error = check_data()
        if error is not None:
            return error

        master_instance = csv.get()
        if master_instance is None:
            return ui.p(ERROR_MESSAGES["nothing_selected"])

        plot_args = [
            csv_filtered.get(),
            input.plot_width(),
            input.plot_height() if (is_main or is_network) else input.plot_height() // 1.5,
            master_instance.axisx,
            master_instance.axisy
        ]

        if is_main:
            plot_args.append(input.show_prevalence())

        if is_network:
            net = create_plot_func(*plot_args)
            if net is None:
                return ui.p(ERROR_MESSAGES["nothing_selected"])

            unique_id = f"network_{hash(str(csv_filtered.get()))}"
            temp_file = os.path.join(tempfile.gettempdir(),
                                     f"{unique_id}.html")

            net.save_graph(temp_file)

            with open(temp_file, 'r', encoding='utf-8') as f:
                html_content = f.read()

            enhanced_content = html_content.replace('</head>',
                                                    """
                                                    <script>
                                                    document.addEventListener('DOMContentLoaded', function() {
                                                        setTimeout(function() {
                                                            if (typeof network !== 'undefined') {
                                                                network.fit();
                                                                network.redraw();
                                                            }
                                                        }, 100);
                                                    });
                                                    </script>
                                                    </head>
                                                    """)
            try:
                os.remove(temp_file)
            except:
                pass

            return ui.tags.div(
                ui.HTML(enhanced_content),
                id=unique_id,
                style=f"width: 100%; height: {input.plot_height()}px; border: none;"
            )
        else:
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
    @reactive.event(input.plot_button, input.transpose_button)
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

    @output
    @render.ui
    @reactive.event(input.plot_button)
    def lifetime_plot():
        """Render the lifetime boxplot."""
        error = check_data()
        if error is not None:
            return error

        master_instance = csv.get()
        if master_instance is None:
            return ui.p(ERROR_MESSAGES["nothing_selected"])

        fig = create_lifetime_plot(csv_filtered.get(),
                                   input.plot_width(),
                                   input.plot_height(),
                                   master_instance.axisx,
                                   master_instance.axisy,
                                   input.show_prevalence())

        if fig is None:
            return ui.p(ERROR_MESSAGES["nothing_selected"])

        plot_html = fig.to_html(
            include_plotlyjs="cdn",
            full_html=False,
            config={"responsive": True}
        )

        return ui.tags.div(
            ui.HTML(plot_html),
            style=f"width:100%; height:{input.plot_height()}px;"
        )

    @output
    @render.ui
    @reactive.event(input.plot_button)
    def network_plot():
        """Render the network visualization plot."""
        return render_plot(create_network_plot, is_network=True)

    @output
    @render.ui
    @reactive.event(input.plot_button, input.show_protein,
                    input.show_interactions, input.molecule_style,
                    input.interaction_sphere_size,
                    input.interaction_line_width,
                    input.update_frame_button)
    def molecule_3d_view():
        """Render the 3D molecular visualization."""
        try:
            master_instance = csv.get()
            if master_instance is None:
                return ui.p("No data loaded")

            # Get current frame
            frame_index = input.frame_number() if input.frame_number() is not None else 0

            # Filter interactions by frame
            df_filtered = filter_interactions_by_frame(csv_filtered.get(),
                                                       frame_index) if csv_filtered.get() is not None else None

            if df_filtered is None or df_filtered.empty:
                return ui.p(f"No interactions found for frame {frame_index}")

            view = create_3d_view(
                df_filtered,
                master_instance.universe,
                width=900,
                height=760,
                show_protein=input.show_protein(),
                show_interactions=input.show_interactions(),
                molecule_style=input.molecule_style(),
                sphere_size=input.interaction_sphere_size(),
                line_width=input.interaction_line_width(),
                frame_index=frame_index
            )

            if view is None:
                return ui.p("Failed to create visualization")

            html_content = f"""
            <div id="3dmol-viewer-{frame_index}" style="width: 100%; height: 600px; position: relative;">
                <script src="https://3dmol.org/build/3Dmol-min.js"></script>
                {view._make_html()}
            </div>
            """
            return ui.HTML(html_content)

        except Exception as e:
            print(f"Error in molecule_3d_view: {e}")
            return ui.p(f"Error: {str(e)}")

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
