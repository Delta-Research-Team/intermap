"""
Main application file for InterMap Visualizations.
Integrates UI components with server-side logic.
"""

import os
import shutil
import tempfile
import tkinter as tk
from datetime import datetime
from pathlib import Path
from tkinter import filedialog

from shiny import App, reactive, render, ui

from .css import ERROR_MESSAGES
from .icsv import CSVFilter, transpose
from .plots import (create_interactions_over_time_plot,
                    create_lifetime_plot, create_network_plot, create_plot,
                    create_sel1_interactions_plot,
                    create_sel2_interactions_plot)
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
    simplify_labels = reactive.Value(False)
    custom_x_axis = reactive.Value(None)
    custom_y_axis = reactive.Value(None)

    heatmap_title = reactive.Value("Interaction Heatmap")
    prevalence_title = reactive.Value("Interaction Prevalence Analysis")
    lifetime_title = reactive.Value("Interaction Lifetimes Distribution")
    timeseries_title = reactive.Value("Interactions Over Time")
    network_title = reactive.Value("Interaction Network")

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
                fig1 = create_sel1_interactions_plot(csv_filtered.get(),
                                                     input.plot_width(),
                                                     input.plot_height() // 1.5,
                                                     master_instance.axisx,
                                                     master_instance.axisy)
                filename1 = f"selection1_prevalence_{timestamp}.html"
                full_path1 = save_dir / filename1
                fig1.write_html(str(full_path1))

                fig2 = create_sel2_interactions_plot(csv_filtered.get(),
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


    # =========================================================================
    # Rendering helper functions
    # =========================================================================
    def check_data():
        """Helper function to check if filtered data is available."""
        if csv_filtered.get() is None:
            return ui.p(ERROR_MESSAGES["nothing_selected"])
        return None

    def process_axis_labels(df, simplify=False):
        """
        Process axis labels based on user preferences.

        Args:
            df (pd.DataFrame): DataFrame to process
            simplify (bool): Whether to simplify axis labels

        Returns:
            pd.DataFrame: Processed DataFrame
        """
        if df is None or not simplify:
            return df

        processed_df = df.copy()

        # Simplify selection columns (main axes)
        for col in ['sel1', 'sel2']:
            if col in processed_df.columns:
                processed_df[col] = processed_df[col].apply(
                    lambda x: '_'.join(x.split('_')[:-1]) if isinstance(x,
                                                                        str) and '_' in x else x
                )

        # Handle pair column for lifetime plot
        if 'pair' in processed_df.columns:
            processed_df['pair'] = processed_df['pair'].apply(
                lambda x: simplify_pair_label(x) if isinstance(x, str) else x
            )

        # Handle selection_pair column for time series plot
        if 'selection_pair' in processed_df.columns:
            processed_df['selection_pair'] = processed_df[
                'selection_pair'].apply(
                lambda x: simplify_pair_label(x) if isinstance(x, str) else x
            )

        return processed_df

    def simplify_pair_label(pair_text):
        """
        Simplify a pair label like "sel1 - sel2 (interaction)"
        by removing the last part of each selection.
        """
        if not isinstance(pair_text, str) or ' - ' not in pair_text:
            return pair_text

        # Split into selections and interaction part
        parts = pair_text.split(' (')
        if len(parts) != 2:
            return pair_text

        selections = parts[0]
        interaction_part = f" ({parts[1]}" if len(parts) > 1 else ""

        # Split the selections
        sel_parts = selections.split(' - ')
        if len(sel_parts) != 2:
            return pair_text

        # Simplify each selection
        simplified_sels = []
        for sel in sel_parts:
            if '_' in sel:
                simplified_sels.append('_'.join(sel.split('_')[:-1]))
            else:
                simplified_sels.append(sel)

        # Combine back
        return f"{simplified_sels[0]} - {simplified_sels[1]}{interaction_part}"

    def render_plot(create_plot_func, is_main=False, is_network=False):
        """Helper function to render plots with common logic."""
        error = check_data()
        if error is not None:
            return error

        master_instance = csv.get()
        if master_instance is None:
            return ui.p(ERROR_MESSAGES["nothing_selected"])

        # Get the current active tab
        active_tab = input.plot_tabs()

        # Get the DataFrame
        df_to_use = csv_filtered.get()

        # Apply label simplification to all plots if enabled
        if simplify_labels.get():
            df_to_use = process_axis_labels(df_to_use, True)

        # Get axis titles
        x_axis_title = master_instance.axisx
        y_axis_title = master_instance.axisy

        # Apply custom axis titles if appropriate for this tab
        if active_tab in ["Sele1 vs Sele2", "Prevalence"]:
            if custom_x_axis.get() is not None:
                x_axis_title = custom_x_axis.get()
            if custom_y_axis.get() is not None:
                y_axis_title = custom_y_axis.get()

        # Determine which title to use based on the active tab
        plot_title = None
        if active_tab == "Sele1 vs Sele2":
            plot_title = heatmap_title.get()
        elif active_tab == "Prevalence":
            plot_title = prevalence_title.get()
        elif active_tab == "Life Time":
            plot_title = lifetime_title.get()
        elif active_tab == "Time Series":
            plot_title = timeseries_title.get()
        elif active_tab == "Network":
            plot_title = network_title.get()

        plot_args = [
            df_to_use,
            input.plot_width(),
            input.plot_height() if (
                        is_main or is_network) else input.plot_height() // 1.5,
            x_axis_title,
            y_axis_title
        ]

        if is_main:
            plot_args.append(input.show_prevalence())

        # Create the figure
        if is_network:
            net = create_plot_func(*plot_args)
            if net is None:
                return ui.p(ERROR_MESSAGES["nothing_selected"])

            # Apply network title
            unique_id = f"network_{hash(str(csv_filtered.get()))}"
            temp_file = os.path.join(tempfile.gettempdir(),
                                     f"{unique_id}.html")
            net.save_graph(temp_file)

            with open(temp_file, 'r', encoding='utf-8') as f:
                html_content = f.read()

            # Insert title into the HTML
            title_html = f'<h3 style="text-align:center;font-family:Roboto;margin-bottom:20px;">{plot_title}</h3>'
            enhanced_content = html_content.replace('<body>',
                                                    f'<body>{title_html}')

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

            # Apply plot title if available
            if plot_title:
                fig.update_layout(title={
                    'text': plot_title,
                    'x': 0.5,
                    'font': {'family': "Roboto", 'size': 20}
                })

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
    @reactive.event(input.plot_button, input.transpose_button,
                    input.simplify_axis_labels, input.apply_axis_names,
                    input.apply_plot_titles)
    def interaction_plot():
        """Render the main interaction heatmap plot."""
        return render_plot(create_plot, is_main=True)

    @output
    @render.ui
    @reactive.event(input.plot_button, input.simplify_axis_labels,
                    input.apply_axis_names, input.apply_plot_titles)
    def sel1_interactions_plot():
        """Render the sel1 interactions plot."""
        return render_plot(create_sel1_interactions_plot)

    @output
    @render.ui
    @reactive.event(input.plot_button, input.simplify_axis_labels,
                    input.apply_axis_names, input.apply_plot_titles)
    def sel2_interactions_plot():
        """Render the sel2 interactions plot."""
        return render_plot(create_sel2_interactions_plot)

    @output
    @render.ui
    @reactive.event(input.plot_button, input.simplify_axis_labels,
                    input.apply_plot_titles)
    def interactions_over_time_plot():
        """Render the interactions over time plot."""
        return render_plot(create_interactions_over_time_plot)

    @output
    @render.ui
    @reactive.event(input.plot_button, input.simplify_axis_labels,
                    input.apply_plot_titles)
    def lifetime_plot():
        """Render the lifetime boxplot."""
        error = check_data()
        if error is not None:
            return error

        master_instance = csv.get()
        if master_instance is None:
            return ui.p(ERROR_MESSAGES["nothing_selected"])

        # Apply label simplification if enabled
        df_to_use = csv_filtered.get()
        if simplify_labels.get():
            df_to_use = process_axis_labels(df_to_use, True)

        fig = create_lifetime_plot(df_to_use,
                                   input.plot_width(),
                                   input.plot_height(),
                                   master_instance.axisx,
                                   master_instance.axisy,
                                   input.show_prevalence())

        if fig is None:
            return ui.p(ERROR_MESSAGES["nothing_selected"])

        # Apply title
        fig.update_layout(title={
            'text': lifetime_title.get(),
            'x': 0.5,
            'font': {'family': "Roboto", 'size': 20}
        })

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
    @reactive.event(input.plot_button, input.simplify_axis_labels,
                    input.apply_plot_titles)
    def network_plot():
        """Render the network visualization plot."""
        return render_plot(create_network_plot, is_network=True)

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

    @reactive.Effect
    @reactive.event(input.simplify_axis_labels)
    def update_simplify_labels():
        """Update axis label simplification setting."""
        simplify_labels.set(input.simplify_axis_labels())
        if filtered_idx.get() is not None:
            session.send_custom_message("refresh-plots", {})

    @reactive.Effect
    @reactive.event(input.apply_axis_names)
    def update_custom_axes():
        """Update custom axis names."""
        custom_x_axis.set(
            input.custom_x_axis() if input.custom_x_axis() else None)
        custom_y_axis.set(
            input.custom_y_axis() if input.custom_y_axis() else None)
        if filtered_idx.get() is not None:
            session.send_custom_message("refresh-plots", {})

    @reactive.Effect
    @reactive.event(input.apply_plot_titles)
    def update_plot_titles():
        """Update custom plot titles."""
        heatmap_title.set(
            input.heatmap_plot_title() if input.heatmap_plot_title() else "Interaction Heatmap")
        prevalence_title.set(
            input.prevalence_plot_title() if input.prevalence_plot_title() else "Interaction Prevalence Analysis")
        lifetime_title.set(
            input.lifetime_plot_title() if input.lifetime_plot_title() else "Interaction Lifetimes Distribution")
        timeseries_title.set(
            input.timeseries_plot_title() if input.timeseries_plot_title() else "Interactions Over Time")
        network_title.set(
            input.network_plot_title() if input.network_plot_title() else "Interaction Network")

        if filtered_idx.get() is not None:
            session.send_custom_message("refresh-plots", {})


app = App(app_ui, server)
