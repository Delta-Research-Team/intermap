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


# =========================================================================
# Pure Data Processing Functions (Extracted for Testing)
# =========================================================================

def get_axis_names(master_instance):
    """Get axis names with fallback to 'sel1' and 'sel2'."""
    if master_instance is None:
        return "sel1", "sel2"
    x_axis = getattr(master_instance, 'axisx', None) or "sel1"
    y_axis = getattr(master_instance, 'axisy', None) or "sel2"
    return x_axis.strip() or "sel1", y_axis.strip() or "sel2"


def simplify_pair_label(pair_text):
    """Simplify a pair label like 'sel1 - sel2 (interaction)'."""
    if not isinstance(pair_text, str) or ' - ' not in pair_text:
        return pair_text
    parts = pair_text.split(' (')
    if len(parts) != 2:
        return pair_text
    selections = parts[0]
    interaction_part = f" ({parts[1]}" if len(parts) > 1 else ""
    sel_parts = selections.split(' - ')
    if len(sel_parts) != 2:
        return pair_text
    simplified_sels = []
    for sel in sel_parts:
        if '_' in sel:
            simplified_sels.append('_'.join(sel.split('_')[:-1]))
        else:
            simplified_sels.append(sel)
    return f"{simplified_sels[0]} - {simplified_sels[1]}{interaction_part}"


def process_axis_labels(df, simplify=False):
    """Process axis labels based on user preferences."""
    if df is None or not simplify:
        return df
    processed_df = df.copy()
    for col in ['sel1', 'sel2']:
        if col in processed_df.columns:
            processed_df[col] = processed_df[col].apply(
                lambda x: '_'.join(x.split('_')[:-1]) if isinstance(x,
                                                                    str) and '_' in x else x
            )
    if 'pair' in processed_df.columns:
        processed_df['pair'] = processed_df['pair'].apply(
            lambda x: simplify_pair_label(x) if isinstance(x, str) else x
        )
    if 'selection_pair' in processed_df.columns:
        processed_df['selection_pair'] = processed_df['selection_pair'].apply(
            lambda x: simplify_pair_label(x) if isinstance(x, str) else x
        )
    return processed_df


def apply_organization_to_figure(fig, method, df=None):
    """Apply organization to a plotly figure's x-axis."""
    if fig is None or getattr(fig, "data", None) is None or len(
            fig.data) == 0 or not hasattr(fig.data[0], 'x') or fig.data[
        0].x is None:
        return fig
    try:
        x_values = []
        for trace in fig.data:
            if hasattr(trace, 'x') and trace.x is not None:
                x_values.extend(trace.x)
        x_values = list(dict.fromkeys(x_values))
        if not x_values:
            return fig

        x_components = {}
        for x in x_values:
            if isinstance(x, str):
                try:
                    if '__' in x:
                        first_part, rest = x.split('__', 1)
                        if '_' in first_part:
                            parts = first_part.split('_')
                            resname = parts[0]
                            resid = int(parts[1]) if len(parts) > 1 and parts[
                                1].isdigit() else 0
                        else:
                            resname = first_part
                            resid = 0
                    else:
                        parts = x.split('_')
                        resname = parts[0] if parts else ""
                        resid = int(parts[1]) if len(parts) > 1 and parts[
                            1].isdigit() else 0
                    x_components[x] = {'resname': resname, 'resid': resid,
                                       'full': x}
                except Exception:
                    x_components[x] = {'resname': x, 'resid': 0, 'full': x}

        if not x_components:
            return fig

        if method == "resname":
            sorted_x = sorted(x_components.keys(),
                              key=lambda x: x_components[x]['resname'])
        elif method == "resnum":
            sorted_x = sorted(x_components.keys(),
                              key=lambda x: x_components[x]['resid'])
        elif method == "annotation":
            if df is not None and 'note1' in df.columns and 'sel1' in df.columns:
                sel1_to_note1 = dict(zip(df['sel1'], df['note1']))
                sorted_x = sorted(x_components.keys(),
                                  key=lambda x: str(sel1_to_note1.get(x, "")))
            else:
                sorted_x = x_values
        else:
            sorted_x = sorted(x_components.keys(),
                              key=lambda x: x_components[x]['resname'])

        fig.update_layout(
            xaxis=dict(categoryorder='array', categoryarray=sorted_x))
        return fig
    except Exception:
        return fig


# =========================================================================
# Shiny Server Logic
# =========================================================================

def server(input, output, session):
    csv = reactive.Value(None)
    filtered_idx = reactive.Value(None)
    csv_filtered = reactive.Value(None)
    transpose_state = reactive.Value(False)
    simplify_labels = reactive.Value(False)
    custom_x_axis = reactive.Value(None)
    custom_y_axis = reactive.Value(None)
    organization_method = reactive.Value("resname")

    network_params = reactive.Value({
        'gravity': -200, 'central_gravity': 0.005, 'spring_length': 200,
        'spring_constant': 0.5, 'avoid_overlap': 0.8,
        'stabilization_iterations': 1000,
        'physics_enabled': True, 'min_node_size': 20, 'max_node_size': 50,
        'min_edge_width': 5, 'max_edge_width': 15
    })

    heatmap_title = reactive.Value("Interaction Heatmap")
    prevalence_title = reactive.Value("Interaction Prevalence Analysis")
    lifetime_title = reactive.Value("Interaction Lifetimes Distribution")
    timeseries_title = reactive.Value("Interactions Over Time")
    network_title = reactive.Value("Interaction Network")

    @reactive.Effect
    @reactive.event(input.pickle_file, input.config_file)
    def initialize_filter():
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
            master_instance = CSVFilter(pickle=pickle_path, cfg=config_path)
            csv.set(master_instance)
            ui.notification_show("Files loaded successfully!", duration=3000,
                                 type="message")
            ui.update_navs("plot_tabs", selected="Sele1 vs Sele2")
        except Exception as e:
            csv.set(None)
            filtered_idx.set(None)
            csv_filtered.set(None)
            ui.notification_show(f"Error loading files: {str(e)}",
                                 duration=5000, type="error")

    @reactive.Effect
    @reactive.event(input.transpose_button)
    def handle_transpose():
        transpose_state.set(not transpose_state.get())
        session.send_custom_message("update-transpose-button",
                                    {"active": transpose_state.get()})
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
    @reactive.event(input.organization_method)
    def update_organization_method():
        organization_method.set(input.organization_method())
        if filtered_idx.get() is not None:
            session.send_custom_message("refresh-plots", {})

    @reactive.Effect
    @reactive.event(input.plot_button, input.selected_interactions,
                    input.selected_annotations, input.prevalence_threshold)
    def update_filtered_idx():
        if csv.get() is None:
            filtered_idx.set(None)
            csv_filtered.set(None)
            return
        master_instance = csv.get()
        try:
            mda_idx, mda_status = master_instance.by_mda(
                input.mda_selection() if input.mda_selection() else 'all')
            prev_idx, prev_status = master_instance.by_prevalence(
                input.prevalence_threshold())
            inter_idx, inter_status = master_instance.by_inters(tuple(
                input.selected_interactions()) if input.selected_interactions() else input.selected_interactions())
            annot_idx, annot_status = master_instance.by_notes(tuple(
                input.selected_annotations()) if input.selected_annotations() else input.selected_annotations())

            if mda_status == -1:
                ui.notification_show(
                    "No information for this MDAnalysis selection",
                    type="message", duration=5)
            elif prev_status == -1:
                ui.notification_show(
                    "No interactions with this prevalence or higher.",
                    type="message", duration=5)
            elif inter_status == -1:
                ui.notification_show(
                    "No information for this Interaction selection",
                    type="message", duration=5)
            elif annot_status == -1:
                ui.notification_show(
                    "No information for this Annotation selection",
                    type="message", duration=5)

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
            filtered_idx.set(None)
            csv_filtered.set(None)
            ui.notification_show(f"Error updating filters: {str(e)}",
                                 duration=5000, type="error")

    @reactive.Effect
    @reactive.event(input.save_plot_trigger)
    def handle_download():
        if csv_filtered.get() is None:
            ui.notification_show(
                "No plot to save. Please generate a plot first.",
                duration=3000, type="warning")
            return
        try:
            root = tk.Tk()
            root.withdraw()
            save_dir = filedialog.askdirectory(
                title="Select Directory to Save Plot")
            if not save_dir: return
            save_dir = Path(save_dir)
            active_tab = input.plot_tabs()
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            master_instance = csv.get()
            axisx, axisy = get_axis_names(master_instance)

            df_to_use = csv_filtered.get()
            org_method = organization_method.get()

            if active_tab == "Life Time":
                fig = create_lifetime_plot(df_to_use, input.plot_width(),
                                           input.plot_height(), axisx, axisy,
                                           input.show_prevalence())
                if org_method != "default": fig = apply_organization_to_figure(
                    fig, org_method, df_to_use)
                filename = f"lifetime_plot_{timestamp}.html"
                fig.write_html(str(save_dir / filename))

            elif active_tab == "Sele1 vs Sele2":
                fig = create_plot(df_to_use, input.plot_width(),
                                  input.plot_height(), axisx, axisy,
                                  input.show_prevalence())
                if org_method != "default": fig = apply_organization_to_figure(
                    fig, org_method, df_to_use)
                filename = f"sele1_vs_sele2_plot_{timestamp}.html"
                fig.write_html(str(save_dir / filename))

            elif active_tab == "Prevalence":
                fig1 = create_sel1_interactions_plot(df_to_use,
                                                     input.plot_width(),
                                                     input.plot_height() // 1.5,
                                                     axisx, axisy)
                if org_method != "default": fig1 = apply_organization_to_figure(
                    fig1, org_method, df_to_use)
                filename1 = f"selection1_prevalence_{timestamp}.html"
                fig1.write_html(str(save_dir / filename1))

                fig2 = create_sel2_interactions_plot(df_to_use,
                                                     input.plot_width(),
                                                     input.plot_height() // 1.5,
                                                     axisx, axisy)
                if org_method != "default": fig2 = apply_organization_to_figure(
                    fig2, org_method, df_to_use)
                filename2 = f"selection2_prevalence_{timestamp}.html"
                fig2.write_html(str(save_dir / filename2))
                ui.notification_show(
                    f"Plots saved as:\n{filename1}\n{filename2}",
                    duration=5000, type="message")
                return

            elif active_tab == "Time Series":
                fig = create_interactions_over_time_plot(df_to_use,
                                                         input.plot_width(),
                                                         input.plot_height(),
                                                         axisx, axisy)
                if org_method != "default": fig = apply_organization_to_figure(
                    fig, org_method, df_to_use)
                filename = f"time_series_plot_{timestamp}.html"
                fig.write_html(str(save_dir / filename))

            elif active_tab == "Network":
                net = create_network_plot(df_to_use, input.plot_width(),
                                          input.plot_height(), axisx, axisy)
                filename = f"network_plot_{timestamp}.html"
                net.save_graph(str(save_dir / filename))

            ui.notification_show(f"Saved as:\n{filename}", duration=5000,
                                 type="message")

        except Exception as e:
            ui.notification_show(f"Error saving plot: {str(e)}", duration=5000,
                                 type="error")

    @reactive.Effect
    @reactive.event(input.export_csv_trigger)
    def handle_export_csv():
        ui.notification_show("Preparing to export data...", duration=5,
                             type="message")
        if csv_filtered.get() is None:
            ui.notification_show(
                "No data to export. Please generate a plot first.", duration=5,
                type="warning")
            return
        try:
            root = tk.Tk()
            root.withdraw()
            root.attributes("-topmost", True)
            save_path = filedialog.asksaveasfilename(
                title="Save CSV File", defaultextension=".csv",
                filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
                initialfile=f"intermap_data_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            )
            if not save_path: return
            df_to_export = csv_filtered.get()
            if simplify_labels.get():
                df_to_export = process_axis_labels(df_to_export.copy(), True)
            df_to_export.to_csv(save_path, index=False)
            ui.notification_show(
                f"Data exported successfully as:\n{os.path.basename(save_path)}",
                duration=5, type="message")
        except Exception as e:
            ui.notification_show(f"Error exporting CSV: {str(e)}", duration=5,
                                 type="error")
        finally:
            if 'root' in locals(): root.destroy()

    def check_data():
        if csv_filtered.get() is None: return ui.p(
            ERROR_MESSAGES["nothing_selected"])
        return None

    def render_plot(create_plot_func, is_main=False, is_network=False):
        error = check_data()
        if error is not None: return error
        master_instance = csv.get()
        if master_instance is None: return ui.p(
            ERROR_MESSAGES["nothing_selected"])

        active_tab = input.plot_tabs()
        df_to_use = csv_filtered.get()

        if simplify_labels.get():
            df_to_use = process_axis_labels(df_to_use, True)

        x_axis_title, y_axis_title = get_axis_names(master_instance)
        if active_tab in ["Sele1 vs Sele2", "Prevalence"]:
            if custom_x_axis.get() is not None: x_axis_title = custom_x_axis.get()
            if custom_y_axis.get() is not None: y_axis_title = custom_y_axis.get()

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

        plot_args = [df_to_use, input.plot_width(), input.plot_height() if (
                    is_main or is_network) else input.plot_height() // 1.5,
                     x_axis_title, y_axis_title]
        if is_main: plot_args.append(input.show_prevalence())

        if is_network:
            net = create_plot_func(*plot_args)
            if net is None: return ui.p(ERROR_MESSAGES["nothing_selected"])
            unique_id = f"network_{hash(str(csv_filtered.get()))}"
            temp_file = os.path.join(tempfile.gettempdir(),
                                     f"{unique_id}.html")
            net.save_graph(temp_file)
            with open(temp_file, 'r', encoding='utf-8') as f:
                html_content = f.read()
            title_html = f'<h3 style="text-align:center;font-family:Roboto;margin-bottom:20px;">{plot_title}</h3>'
            enhanced_content = html_content.replace('<body>',
                                                    f'<body>{title_html}')
            try:
                os.remove(temp_file)
            except:
                pass
            return ui.tags.div(ui.HTML(enhanced_content), id=unique_id,
                               style=f"width: 100%; height: {input.plot_height()}px; border: none;")
        else:
            fig = create_plot_func(*plot_args)
            if fig is None: return ui.p(ERROR_MESSAGES["nothing_selected"])
            if plot_title:
                fig.update_layout(title={'text': plot_title, 'x': 0.5,
                                         'font': {'family': "Roboto",
                                                  'size': 20}})
            if organization_method.get() != "default":
                fig = apply_organization_to_figure(fig,
                                                   organization_method.get(),
                                                   df_to_use)
            plot_html = fig.to_html(include_plotlyjs="cdn", full_html=False,
                                    config={"responsive": True})
            style = "width:100%; height:800px; border:1px solid white; margin: 10px 0;" if is_main else f"width:100%; height:{input.plot_height() // 2}px;"
            return ui.tags.div(ui.HTML(plot_html), style=style)

    def create_checkbox_group(input_id, data_keys_getter):
        if csv.get() is None: return ui.input_checkbox_group(input_id, "", [])
        choices = list(data_keys_getter(csv.get()))
        return ui.input_checkbox_group(input_id, "", choices, selected=choices)

    @output
    @render.ui
    @reactive.event(input.plot_button, input.transpose_button,
                    input.simplify_axis_labels, input.apply_axis_names,
                    input.apply_plot_titles, input.organization_method)
    def interaction_plot():
        return render_plot(create_plot, is_main=True)

    @output
    @render.ui
    @reactive.event(input.plot_button, input.simplify_axis_labels,
                    input.apply_axis_names, input.apply_plot_titles,
                    input.organization_method)
    def sel1_interactions_plot():
        return render_plot(create_sel1_interactions_plot)

    @output
    @render.ui
    @reactive.event(input.plot_button, input.simplify_axis_labels,
                    input.apply_axis_names, input.apply_plot_titles,
                    input.organization_method)
    def sel2_interactions_plot():
        return render_plot(create_sel2_interactions_plot)

    @output
    @render.ui
    @reactive.event(input.plot_button, input.simplify_axis_labels,
                    input.apply_plot_titles, input.organization_method)
    def interactions_over_time_plot():
        return render_plot(create_interactions_over_time_plot)

    @output
    @render.ui
    @reactive.event(input.plot_button, input.simplify_axis_labels,
                    input.apply_plot_titles, input.organization_method)
    def lifetime_plot():
        error = check_data()
        if error is not None: return error
        master_instance = csv.get()
        if master_instance is None: return ui.p(
            ERROR_MESSAGES["nothing_selected"])
        df_to_use = csv_filtered.get()
        if simplify_labels.get(): df_to_use = process_axis_labels(df_to_use,
                                                                  True)
        axisx, axisy = get_axis_names(master_instance)
        fig = create_lifetime_plot(df_to_use, input.plot_width(),
                                   input.plot_height(), axisx, axisy,
                                   input.show_prevalence())
        if fig is None: return ui.p(ERROR_MESSAGES["nothing_selected"])
        fig.update_layout(title={'text': lifetime_title.get(), 'x': 0.5,
                                 'font': {'family': "Roboto", 'size': 20}})
        if organization_method.get() != "default":
            fig = apply_organization_to_figure(fig, organization_method.get(),
                                               df_to_use)
        plot_html = fig.to_html(include_plotlyjs="cdn", full_html=False,
                                config={"responsive": True})
        return ui.tags.div(ui.HTML(plot_html),
                           style=f"width:100%; height:{input.plot_height()}px;")

    @output
    @render.ui
    @reactive.event(input.plot_button, input.simplify_axis_labels,
                    input.apply_plot_titles, input.apply_network_settings,
                    input.organization_method)
    def network_plot():
        error = check_data()
        if error is not None: return error
        master_instance = csv.get()
        if master_instance is None: return ui.p(
            ERROR_MESSAGES["nothing_selected"])
        df_to_use = csv_filtered.get()
        if simplify_labels.get(): df_to_use = process_axis_labels(df_to_use,
                                                                  True)
        axisx, axisy = get_axis_names(master_instance)
        net = create_network_plot(df_to_use, input.plot_width(),
                                  input.plot_height(), axisx, axisy,
                                  network_params.get())
        if net is None: return ui.p(ERROR_MESSAGES["nothing_selected"])
        unique_id = f"network_{hash(str(csv_filtered.get()))}"
        temp_file = os.path.join(tempfile.gettempdir(), f"{unique_id}.html")
        net.save_graph(temp_file)
        with open(temp_file, 'r', encoding='utf-8') as f:
            html_content = f.read()
        title_html = f'<h3 style="text-align:center;font-family:Roboto;margin-bottom:20px;">{network_title.get()}</h3>'
        enhanced_content = html_content.replace('<body>',
                                                f'<body>{title_html}')
        try:
            os.remove(temp_file)
        except:
            pass
        return ui.tags.div(ui.HTML(enhanced_content), id=unique_id,
                           style=f"width: 100%; height: {input.plot_height()}px; border: none;")

    @output
    @render.ui
    def interaction_checkboxes():
        return create_checkbox_group("selected_interactions",
                                     lambda x: x.inters2df.keys())

    @output
    @render.ui
    def annotation_checkboxes():
        return create_checkbox_group("selected_annotations",
                                     lambda x: x.notes2df.keys())

    @reactive.Effect
    @reactive.event(input.select_all_interactions)
    def update_interaction_selections():
        if csv.get() is None: return
        ui.update_checkbox_group("selected_interactions", selected=list(
            csv.get().inters2df.keys()) if input.select_all_interactions() else [])

    @reactive.Effect
    @reactive.event(input.select_all_annotations)
    def update_annotation_selections():
        if csv.get() is None: return
        ui.update_checkbox_group("selected_annotations", selected=list(
            csv.get().notes2df.keys()) if input.select_all_annotations() else [])

    @reactive.Effect
    @reactive.event(input.simplify_axis_labels)
    def update_simplify_labels():
        simplify_labels.set(input.simplify_axis_labels())
        if filtered_idx.get() is not None: session.send_custom_message(
            "refresh-plots", {})

    @reactive.Effect
    @reactive.event(input.apply_axis_names)
    def update_custom_axes():
        custom_x_axis.set(
            input.custom_x_axis() if input.custom_x_axis() else None)
        custom_y_axis.set(
            input.custom_y_axis() if input.custom_y_axis() else None)
        if filtered_idx.get() is not None: session.send_custom_message(
            "refresh-plots", {})

    @reactive.Effect
    @reactive.event(input.apply_plot_titles)
    def update_plot_titles():
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
        if filtered_idx.get() is not None: session.send_custom_message(
            "refresh-plots", {})

    @reactive.Effect
    @reactive.event(input.apply_network_settings)
    def update_network_params():
        network_params.set({
            'gravity': input.network_gravity(),
            'central_gravity': input.network_central_gravity(),
            'spring_length': input.network_spring_length(),
            'spring_constant': input.network_spring_constant(),
            'avoid_overlap': input.network_avoid_overlap(),
            'stabilization_iterations': input.network_stabilization_iterations(),
            'physics_enabled': input.network_physics_enabled(),
            'min_node_size': input.network_min_node_size(),
            'max_node_size': input.network_max_node_size(),
            'min_edge_width': input.network_min_edge_width(),
            'max_edge_width': input.network_max_edge_width()
        })
        if csv_filtered.get() is not None: session.send_custom_message(
            "refresh-plots", {})


app = App(app_ui, server)
