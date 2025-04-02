"""
Main application file for InterMap Visualizations.
Integrates UI components with server-side logic.

"""

import pandas as pd
from shiny import App, render, reactive, ui, experimental
import os

from .components.ui import app_ui
from .components.plots import (
    create_plot,
    create_ligand_interactions_plot,
    create_receptor_interactions_plot,
    create_interactions_over_time_plot,
    initialize_search_state
)
from .utils.helpers import (
    generate_interaction_choices,
    find_topology_file,
    validate_mda_selection
)
from .components.sele_shiny import TSVFilter
from .config import all_interactions_colors, ERROR_MESSAGES

def server(input, output, session):
    """Define server logic for the InterMap Visualizations app."""

    data_store = reactive.Value({'df': None})
    search_state = reactive.Value({'search_term': '', 'active': False})
    show_state = reactive.Value(False)
    topology_state = reactive.Value({'has_topology': False, 'path': None})

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
            ui.notification_show(
                f"Error: {str(e)}",
                type="error",
                duration=None
            )

    @reactive.Effect
    @reactive.event(input.file)
    async def check_topology():
        """Check for topology file when TSV is uploaded."""
        file_infos = input.file()
        if not file_infos:
            topology_state.set({'has_topology': False, 'path': None})
            await session.send_custom_message(
                'updateTopologyIndicator',
                {'hasTopology': False}
            )
            return

        temp_tsv_path = file_infos[0]['datapath']

        original_filename = file_infos[0]['name']

        current_dir = os.getcwd()

        original_path = os.path.join(current_dir, original_filename)


        has_topo, topo_path = find_topology_file(temp_tsv_path, original_path)
        topology_state.set({'has_topology': has_topo, 'path': topo_path})

        await session.send_custom_message(
            'updateTopologyIndicator',
            {'hasTopology': has_topo}
        )

        if has_topo:
            ui.notification_show(
                "Topology file found!",
                type="message",
                duration=3
            )
        else:
            ui.notification_show(
                ERROR_MESSAGES['no_topology'],
                type="warning",
                duration=3
            )

    @reactive.Effect
    @reactive.event(input.file)
    def validate_file():
        """Validate file and enable/disable show button."""
        file_infos = input.file()
        if file_infos:
            file_name = file_infos[0]['name']
            if not file_name.lower().endswith('.tsv'):
                ui.notification_show(
                    ERROR_MESSAGES['invalid_file'],
                    type="error"
                )
                ui.update_action_button("show_plots", disabled=True)
                return
            ui.update_action_button("show_plots", disabled=False)
        else:
            ui.update_action_button("show_plots", disabled=True)
        show_state.set(False)

    @reactive.Effect
    @reactive.event(input.show_plots)
    def handle_show():
        """Handle show button clicks and validate MDAnalysis selection if present."""
        if topology_state.get()['has_topology'] and input.mda_selection():
            is_valid, error_msg = validate_mda_selection(input.mda_selection())
            if not is_valid:
                ui.notification_show(
                    ERROR_MESSAGES['invalid_selection'].format(error_msg),
                    type="error"
                )
                return
        show_state.set(True)

    @reactive.Effect
    @reactive.event(input.search_button)
    def handle_search():
        """Update search state when search button is clicked."""
        search_state.set({
            'search_term': input.residue_filter(),
            'active': True
        })

    @output
    @render.ui
    def residue_not_found():
        """Show error message when residue is not found."""
        if not search_state.get()['active'] or data_store.get()['df'] is None:
            return ui.div()

        df = data_store.get()['df']
        search_term = search_state.get()['search_term'].upper()

        if not search_term:
            return ui.div()

        sel1_matches = df['sel1_atom'].str.contains(search_term, case=False).any()
        sel2_matches = df['sel2_atom'].str.contains(search_term, case=False).any()

        if not sel1_matches and not sel2_matches:
            return ui.p("Not Found", style="color: red; margin-top: 5px;")
        return ui.div()

    @output
    @render.ui
    @reactive.event(input.file, input.plot_width, input.plot_height,
                   input.selected_interactions, input.prevalence_threshold,
                   input.show_prevalence, input.search_button, input.show_plots,
                   input.mda_selection)
    def interaction_plot():
        """Render the main interaction heatmap plot."""
        if not show_state():
            return None

        file_infos = input.file()
        if not file_infos:
            return ui.p(
                ERROR_MESSAGES['no_file'],
                style="font-family: Ubuntu Mono; font-style: italic;"
            )

        try:
            if data_store.get()['df'] is None:
                tsv_path = file_infos[0]['datapath']
                selection = input.mda_selection()
                topo_state = topology_state.get()

                if topo_state['has_topology'] and selection:
                    # Validate MDAnalysis selection
                    is_valid, error_msg = validate_mda_selection(selection)
                    if not is_valid:
                        ui.notification_show(
                            ERROR_MESSAGES['invalid_selection'].format(error_msg),
                            type="error"
                        )
                        return ui.p(ERROR_MESSAGES['invalid_selection'].format(error_msg))

                    try:
                        # Use TSVFilter if topology and selection are present
                        filter_obj = TSVFilter(topo_state['path'], selection, tsv_path)
                        df = filter_obj.df
                    except Exception as e:
                        ui.notification_show(str(e), type="error")
                        return ui.p(f"Error: {str(e)}")
                else:
                    # Normal TSV reading
                    try:
                        df = pd.read_csv(tsv_path, sep=None, engine='python')
                    except:
                        try:
                            df = pd.read_csv(tsv_path, sep='\t', encoding='utf-8')
                        except:
                            try:
                                df = pd.read_csv(tsv_path, sep='\t', encoding='latin-1')
                            except Exception as e:
                                print(f"Error reading file: {str(e)}")
                                raise

                # Clean column names
                df.columns = df.columns.str.strip()

                # Verify required columns
                required_columns = ['sel1_atom', 'sel2_atom', 'interaction_name', 'prevalence', 'timeseries']
                missing_columns = [col for col in required_columns if col not in df.columns]

                if missing_columns:
                    raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

                # Process timeseries
                frames = df['timeseries'].iloc[0].strip() if not df.empty else ''
                num_frames = len(frames)

                # Create frame columns
                for i in range(num_frames):
                    df[str(i)] = df['timeseries'].str.strip().apply(lambda x: int(x[i]))

                data_store.set({'df': df})
            else:
                df = data_store.get()['df']

            fig = create_plot(
                df,
                input.plot_width(),
                input.plot_height(),
                input.selected_interactions(),
                input.prevalence_threshold(),
                input.show_prevalence()
            )

            if fig is None:
                return ui.p(ERROR_MESSAGES['no_interactions'])

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
            error_message = ERROR_MESSAGES['plot_error'].format(str(e))
            return ui.div(
                ui.p(
                    error_message,
                    style="color:red; font-family: Ubuntu Mono;"
                ),
                ui.pre(
                    traceback.format_exc(),
                    style="font-family: Ubuntu Mono; font-size: 12px; background-color: #f8f9fa; padding: 10px; border-radius: 4px;"
                )
            )

    @output
    @render.ui
    @reactive.event(input.file, input.plot_width, input.plot_height,
                   input.selected_interactions, input.prevalence_threshold,
                   input.search_button, input.show_plots)
    def ligand_interactions_plot():
        """Render the ligand interactions plot."""
        if not show_state():
            return None

        data = data_store.get()
        if data['df'] is None:
            return ui.p(ERROR_MESSAGES['no_file'])

        try:
            fig = create_ligand_interactions_plot(
                data['df'],
                input.plot_width(),
                input.plot_height() // 2,
                input.selected_interactions(),
                input.prevalence_threshold()
            )

            if fig is None:
                return ui.p(
                    "Loading ..."
                )

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
            error_message = ERROR_MESSAGES['plot_error'].format(str(e))
            return ui.div(
                ui.p(
                    error_message,
                    style="color:red; font-family: Ubuntu Mono;"
                ),
                ui.pre(
                    traceback.format_exc(),
                    style="font-family: Ubuntu Mono; font-size: 12px;"
                )
            )

    @output
    @render.ui
    @reactive.event(input.file, input.plot_width, input.plot_height,
                   input.selected_interactions, input.prevalence_threshold,
                   input.search_button, input.show_plots)
    def receptor_interactions_plot():
        """Render the receptor interactions plot."""
        if not show_state():
            return None

        data = data_store.get()
        if data['df'] is None:
            return ui.p(ERROR_MESSAGES['no_file'])

        try:
            fig = create_receptor_interactions_plot(
                data['df'],
                input.plot_width(),
                input.plot_height() // 2,
                input.selected_interactions(),
                input.prevalence_threshold()
            )

            if fig is None:
                return ui.p(
                    "Loading ..."
                )

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
            error_message = ERROR_MESSAGES['plot_error'].format(str(e))
            return ui.div(
                ui.p(
                    error_message,
                    style="color:red; font-family: Ubuntu Mono;"
                ),
                ui.pre(
                    traceback.format_exc(),
                    style="font-family: Ubuntu Mono; font-size: 12px;"
                )
            )

    @output
    @render.ui
    @reactive.event(input.file, input.plot_width, input.plot_height,
                   input.selected_interactions, input.prevalence_threshold,
                   input.search_button, input.show_plots)
    def interactions_over_time_plot():
        """Render the interactions over time plot."""
        if not show_state():
            return None

        data = data_store.get()
        if data['df'] is None:
            return ui.p(ERROR_MESSAGES['no_file'])

        try:
            fig = create_interactions_over_time_plot(
                data['df'],
                input.plot_width(),
                input.plot_height() // 2,
                input.selected_interactions(),
                input.prevalence_threshold()
            )

            if fig is None:
                return ui.p(
                    "Loading ..."
                )

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
            error_message = ERROR_MESSAGES['plot_error'].format(str(e))
            return ui.div(
                ui.p(
                    error_message,
                    style="color:red; font-family: Ubuntu Mono;"
                ),
                ui.pre(
                    traceback.format_exc(),
                    style="font-family: Ubuntu Mono; font-size: 12px;"
                )
            )

    @output
    @render.ui
    def interaction_checkboxes():
        """Render the interaction checkboxes."""
        df = data_store.get()['df']
        choices = generate_interaction_choices(df)

        return ui.input_checkbox_group(
            "selected_interactions",
            "",
            choices=choices,
            selected=list(all_interactions_colors.keys())
        )


# Create the Shiny app instance
app = App(app_ui, server)
