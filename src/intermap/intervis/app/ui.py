"""
User interface components for the InterMap Visualizations app.
"""
from screeninfo import get_monitors
from shiny import ui

from intermap.intervis.app.css import CSS_STYLES
from intermap.intervis.app.helpers import get_image_base64

for m in get_monitors():
    personal_width = m.width
    personal_height = m.height
from os.path import dirname, join, abspath

current_dir = dirname(abspath(__file__))
app_dir = dirname(current_dir)
proj_dir = dirname(app_dir)

logo_path = join(proj_dir, "intervis", "statics", "Untitled.png")
favicon_path = join(proj_dir, "statics", "favicon-32x32.png")

def create_network_controls_panel():
    """Create the network controls panel for Tab 5."""
    return ui.div(
        {
            "id": "network_controls_panel",
            "class": "network-controls-panel",
            "style": "display: none; position: fixed; left: 10px; top: 150px; width: 280px; z-index: 100; background: white; padding: 15px; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1); max-height: 70vh; overflow-y: auto;"
        },
        ui.h5("Network Controls", style="margin-bottom: 15px; color: #4051b5; font-family: Roboto;"),

        # Physics Controls
        ui.div(
            {"class": "control-group", "style": "margin-bottom: 15px;"},
            ui.h6("Physics", style="font-weight: 500; margin-bottom: 8px; color: #333;"),
            ui.input_slider("network_gravity", "Gravitational Constant",
                           min=-500, max=0, value=-200, step=10),
            ui.input_slider("network_central_gravity", "Central Gravity",
                           min=0, max=0.1, value=0.005, step=0.001),
            ui.input_slider("network_spring_length", "Spring Length",
                           min=50, max=500, value=200, step=10),
            ui.input_slider("network_spring_constant", "Spring Constant",
                           min=0.1, max=2.0, value=0.5, step=0.1),
            ui.input_slider("network_avoid_overlap", "Avoid Overlap",
                           min=0, max=1, value=0.8, step=0.1),
        ),

        # Node Controls
        ui.div(
            {"class": "control-group", "style": "margin-bottom: 15px;"},
            ui.h6("Nodes", style="font-weight: 500; margin-bottom: 8px; color: #333;"),
            ui.input_slider("network_min_node_size", "Min Node Size",
                           min=10, max=30, value=20, step=2),
            ui.input_slider("network_max_node_size", "Max Node Size",
                           min=40, max=100, value=50, step=5),
        ),

        # Edge Controls
        ui.div(
            {"class": "control-group", "style": "margin-bottom: 15px;"},
            ui.h6("Edges", style="font-weight: 500; margin-bottom: 8px; color: #333;"),
            ui.input_slider("network_min_edge_width", "Min Edge Width",
                           min=1, max=10, value=5, step=1),
            ui.input_slider("network_max_edge_width", "Max Edge Width",
                           min=10, max=30, value=15, step=1),
        ),

        # Simulation Controls
        ui.div(
            {"class": "control-group"},
            ui.h6("Simulation", style="font-weight: 500; margin-bottom: 8px; color: #333;"),
            ui.input_slider("network_stabilization_iterations", "Stabilization Iterations",
                           min=100, max=2000, value=1000, step=100),
            ui.input_switch("network_physics_enabled", "Physics Enabled", value=True),
        ),

        # Apply button
        ui.div(
            {"style": "margin-top: 15px; text-align: center;"},
            ui.input_action_button("apply_network_settings", "Apply Settings",
                                 class_="search-button",
                                 style="width: 100%;")
        )
    )

def create_app_ui():
    """Create the main UI of the app."""
    return ui.page_fluid(
        ui.tags.head(
            ui.tags.title("InterVis"),
            ui.tags.link(rel="icon", type="image/png", href=favicon_path),
            ui.tags.link(rel="stylesheet",
                         href="https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&display=swap",
                         ),
            ui.tags.script(src="https://cdn.plot.ly/plotly-latest.min.js"),
            ui.tags.style("\n".join(CSS_STYLES.values())),
        ),
        create_welcome_section(),
        create_network_controls_panel(),
        ui.row(create_plots_section(),
               create_filters_section()),
        create_footer(),
    )

def create_welcome_section():
    """Create the welcome section of the app."""
    return ui.div(
        {
            "class": "welcome-section",
            "style": "box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);"
        },
        ui.div(
            {"class": "welcome-text"},
            ui.h1("Welcome to Intervis!", {"class": "welcome-title"}),
            ui.p("InterMap's Visualization App",
                {"class": "welcome-subtitle", "style": "font-style: italic;"}),
        ),
        ui.a(
            ui.img(
                src=get_image_base64(logo_path),
                class_="welcome-image",
                alt="InterMap Logo",
            ),
            href="https://rglez.github.io/intermap/",
            target="_blank",
        ),
    )

def create_file_input_section():
    """Create the file input section."""
    return ui.div(
        {"class": "file-input-container"},
        ui.div(
            {"class": "file-browse-container"},
            # Pickle file input
            ui.input_file("pickle_file", "Upload Pickle file",
                accept=[".pickle"],
                button_label="Browse Pickle",
                placeholder="No pickle file selected"),
            # Config file input
            ui.input_file("config_file", "Upload Config file",
                accept=[".cfg"],
                button_label="Browse Config",
                placeholder="No config file selected"),
        ),
        ui.hr(),
        ui.div(
            {"class": "mda-selection-container"},
            ui.input_text("mda_selection", "Atomic Selection (MDAnalysis Syntax)",
                placeholder="e.g., resname ALA or protein",
                value="", width="100%")
        )
    )

def create_filters_section():
    """Create the filters section of the app."""
    return ui.column(
        3,
        ui.div(
            {"class": "interaction-filter"},
            # Panel 1
            ui.h4("Data Input", style="font-family: Roboto;"),
            create_file_input_section(),
            ui.hr(),

            # Panel 2 and 3
            ui.layout_columns(
                ui.h4("Interactions"),
                ui.h4("Annotations"),
            ),
            ui.layout_columns(
                ui.div(
                    ui.div(
                        {"class": "select-all-container"},
                        ui.input_checkbox(
                            "select_all_interactions",
                            "Select All",
                            value=True
                        )
                    ),
                    ui.div(
                        {"class": "checkbox-group"},
                        ui.output_ui("interaction_checkboxes")
                    )
                ),
                ui.div(
                    ui.div(
                        {"class": "select-all-container"},
                        ui.input_checkbox(
                            "select_all_annotations",
                            "Select All",
                            value=True
                        )
                    ),
                    ui.div(
                        {"class": "checkbox-group"},
                        ui.output_ui("annotation_checkboxes")
                    )
                ),
            ),

            # Panel 4
            ui.h4("Prevalence"),
            ui.layout_columns(
                ui.input_slider("prevalence_threshold", "",
                                min=0, max=100, value=30,
                                step=1),
                ui.input_switch("show_prevalence", "Show",
                                value=False),
                col_widths=[9, 3],
            ),
            ui.hr(),

            # Panel 5
            ui.h4("Plot Settings"),
            ui.layout_columns(
                ui.input_numeric("plot_width", "Width",
                                 value=personal_width * 0.67),
                ui.input_numeric("plot_height", "Height",
                                 value=personal_height * 0.75),
            ),
            ui.hr(),

            # Panel 6 - Data Sorting
            ui.h4("Data Sorting"),
            ui.div(
                {"class": "organization-container",
                 "style": "margin-bottom: 15px;"},
                ui.layout_columns(
                    ui.p("Organize labels by:",
                         style="margin-top: 8px; margin-bottom: 0; font-weight: 500;"),
                    ui.input_radio_buttons(
                        "organization_method",
                        "",
                        choices={
                            "resname": "Residue Name",
                            "resnum": "Residue Number",
                            "annotation": "Annotations"
                        },
                        selected="resname",
                        inline=True,
                        width="100%"
                    ),
                    col_widths=[4, 8]
                )
            ),
            ui.hr(),

            # Panel 7 - Axis Settings section
            ui.h4("Axis Settings and Plot Title"),
            ui.div(
                {"class": "axis-settings-container"},
                # Option to simplify axis labels
                ui.input_checkbox("simplify_axis_labels",
                                  "Show simplified axis labels", value=False),

                # Custom axis names
                ui.input_action_button(
                    "rename_axes_button",
                    "Rename Axes",
                    class_="search-button"
                ),

                # Containers for custom axis names (initially hidden)
                ui.div(
                    {"id": "custom_axes_inputs",
                     "style": "display: none; margin-top: 10px;"},
                    ui.input_text("custom_x_axis", "X-Axis Title",
                                  placeholder="Custom X-Axis name"),
                    ui.input_text("custom_y_axis", "Y-Axis Title",
                                  placeholder="Custom Y-Axis name"),
                    ui.div(
                        {
                            "style": "display: flex; justify-content: flex-end; margin-top: 8px;"},
                        ui.input_action_button("apply_axis_names", "Apply",
                                               class_="search-button")
                    )
                ),
                ui.p(
                    "Note: Custom axis names apply only to Heatmap and Prevalence plots",
                    style="font-size: 12px; font-style: italic; margin-top: 8px; color: #666;"),

                # Custom plot titles
                ui.input_action_button(
                    "customize_plot_titles_button",
                    "Customize Plot Titles",
                    class_="search-button"
                ),

                # Containers for custom plot titles (initially hidden)
                ui.div(
                    {"id": "custom_plot_titles_inputs",
                     "style": "display: none; margin-top: 10px;"},
                    ui.input_text("heatmap_plot_title", "Heatmap Plot Title",
                                  value="Interaction Heatmap"),
                    ui.input_text("prevalence_plot_title",
                                  "Prevalence Plot Title",
                                  value="Interaction Prevalence Analysis"),
                    ui.input_text("lifetime_plot_title", "Lifetime Plot Title",
                                  value="Interaction Lifetimes Distribution"),
                    ui.input_text("timeseries_plot_title",
                                  "Time Series Plot Title",
                                  value="Interactions Over Time"),
                    ui.input_text("network_plot_title", "Network Plot Title",
                                  value="Interaction Network"),
                    ui.div(
                        {
                            "style": "display: flex; justify-content: flex-end; margin-top: 8px;"},
                        ui.input_action_button("apply_plot_titles", "Apply",
                                               class_="search-button")
                    )
                ),
            ),
            # Transpose button
            ui.div(
                {"class": "transpose-button-container"},
                ui.input_action_button(
                    "transpose_button",
                    ui.div(
                        {"class": "transpose-button-content"},
                        ui.tags.i({"class": "fas fa-exchange-alt"}),
                        "Transpose Axes"
                    ),
                    class_="transpose-button"
                )
            ),
            ui.hr(),

            ui.div({"id": "directory-picker-container"}),

            # JavaScript handlers
            ui.tags.script("""
                $(document).ready(function() {
                    // Handlers for transpose button
                    Shiny.addCustomMessageHandler('update-transpose-button', function(message) {
                        const button = document.querySelector('.transpose-button');
                        if (message.active) {
                            button.classList.add('active');
                        } else {
                            button.classList.remove('active');
                        }
                    });

                    // Handler for refreshing plots
                    Shiny.addCustomMessageHandler('refresh-plots', function(message) {
                        var plotTabs = ['interaction_plot', 'sel1_interactions_plot', 
                                        'sel2_interactions_plot', 'interactions_over_time_plot', 
                                        'lifetime_plot', 'network_plot'];
                        plotTabs.forEach(function(plotId) {
                            if (Shiny.shinyapp.$bindings[plotId]) {
                                Shiny.shinyapp.$bindings[plotId].invalidate();
                            }
                        });
                    });

                    // Network controls panel visibility
                    $(document).on('click', 'a[data-value]', function() {
                        var tabValue = $(this).data('value');
                        var networkPanel = $('#network_controls_panel');
                        
                        if (tabValue === 'Network') {
                            networkPanel.show();
                        } else {
                            networkPanel.hide();
                        }
                    });

                    // Toggle custom axis inputs
                    $('#rename_axes_button').on('click', function() {
                        $('#custom_axes_inputs').toggle();
                    });

                    // Apply button handler for axis names
                    $('#apply_axis_names').on('click', function() {
                        Shiny.setInputValue('custom_axes_applied', Date.now());
                        $('#custom_axes_inputs').hide();
                    });

                    // Toggle custom plot title inputs
                    $('#customize_plot_titles_button').on('click', function() {
                        $('#custom_plot_titles_inputs').toggle();
                    });

                    // Apply button handler for plot titles
                    $('#apply_plot_titles').on('click', function() {
                        Shiny.setInputValue('custom_plot_titles_applied', Date.now());
                        $('#custom_plot_titles_inputs').hide();
                    });

                    // Button click handlers
                    $('#download_plot_button').on('click', function() {
                        Shiny.setInputValue('save_plot_trigger', Date.now());
                    });

                    $('#export_csv_button').on('click', function() {
                        console.log("Export CSV button clicked");
                        Shiny.setInputValue('export_csv_trigger', Date.now());
                    });
                });
            """),

            # Action buttons
            ui.div(
                {"class": "action-buttons-container"},
                ui.input_action_button("plot_button", "PLOT"),
                ui.input_action_button("download_plot_button", "SAVE PLOT"),
                ui.input_action_button("export_csv_button", "EXPORT CSV",
                                       class_="export-csv-button"),
            )
        )
    )

def create_plots_section():
    """Create the plots section of the app with styled tabbed navigation."""
    return ui.column(
        9,
        ui.div(
            {"class": "plots-main-container"},
            ui.div(
                {"class": "custom-tabs-container"},
                ui.navset_tab(
                    # Tab 1: Interaction Heatmap
                    ui.nav_panel(
                        "Sele1 vs Sele2",
                        ui.div(
                            {"class": "plot-tab-content"},
                            ui.output_ui("interaction_plot"),
                        ),
                    ),

                    # Tab 2: Prevalence
                    ui.nav_panel(
                        "Prevalence",
                        ui.div(
                            {"class": "prevalence-container"},
                            ui.div(
                                {"class": "prevalence-plots"},
                                ui.div(
                                    {"class": "prevalence-plot-item"},
                                    ui.h3("", {"class": "prevalence-plot-title"}),
                                    ui.div(
                                        {"class": "prevalence-plot-content"},
                                        ui.output_ui("sel1_interactions_plot"),
                                    ),
                                ),
                                ui.div(
                                    {"class": "prevalence-plot-item"},
                                    ui.h3("", {"class": "prevalence-plot-title"}),
                                    ui.div(
                                        {"class": "prevalence-plot-content"},
                                        ui.output_ui("sel2_interactions_plot"),
                                    ),
                                ),
                            ),
                        ),
                    ),

                    # Tab 3: Life Time
                    ui.nav_panel(
                        "Life Time",
                        ui.div(
                            {"class": "plot-tab-content"},
                            ui.output_ui("lifetime_plot"),
                        ),
                    ),

                    # Tab 4: Time Series
                    ui.nav_panel(
                        "Time Series",
                        ui.div(
                            {"class": "plot-tab-content"},
                            ui.output_ui("interactions_over_time_plot"),
                        ),
                    ),

                    # Tab 5: Network
                    ui.nav_panel(
                        "Network",
                        ui.div(
                            {"class": "plot-tab-content",
                             "style": "width: 80%; max-width: 80%; margin: 20px auto; margin-left: 320px;"},
                            ui.output_ui("network_plot")
                        )
                    ),

                    id="plot_tabs",
                    selected="Interaction Matrix",
                ),
            ),
        ),
        ui.tags.script("""
            $(document).ready(function() {
                $(document).on('click', 'a[data-value="Network"]', function() {
                    setTimeout(function() {
                        $("#plot_button").click();
                }, 100);
                });
            });
        """)
    )

def create_footer():
    """Create a footer with documentation and GitHub links."""
    return ui.div(
        {"class": "footer-container"},
        # First section - Docs & Tutorials
        ui.a(
            ui.div(
                {"class": "footer-link"},
                ui.tags.i({"class": "fas fa-book-open"}),
                "Docs & Tutorials",
            ),
            href="https://rglez.github.io/intermap/",
            target="_blank",
            style="text-decoration: none;",
        ),

        ui.div("|", {"class": "footer-divider"}),

        # Second section - InterMap's Paper
        ui.a(
            ui.div(
                {"class": "footer-link"},
                ui.tags.i({"class": "fas fa-scroll"}),
                "InterMap's Paper",
            ),
            href="https://scholar.google.com/citations?user=kSAc11cAAAAJ&hl=es&oi=ao",
            target="_blank",
            style="text-decoration: none;",
        ),

        ui.div("|", {"class": "footer-divider"}),

        # Third section - Delta-Research
        ui.a(
            ui.div(
                {"class": "footer-link"},
                ui.tags.i({"class": "fab fa-github"}),
                "Delta-Research-Team/intermap",
            ),
            href="https://github.com/Delta-Research-Team/intermap",
            target="_blank",
            style="text-decoration: none;",
        ),
        ui.tags.link(
            rel="stylesheet",
            href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.2/css/all.min.css"
        ),
    )

app_ui = create_app_ui()
