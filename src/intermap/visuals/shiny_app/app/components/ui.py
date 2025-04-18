"""
User interface components for the InterMap Visualizations app.

"""

from shiny import ui
from ..config import CSS_STYLES
from ..utils.helpers import get_image_base64


for m in get_monitors():
    personal_width = m.width
    personal_height = m.height


def create_file_input_section():
    """Create the file input section with topology indicator."""
    return ui.div(
        {"class": "file-input-container"},
        ui.div(
            {"class": "file-browse-container"},
            ui.input_file(
                "file",
                "Upload CSV file:",
                accept=[".csv"],
                button_label="Browse",
                placeholder="No file selected"
            ),
            ui.div(
                {"class": "topology-indicator", "id": "topology-indicator"},
                ui.tags.i({"class": "fas fa-circle"}),
                {"title": "No topology file found"}
            )
        ),
        ui.div(
            {"class": "mda-selection-container"},
            ui.input_text(
                "mda_selection",
                "MDAnalysis Selection:",
                placeholder="e.g., resname ALA or protein",
                value="",
                width="100%"
            )
        )
    )


def create_filters_section():
    """Create the filters section of the app."""
    return ui.column(
        3,
        ui.div(
            {"class": "interaction-filter"},
            # File Input Section
            ui.h4("Data Input"),
            create_file_input_section(),
            ui.input_action_button(
                "show_plots",
                "Show",
                style="margin-top: 10px; width: 100%; font-family: Ubuntu Mono;",
                disabled=True
            ),
            ui.hr(),

            # Prevalence Filter
            ui.h4("Prevalence Filter"),
            ui.input_slider(
                "prevalence_threshold",
                "",
                min=0,
                max=100,
                value=30,
                step=1
            ),
            ui.input_switch(
                "show_prevalence",
                "Show Prevalence Values",
                value=False
            ),
            ui.hr(),

            # Residue Filter
            ui.h4("Filter"),
            ui.div(
                {"class": "search-container"},
                ui.input_text(
                    "residue_filter",
                    "",
                    placeholder="Enter residue name ...",
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

            # Interaction Filter
            ui.h4("Interaction Type Filter"),
            ui.div(
                {"class": "checkbox-group"},
                ui.output_ui("interaction_checkboxes")
            ),
            ui.hr(),

            # Plot Settings
            ui.h4("Plot Settings"),
            ui.layout_columns(
                ui.input_numeric(
                    "plot_width",
                    "Plot Width:",
                    value=personal_width*0.67
                ),
                ui.input_numeric(
                    "plot_height",
                    "Plot Height:",
                    value=personal_height*0.75
                )
            ),
            ui.hr(),
        )
    )


def create_welcome_section():
    """Create the welcome section of the app."""
    return ui.div(
        {"class": "welcome-section"},
        ui.div(
            {"class": "welcome-text"},
            ui.h1("Welcome to InterMap Visualizations!",
                  {"class": "welcome-title"}),
            ui.p(
                "Revolutionizing Molecular Interaction Analysis with High-Speed Computation",
                {"class": "welcome-subtitle", "style": "font-style: italic;"})
        ),
        ui.img(
            {"src": get_image_base64(
                "/home/fajardo01/03_Fajardo_Hub/02_InterMap/visualizations/statics/image/Untitled.png"),
             "class": "welcome-image",
             "alt": "InterMap Logo"}
        )
    )


def create_plots_section():
    """Create the plots section of the app with styled tabbed navigation."""
    return ui.column(
        9,
        ui.div(
            {
                "style": "display: flex; flex-direction: column; align-items: center; width: 100%;"
            },
            ui.div(
                {"class": "custom-tabs-container"},
                ui.navset_tab(
                    # Tab 1: Interaction Heatmap
                    ui.nav_panel(
                        "Interaction Matrix",
                        ui.div(
                            {
                                "style": "width: 100%; max-width: 90%; margin: 20px auto;"},
                            ui.h4(
                                "Interaction Heatmap",
                                style="background-color: #f5f5f5; padding: 10px; border-radius: 4px; text-align: center; font-family: 'Ubuntu Mono';"
                            ),
                            ui.output_ui("interaction_plot")
                        )
                    ),

                    # Tab 2: Ligand Interactions
                    ui.nav_panel(
                        "Ligand Analysis",
                        ui.div(
                            {
                                "style": "width: 100%; max-width: 90%; margin: 20px auto;"},
                            ui.h4(
                                "Ligand Atoms Interaction Analysis",
                                style="background-color: #f5f5f5; padding: 10px; border-radius: 4px; text-align: center; font-family: 'Ubuntu Mono';"
                            ),
                            ui.output_ui("ligand_interactions_plot")
                        )
                    ),

                    # Tab 3: Receptor Interactions
                    ui.nav_panel(
                        "Receptor Analysis",
                        ui.div(
                            {
                                "style": "width: 100%; max-width: 90%; margin: 20px auto;"},
                            ui.h4(
                                "Receptor Atoms Interaction Analysis",
                                style="background-color: #f5f5f5; padding: 10px; border-radius: 4px; text-align: center; font-family: 'Ubuntu Mono';"
                            ),
                            ui.output_ui("receptor_interactions_plot")
                        )
                    ),

                    # Tab 4: Interactions Over Time
                    ui.nav_panel(
                        "Time Analysis",
                        ui.div(
                            {
                                "style": "width: 100%; max-width: 90%; margin: 20px auto;"},
                            ui.h4(
                                "Interaction Types Over Time",
                                style="background-color: #f5f5f5; padding: 10px; border-radius: 4px; text-align: center; font-family: 'Ubuntu Mono';"
                            ),
                            ui.output_ui("interactions_over_time_plot")
                        )
                    ),

                    id="plot_tabs",
                    selected="Interaction Matrix"
                )
            )
        )
    )


def create_app_ui():
    """Create the main UI of the app."""
    return ui.page_fluid(
        # Head section with styles and scripts
        ui.tags.head(
            ui.tags.link(
                rel="stylesheet",
                href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"
            ),
            ui.tags.script(src="https://cdn.plot.ly/plotly-latest.min.js"),
            ui.tags.style("\n".join(CSS_STYLES.values())),
            # Add custom JavaScript
            ui.tags.script("""
                window.Shiny.addCustomMessageHandler('updateTopologyIndicator', function(message) {
                    const indicator = document.getElementById('topology-indicator');
                    if (!indicator) return;

                    if (message.hasTopology) {
                        indicator.classList.add('active');
                        indicator.title = 'Topology file found';
                    } else {
                        indicator.classList.remove('active');
                        indicator.title = 'No topology file found';
                    }
                });
            """)
        ),
        # Welcome section
        create_welcome_section(),
        # Main content row
        ui.row(
            create_filters_section(),
            create_plots_section()
        ),
        # Footer section
        create_footer()
    )

def create_footer():
    """Create a footer with documentation and GitHub links."""
    return ui.div(
        {
            "style": """
                position: fixed;
                bottom: 0;
                left: 0;
                right: 0;
                background-color: #f8f9fa;
                padding: 10px 20px;
                border-top: 1px solid #dee2e6;
                display: flex;
                justify-content: space-between;
                align-items: center;
                font-family: 'Ubuntu Mono';
                z-index: 1000;
            """
        },
        # Left side - Documentation link
        ui.a(
            ui.div(
                {
                    "style": """
                        display: flex;
                        align-items: center;
                        gap: 8px;
                        color: #495057;
                        text-decoration: none;
                        transition: color 0.3s ease;
                    """
                },
                ui.tags.i({"class": "fas fa-book"}),  # Book icon
                "Documentation"
            ),
            href="https://rglez.github.io/intermap/",
            target="_blank",
            style="text-decoration: none;"
        ),

        # Right side - GitHub link
        ui.a(
            ui.div(
                {
                    "style": """
                        display: flex;
                        align-items: center;
                        gap: 8px;
                        color: #495057;
                        text-decoration: none;
                        transition: color 0.3s ease;
                    """
                },
                ui.tags.i({"class": "fab fa-github"}),  # GitHub icon
                "/rglez/intermap"
            ),
            href="https://github.com/rglez/intermap",
            target="_blank",
            style="text-decoration: none;"
        )
    )


# Create the UI instance
app_ui = create_app_ui()
