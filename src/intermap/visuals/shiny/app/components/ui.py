"""
User interface components for the InterMap Visualizations app.

"""

from shiny import ui
from screeninfo import get_monitors

from ..css import CSS_STYLES
from ..utils.helpers import get_image_base64
from pathlib import Path


# todo: Change absolute path @l142

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
            ui.h4("Data Input", style="font-family: Roboto;"),
            create_file_input_section(),
            ui.input_action_button(
                "show_plots",
                "Show",
                style="margin-top: 10px; width: 100%; font-family: Roboto;",
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
                    value=personal_width * 0.67
                ),
                ui.input_numeric(
                    "plot_height",
                    "Plot Height:",
                    value=personal_height * 0.75
                )
            ),
            ui.hr(),
        )
    )

def get_image_path():
    """Get absolute path for static image regardless of OS."""
    return str(Path(__file__).resolve().parent / 'static' / 'image' / 'Untitled.png')

def create_welcome_section():
    """Create the welcome section of the app."""
    return ui.div(
        {"class": "welcome-section"},
        ui.div(
            {"class": "welcome-text"},
            ui.h1("Welcome to Intervis!",
                  {"class": "welcome-title"}),
            ui.p(
                "InterMap: Accelerated Detection of Interaction Fingerprints on Large-Scale Molecular Ensembles",
                {"class": "welcome-subtitle", "style": "font-style: italic;"})
        ),
        ui.img(
            {"src": get_image_base64(get_image_path()),
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
                            ui.output_ui("interaction_plot")
                        )
                    ),

                    # Tab 2: Ligand Interactions
                    ui.nav_panel(
                        "Ligand Analysis",
                        ui.div(
                            {
                                "style": "width: 100%; max-width: 90%; margin: 20px auto;"},
                            ui.output_ui("ligand_interactions_plot")
                        )
                    ),

                    # Tab 3: Receptor Interactions
                    ui.nav_panel(
                        "Receptor Analysis",
                        ui.div(
                            {
                                "style": "width: 100%; max-width: 90%; margin: 20px auto;"},
                            ui.output_ui("receptor_interactions_plot")
                        )
                    ),

                    # Tab 4: Interactions Over Time
                    ui.nav_panel(
                        "Time Analysis",
                        ui.div(
                            {
                                "style": "width: 100%; max-width: 90%; margin: 20px auto;"},
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
                href="https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&display=swap"
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
                background-color: #4051b5ff;
                padding: 10px 20px;
                border-top: 1px solid #dee2e6;
                display: flex;
                justify-content: space-between;
                align-items: center;
                font-family: 'Roboto';
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
                        color: white;
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
                        color: white;
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
