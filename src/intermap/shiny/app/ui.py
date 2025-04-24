"""
User interface components for the InterMap Visualizations app.
"""
from screeninfo import get_monitors
from shiny import ui

from intermap.shiny.app.css import CSS_STYLES
from intermap.shiny.app.helpers import get_image_base64

for m in get_monitors():
    personal_width = m.width
    personal_height = m.height
from os.path import dirname, join, abspath


current_dir = dirname(abspath(__file__))
app_dir = dirname(current_dir)
proj_dir = dirname(app_dir)

logo_path = join(proj_dir, "shiny", "statics", "Untitled.png")
favicon_path = join(proj_dir, "statics", "favicon-32x32.png")

"""
   UI: create_app_ui
  #####################################
  #          Welcome Section          #
  #####################################
  #                   #   create_file #
  #      create       # input_section #
  #      plots        #---------------#  
  #      section      #      create   # 
  #                   #      filters  #
  #                   #      section  #
  #####################################
  #              Footer               #
  #####################################
"""


def create_app_ui():
    """Create the main UI of the app."""
    return ui.page_fluid(
        ui.tags.head(
            ui.tags.title("InterMap Visualizations"),
            ui.tags.link(rel="icon", type="image/png", href=favicon_path),
            ui.tags.link(rel="stylesheet",
                href="https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&display=swap",
            ),
            ui.tags.script(src="https://cdn.plot.ly/plotly-latest.min.js"),
            ui.tags.style("\n".join(CSS_STYLES.values())),
        ),
        create_welcome_section(),
        ui.row(create_plots_section(),
               create_filters_section()),
        create_footer(),
    )


def create_welcome_section():
    """Create the welcome section of the app."""
    return ui.div({"class": "welcome-section"},
        ui.div({"class": "welcome-text"},
            ui.h1("Welcome to Intervis!", {"class": "welcome-title"}),
            ui.p("InterMap's Visualization Hub",
                {"class": "welcome-subtitle", "style": "font-style: italic;"},),),
        ui.a(ui.img(
                src=get_image_base64(logo_path),
                class_="welcome-image", alt="InterMap Logo",),
            href="https://rglez.github.io/intermap/",
            target="_blank",
        ),
    )


def create_file_input_section():
    """Create the file input section."""
    return ui.div({"class": "file-input-container"},
        ui.div({"class": "file-browse-container"},
            # CSV file input
            ui.input_file("csv_file", "Upload CSV file:",
                accept=[".csv"],
                button_label="Browse CSV",
                placeholder="No CSV file selected"),
            # Topology file input
            ui.input_file("top_file", "Upload Topology file:",
                accept=['.psf', '.pdb', '.ent', '.pqr', '.pdbqt', '.gro', 'top', '.prmtop', '.parm7',
                        '.dms', '.tpr', '.itp', '.mol2', '.data', '.lammpsdump', '.xyz', '.txyz',
                        '.arc', '.gms', '.log', '.config', '.history', '.xml', '.gsd', '.mmtf', '.in'],
                button_label="Browse Topology",
                placeholder="No topology file selected"),
            ),
        ui.div({"class": "mda-selection-container"},
            ui.input_text("mda_selection", "MDAnalysis Selection:",
                placeholder="e.g., resname ALA or protein (press Enter to apply)",
                value="", width="100%"),
            ui.tags.script(
                """
                document.getElementById('mda_selection').addEventListener('keypress', 
                function(e) {
                    if (e.key === 'Enter') {
                        e.preventDefault();
                        Shiny.setInputValue('mda_selection_submit', 
                        Date.now());}});""")),)

def create_filters_section():
    """Create the filters section of the app."""
    return ui.column(3,
        ui.div({"class": "interaction-filter"},
            # Panel 1
            ui.h4("Data Input", style="font-family: Roboto;"),
            create_file_input_section(),
            ui.hr(),

            # Panel 2 and 3
            ui.layout_columns(
                ui.h4("Interactions"),
                ui.h4("Annotations"),),
            ui.layout_columns(
                ui.div({"class": "checkbox-group"},
                ui.output_ui("interaction_checkboxes")),
                ui.div({"class": "checkbox-group"},
                ui.output_ui("annotation_checkboxes")),),
            ui.hr(),

            # Panel 5
            ui.h4("Prevalence"),
            ui.layout_columns(
                ui.input_slider("prevalence_threshold", "",
                                min=0, max=100, value=30, step=1),
                ui.input_switch("show_prevalence", "Show",
                                value=False), col_widths=[9, 3],),
            ui.hr(),

            # Panel 6
            ui.h4("Plot Settings"),
            ui.layout_columns(
                ui.input_numeric("plot_width", "Plot Width:",
                                 value=personal_width * 0.67),
                ui.input_numeric("plot_height", "Plot Height:",
                                 value=personal_height * 0.75),),
            ui.hr(),),)

def create_plots_section():
    """Create the plots section of the app with styled tabbed navigation."""
    return ui.column(9,
        ui.div({"style": "display: flex; flex-direction: column; "
                         "align-items: center; width: 100%;"},
            ui.div({"class": "custom-tabs-container"},
                ui.navset_tab(
                    # Tab 1: Interaction Heatmap
                    ui.nav_panel("Sele1 vs Sele2",
                        ui.div({"style": "width: 100%; max-width: 90%; "
                                         "margin: 20px auto;"},
                            ui.output_ui("interaction_plot"),),),

                    # Tab 2: Ligand Interactions
                    ui.nav_panel("Details on Sele1",
                        ui.div({"style": "width: 100%; max-width: 90%; "
                                         "margin: 20px auto;"},
                            ui.output_ui("ligand_interactions_plot"),),),

                    # Tab 3: Receptor Interactions
                    ui.nav_panel("Details on Sele2",
                        ui.div({"style": "width: 100%; max-width: 90%; "
                                         "margin: 20px auto;"},
                            ui.output_ui("receptor_interactions_plot"),),),
                    # Tab 4: Interactions Over Time
                    ui.nav_panel("Time Series",
                        ui.div({"style": "width: 100%; max-width: 90%; "
                                         "margin: 20px auto;"},
                            ui.output_ui("interactions_over_time_plot"),),),
                    id="plot_tabs",
                    selected="Interaction Matrix",),),),)


def create_footer():
    """Create a footer with documentation and GitHub links."""
    return ui.div(
        {
            "style": """
                position: fixed;
                bottom: 0;
                left: 0;
                right: 0;
                background-color: black;
                padding: 10px 20px;
                border-top: 1px solid #dee2e6;
                display: flex;
                justify-content: space-between;
                align-items: center;
                font-family: 'Roboto';
                z-index: 1000;
            """
        },
        # First section - Intervis Tutorials
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
                        flex: 1;
                        justify-content: center;
                    """
                },
                ui.tags.i({"class": "fas fa-fingerprint"}),  # Fingerprint icon
                "Intervis Tutorials",
            ),
            href="https://rglez.github.io/intermap/tutorials",
            target="_blank",
            style="text-decoration: none;",
        ),
        # Separator
        ui.div("|", {"style": "color: white; margin: 0 15px;"}),

        # Second section - Documentation
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
                        flex: 1;
                        justify-content: center;
                    """
                },
                ui.tags.i({"class": "fas fa-book-open"}),  # Open book icon
                "InterMap's Documentation",
            ),
            href="https://rglez.github.io/intermap/",
            target="_blank",
            style="text-decoration: none;",
        ),
        # Separator
        ui.div("|", {"style": "color: white; margin: 0 15px;"}),

        # Third section - Intermap Paper
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
                        flex: 1;
                        justify-content: center;
                    """
                },
                ui.tags.i({"class": "fas fa-scroll"}),  # Paper scroll icon
                "Intermap Paper",
            ),
            href="https://scholar.google.com/citations?user=kSAc11cAAAAJ&hl=es&oi=ao",
            target="_blank",
            style="text-decoration: none;",
        ),
        # Separator
        ui.div("|", {"style": "color: white; margin: 0 15px;"}),

        # Fourth section - GitHub
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
                        flex: 1;
                        justify-content: center;
                    """
                },
                ui.tags.i({"class": "fab fa-github"}),  # GitHub icon
                "Delta-Research-Team/intermap",
            ),
            href="https://github.com/Delta-Research-Team/intermap",
            target="_blank",
            style="text-decoration: none;",
        ),
        # Font Awesome
        ui.tags.link(
            rel="stylesheet",
            href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.2/css/all.min.css"
        ),
    )

app_ui = create_app_ui()
