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
            ui.tags.title("InterVis"),
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
    return ui.div({
        "class": "welcome-section",
        "style": "box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);"
    },
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
            ui.input_file("csv_file", "Upload CSV file",
                accept=[".csv"],
                button_label="Browse CSV",
                placeholder="No CSV file selected"),
            # Topology file input
            ui.input_file("top_file", "Upload Topology file",
                accept=['.psf', '.pdb', '.ent', '.pqr', '.pdbqt', '.gro', 'top', '.prmtop', '.parm7',
                        '.dms', '.tpr', '.itp', '.mol2', '.data', '.lammpsdump', '.xyz', '.txyz',
                        '.arc', '.gms', '.log', '.config', '.history', '.xml', '.gsd', '.mmtf', '.in'],
                button_label="Browse Topology",
                placeholder="No topology file selected"),
            ),
        ui.hr(),
        ui.div({"class": "mda-selection-container"},
            ui.input_text("mda_selection", "Atomic Selection (MDAnalysis Syntax)",
                placeholder="e.g., resname ALA or protein",
                value="", width="100%")))


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
                                ui.h4("Annotations"),
                            ),
                            ui.layout_columns(
                                ui.div(
                                    ui.div({"class": "select-all-container"},
                                           ui.input_checkbox(
                                               "select_all_interactions",
                                               "Select All",
                                               value=True
                                           )
                                           ),
                                    ui.div({"class": "checkbox-group"},
                                           ui.output_ui(
                                               "interaction_checkboxes")
                                           )
                                ),
                                ui.div(
                                    ui.div({"class": "select-all-container"},
                                           ui.input_checkbox(
                                               "select_all_annotations",
                                               "Select All",
                                               value=True
                                           )
                                           ),
                                    ui.div({"class": "checkbox-group"},
                                           ui.output_ui(
                                               "annotation_checkboxes")
                                           )
                                ),
                            ),

                            # Panel 5
                            ui.h4("Prevalence"),
                            ui.layout_columns(
                                ui.input_slider("prevalence_threshold", "",
                                                min=0, max=100, value=30,
                                                step=1),
                                ui.input_switch("show_prevalence", "Show",
                                                value=False),
                                col_widths=[9, 3], ),
                            ui.hr(),

                            # Panel 6
                            ui.h4("Plot Settings"),
                            ui.layout_columns(
                                ui.input_numeric("plot_width", "Plot Width",
                                                 value=personal_width * 0.67),
                                ui.input_numeric("plot_height", "Plot Height",
                                                 value=personal_height * 0.75), ),
                            ui.hr(),

                            # Add Sort Options
                            ui.h4("Sort Options"),
                            ui.input_radio_buttons(
                                "sort_by",
                                "",
                                choices={
                                    "name": "Sort by Name",
                                    "number": "Sort by Number",
                                    "annotation": "Sort by Annotation"
                                },
                                selected="name"
                            ),

                            ui.div(
                                {"class": "transpose-button-container"},
                                ui.input_action_button(
                                    "transpose_button",
                                    ui.div(
                                        {"class": "transpose-button-content"},
                                        ui.tags.i(
                                            {"class": "fas fa-exchange-alt"}),
                                        "Transpose Axes"
                                    ),
                                    class_="transpose-button"
                                )
                            ),
                            ui.hr(),


                            ui.tags.style("""
               .transpose-button-container {
                   padding: 10px 15px;
                   margin-bottom: 15px;
               }

               .transpose-button {
                   width: 100%;
                   background-color: #4d4d4dd0;
                   color: white;
                   padding: 12px 20px;
                   border: none;
                   border-radius: 4px;
                   cursor: pointer;
                   transition: all 0.3s ease;
                   display: flex;
                   align-items: center;
                   justify-content: center;
                   font-family: Roboto;
                   font-size: 14px;
               }

               .transpose-button:hover {
                   background-color: #4051b5ff;
                   transform: translateY(-2px);
                   box-shadow: 0 2px 5px rgba(0,0,0,0.2);
               }

               .transpose-button:active {
                   transform: translateY(0);
               }

               .transpose-button-content {
                   display: flex;
                   align-items: center;
                   gap: 8px;
               }

               .transpose-button.active {
                   background-color: #4051b5ff;
               }

               .fa-exchange-alt {
                   transition: transform 0.3s ease;
               }

               .transpose-button.active .fa-exchange-alt {
                   transform: rotate(90deg);
               }
           """),

                            # Script para manejar la actualización del botón
                            ui.tags.script("""
               $(document).ready(function() {
                   Shiny.addCustomMessageHandler('update-transpose-button', function(message) {
                       const button = document.querySelector('.transpose-button');
                       if (message.active) {
                           button.classList.add('active');
                       } else {
                           button.classList.remove('active');
                       }
                   });

                   Shiny.addCustomMessageHandler('refresh-plots', function(message) {
                       var plotTabs = ['interaction_plot', 'ligand_interactions_plot', 
                                     'receptor_interactions_plot', 'interactions_over_time_plot'];
                       plotTabs.forEach(function(plotId) {
                           if (Shiny.shinyapp.$bindings[plotId]) {
                               Shiny.shinyapp.$bindings[plotId].invalidate();
                           }
                       });
                   });
               });
           """),

                            ui.div(
                                {"style": "padding: 0 15px;"},
                                ui.input_action_button(
                                    "plot_button",
                                    "PLOT",
                                    width="100%",
                                    style="""
                                      background-color: #4a4a4a;
                                      color: white;
                                      padding: 20px 35px;  
                                      font-family: Roboto;
                                      font-size: 16px;
                                      border: none;
                                      border-radius: 4px;
                                      cursor: pointer;
                                      transition: background-color 0.3s ease;
                                      margin: 10px 0;
                                  """
                                )
                            ),
                            # Botón de descarga
                            ui.div(
                                {"style": "padding: 0 15px;"},
                                ui.input_action_button(
                                    "download_plot_button",
                                    "SAVE PLOT",
                                    width="100%",
                                    style="""
                              background-color: #4051b5ff;
                              color: white;
                              padding: 20px 35px;
                              font-family: Roboto;
                              font-size: 16px;
                              border: none;
                              border-radius: 4px;
                              cursor: pointer;
                              transition: background-color 0.3s ease;
                              margin: 10px 0;
                          """
                                )
                            ),
                            # Actualizar el CSS
                            ui.tags.style("""
                      #plot_button:hover {
                          background-color: #4051b5ff !important;
                      }
                      #download_plot_button:hover {
                          background-color: #4a4a4a !important;
                      }
                  """),
                            ))



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

                    # Tab 2:
                    ui.nav_panel("Prevalence",
                                 ui.div(
                                     {"style": "width: 100%; max-width: 90%; "
                                               "margin: 20px auto;"},
                                     ui.div({
                                                "style": "display: flex; flex-direction: column; gap: 50px;"},
                                            ui.div({
                                                       "style": "width: 100%; margin-bottom: 40px; padding: 20px; "
                                                                "border-bottom: 1px solid #e0e0e0;"},
                                                   ui.div({
                                                              "style": "height: 400px;"},
                                                          ui.output_ui(
                                                              "ligand_interactions_plot"),
                                                          ),
                                                   ),
                                            ui.div({
                                                       "style": "width: 100%; padding: 20px;"},
                                                   ui.div({
                                                              "style": "height: 400px;"},
                                                          ui.output_ui(
                                                              "receptor_interactions_plot"),
                                                          ),
                                                   ),
                                            ),
                                     ),
                                 ),

                    # Tab 3: Interactions Over Time (antes Tab 4)
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
                        font-size: 18px;
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
                        font-size: 18px;
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
                        font-size: 18px;
                    """
                },
                ui.tags.i({"class": "fas fa-scroll"}),  # Paper scroll icon
                "InterMap's Paper",
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
                        font-size: 18px;
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
