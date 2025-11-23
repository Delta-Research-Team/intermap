``InterMap`` comes with a local, interactive Shiny application named ``InterVis``, which allows users to explore the
interactions computed by the tool. This user-friendly web interface simplifies the analysis and interpretation of
molecular interaction data. No internet connection is required to run the app.

## 1. Launching InterVis

Once InterMap is installed and your virtual environment is activated, run the following command from the terminal:

```bash
intervis
```

This will automatically launch the InterVis server and open your default web browser with the app's interface.

![Image title](assets/global_panel.png){width="75%"}

/// caption
InterVis interface is composed of three main areas: the ``Control Panel`` (right side), the ``Visualization Panel`` (
left side), and the ``Footer`` (bottom bar).
///

The InterVis interface consists of three main areas:

- The ``Control Panel`` (right side) contains tools for loading and filtering data. It is always visible for quick
  access and its effects are reflected in the visualizations after clicking the ``Plot`` button.

- The ``Visualization Panel`` (left side) displays different types of visualizations. It is organized into tabs for
  smooth navigation and allows users to explore interaction data interactively.

- The ``Footer`` (bottom bar) provides links to the InterMap's documentation and tutorials, its scientific paper, and
  the GitHub repository.

## 2. Uploading Data

InterVis needs two files before attempting to visualize the interactions: the `.PICKLE` file containing the interaction
data and the `.CFG` file with the configuration parameters used during the InterMap analysis. You can upload these files
using the ``Control Panel`` on the right side of the interface.

- Click on `Browse PICKLE` / `Browse Config` to upload these file. Alternatively, you can drag and drop them into the
  designated area.

- The time Intervis takes to upload the PICKLE and show the plots depends on the amount of data.

!!! warning "Environment Updates"

    The paths in the `config file` leading to the topology and trajectory must be absolute.

## 3. Configuring Visualizations

Once the files are uploaded, you can configure the visualizations that will render in the different tabs of the
``Visualization Panel`` using the following sections under the  ``Control Panel``.

??? note "Atomic Selections"

    - This textbox allows to filter interactions using the [MDAnalysis selection syntax](https://userguide.mdanalysis.org/stable/selections.html). Only atoms/residues that match the selection will be considered in the visualizations.

    !!! example "MDAnalysis Selection Examples"
    ```
    - "protein"                                          # All protein atoms
    - "resname ALA"                                      # All alanine residues
    - "protein and name CA"                              # Alpha carbon atoms in proteins
    - "resid 10 to 30"                                   # Residues with IDs from 10 to 30
    - "(resname ALA or resname VAL) and not name H*"     # Heavy atoms in ALA or VAL residues
    ```

??? note "Interactions"

    Upon uploading the `.PICKLE` file, InterVis will automatically detect all stored interactions and display them in this panel.

    - You can select one or multiple types. Click ``Select All`` to toggle all options.

??? note "Annotations"

    Upon uploading the `.PICKLE` file, InterVis will automatically detect all stored annotations and display them in this panel.
    
    - You can filter the visualized data based on the custom annotations defined in the configuration file.


    - You can select one or multiple annotations. Click ``Select All`` to toggle all options.

??? note "Prevalence Threshold"

    When the files are loaded, InterVis automatically assigns a default prevalence threshold of ``30%``. You can adjust this value to filter out interactions based on their prevalence (0–100%).

    - The ``Show`` option will allow you to show/hide numeric values on plots.


## 4. Plot Settings and Customization

InterVis provides several advanced configuration options to customize your visualizations:

??? note "Plot Dimensions"

    You can customize the size of your plots to fit your screen or presentation needs:
    
    - ``Width``: Set the plot width in pixels (default: automatically adjusted to your screen)
    - ``Height``: Set the plot height in pixels (default: automatically adjusted to your screen)
    
    These settings apply to all visualization tabs and help optimize the viewing experience for different screen sizes or export requirements.

??? note "Data Sorting"

    Control how labels are organized in your visualizations through the ``Organize labels by`` option:
    
    - ``Residue Name``: Sort alphabetically by residue type (ALA, ARG, ASP, etc.)
    - ``Residue Number``: Sort numerically by residue number
    - ``Annotations``: Sort by custom annotations defined in your configuration file
    
    This feature helps you identify patterns more easily by organizing data in the most meaningful way for your analysis.

??? note "Axis Settings and Plot Titles"

    Customize axis labels and plot titles for publication-ready figures:

    ### Simplify Axis Labels
    
    - Enable this checkbox to show simplified axis labels without detailed atom/residue information
    - Useful for creating cleaner visualizations for presentations
    
    ### Rename Axes
    
    1. Click the ``Rename Axes`` button to reveal custom axis input fields
    2. Enter your desired names for X-Axis and Y-Axis titles
    3. Click ``Apply`` to update the axes
    
    !!! info "Scope"
        Custom axis names apply only to Heatmap and Prevalence plots
    
    ### Customize Plot Titles
    
    1. Click the ``Customize Plot Titles`` button to reveal title customization options
    2. You can customize titles for all five visualization types:
        - ``Heatmap Plot`` (default: "Interaction Heatmap")
        - ``Prevalence Plot`` (default: "Interaction Prevalence Analysis")
        - ``Lifetime Plot`` (default: "Interaction Lifetimes Distribution")
        - ``Time Series Plot`` (default: "Interactions Over Time")
        - ``Network Plot`` (default: "Interaction Network")
    3. Enter your desired titles and click ``Apply``
    
    ### Transpose Axes
    
    - Click the ``Transpose Axes`` button to swap X and Y axes in the Heatmap view
    - The button highlights when transpose is active
    - Useful for alternate viewing perspectives or when one selection has many more elements than the other


## 5. Supported Visualizations

??? note "**Sele1 vs Sele2**"

    In this Tab, an overview of all interaction pairs between two selections is reported.

    ![Image title](assets/Tab_1.png){width="100%"}

    /// caption
    The ``Sele1 vs Sele2`` tab provides an overview of all interaction pairs between two selections.
    ///

    - X-axis``: Atoms/Residues of Selection 1 
 
    - ``Y-axis``: Atoms/Residues of Selection 2

    - ``Cells``
      - *Color:* Interaction type (following an interaction priority system — if multiple interactions occur, the most relevant is shown).
      - *Number:* Interaction prevalence (shown upon activation of the ``Show`` option adjacent to the ``Prevalence`` filter).
 
     **Hover**

      - *Sel1:* Atom/Residue of the first selection
      - *Note1:* Annotation corresponding to Sel 1
      - *Sel2:* Atom/Residue of the second selection
      - *Note2:* Annotation corresponding to Sel 2
      - *Interaction:* Type of interaction
      - *Prevalence:* Percentage of frames with this interaction

    - ``Transpose Button`` to swap X/Y axes for alternate views

??? note "**Interactions Prevalence**"

    In this Tab, the prevalence of interactions is visualized as a bar chart. In the top panel, the interactions of Selection 1 are shown, while the bottom panel focuses on Selection 2.
    
    ![Image title](assets/Tab_2.png){width="100%"}

    /// caption
    The ``Prevalence`` tab visualizes the prevalence of interactions as a bar chart, with Selection 1 in the top panel and Selection 2 in the bottom panel.
    ///


    - ``X-axis``: Components from Selection 1/Selection 2

    - ``Y-axis``: Prevalence (%)

    - ``Bars``:
        - *Color:* Interaction type
        - *Height:* Prevalence percentage

    - **Hover**:
      - Selection_1: Prevalence of interaction in Selection 1
      - Selection_2: Prevalence of interaction in Selection 2
      - Interaction: Type of interaction
      - Prevalence: Percentage of frames with this interaction
      - Annotation: Custom annotation for the interaction

??? note "**Interaction Lifetimes**"

    box plots representing the distribution of interaction lifetimes.

    ![Image title](assets/Tab_3.png){width="100%"}

    /// caption
    The Life Time tab presents box plots displaying the distribution of interaction lifetimes across the trajectory.
    ///

    - ``X-axis``: Interaction pairs (Selection 1 – Selection 2)

    - ``Y-axis``: Duration of interactions (in frames)

    - ``Color``: Interaction type

    - **Hover**:
      - Pair: Selection 1 – Selection 2
      - Interaction: Type of interaction
      - Prevalence: Percentage of frames with this interaction
      - Lifetime: Number of frames the interaction was present
      - Frame Range: Start and end frames where the interaction was observed

??? note "**Time Series**"

    Time series plot showing interaction presence over time.

    ![Image title](assets/Tab_4.png){width="100%"}

    /// caption
    The ``Time Series`` tab displays a time series plot showing the presence of interactions over time.
    ///

    1. ``Main Panel``
        - ``X-axis``: Frame number
        - ``Y-axis``: Interaction pairs
        - ``Dots``: Interaction presence

    2. ``Top Histogram``
        - Number of interactions per frame

    3. ``Side Histogram``
        - Prevalence of interactions across frames


    **Hover**:

    - Frame: Frame number
    - Selection Pair: Interacting residues
    - Interaction: Type of interaction
    - Prevalence: Percentage of frames with this interaction


??? note "**Network**"

    This tab visualizes the interaction network between two selections, where nodes represent components and edges represent interactions.

    ![Image title](assets/Tab_5.png){width="100%"}
    
    /// caption
    The ``Network`` tab visualizes the interaction network between two selections, where nodes represent components and edges represent interactions.
    ///

    ``Interpretation Guide:``

    - ``Nodes``:
        - 🔵 Blue = Selection 1 components
        - 🔴 Red = Selection 2 components
        - Size = Number of connections (degree)

    - ``Links``:
        - ``Color``: Interaction type
        - ``Width``: Interaction prevalence
        - ``Length``: Inversely related to prevalence (stronger interactions appear shorter)

     **Hover**:

        - Node information: Component name, number of connections
        - Edge information: Interaction type, prevalence percentage

    ### Network Controls Panel

    When you navigate to the ``Network`` tab, an additional control panel automatically appears on the left side of the screen. This panel provides advanced physics and visual customization options for fine-tuning the network layout.
    
    /// caption
    The Network Controls panel appears automatically when viewing the Network tab, offering fine-grained control over the network layout and appearance.
    ///

    #### Physics Controls
    
    Fine-tune the force-directed layout algorithm:
    
    - ``Gravitational Constant`` (-500 to 0, default: -200)
        - Controls the attractive force between nodes
        - More negative values = stronger attraction
        - Use lower values (closer to 0) for more spread-out networks
    
    - ``Central Gravity`` (0 to 0.1, default: 0.005)
        - Pulls all nodes toward the center of the visualization
        - Higher values = tighter, more compact networks
        - Lower values = more distributed layouts
    
    - ``Spring Length`` (50 to 500, default: 200)
        - Sets the natural resting distance between connected nodes
        - Larger values = more space between nodes
        - Smaller values = more compact arrangements
    
    - ``Spring Constant`` (0.1 to 2.0, default: 0.5)
        - Controls the stiffness of connections between nodes
        - Higher values = stiffer springs, more rigid structure
        - Lower values = more flexible, fluid layouts
    
    - ``Avoid Overlap`` (0 to 1, default: 0.8)
        - Prevents nodes from overlapping each other
        - Higher values = stronger repulsion, better separation
        - Set to 0 to allow overlapping (not recommended)

    #### Node Controls
    
    Customize node appearance:
    
    - ``Min Node Size`` (10 to 30, default: 20)
        - Minimum size for nodes with few connections
        - Ensures small nodes remain visible
    
    - ``Max Node Size`` (40 to 100, default: 50)
        - Maximum size for highly connected hub nodes
        - Prevents large nodes from dominating the visualization

    #### Edge Controls
    
    Customize edge (link) appearance:
    
    - ``Min Edge Width`` (1 to 10, default: 5)
        - Minimum width for edges with low prevalence
        - Ensures all connections remain visible
    
    - ``Max Edge Width`` (10 to 30, default: 15)
        - Maximum width for edges with high prevalence
        - Visual emphasis on strong interactions

    #### Simulation Controls
    
    Control the physics simulation:
    
    - ``Stabilization Iterations`` (100 to 2000, default: 1000)
        - Number of physics simulation steps before displaying the network
        - Higher values = more stable final layout (but slower)
        - Lower values = faster rendering (but may be less stable)
    
    - ``Physics Enabled`` (Toggle switch)
        - ``On``: Physics simulation runs continuously; nodes respond to forces
        - ``Off``: Disable physics for manual node positioning and static layout

    #### Apply Settings
    
    After adjusting any parameters, click the ``Apply Settings`` button to update the network visualization with your new configuration. The network will re-stabilize using the new physics parameters.

    !!! tip "Network Optimization Tips"
        - ``For large networks``: Reduce ``Stabilization Iterations`` to 200-500 for faster rendering
        - ``For dense networks``: Increase ``Spring Length`` and reduce ``Central Gravity`` to spread nodes apart
        - ``For publication figures``: Disable ``Physics Enabled`` after finding a good layout, then manually adjust node positions
        - ``For exploring dynamics``: Keep ``Physics Enabled`` and adjust forces in real-time to see immediate effects
        - ``For clearer visualization``: Increase ``Avoid Overlap`` to 1.0 to ensure perfect node separation


## 10. Footer Resources

The footer provides quick access to important resources:

- ``Docs & Tutorials``: Complete InterMap documentation and step-by-step tutorials
- ``InterMap's Paper``: Scientific publication with methodology and validation
- ``Delta-Research-Team/intermap``: Source code, issue tracker, and latest releases on GitHub
