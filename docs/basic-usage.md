## 1. Preparing the configuration file

To run InterMap, you must prepare a configuration file specifying the parameters for your analysis. The
configuration file is a regular text file that contains various sections, each corresponding to different aspects of the setup.

In the following example (that can be found under the `example` folder of the GitHub repository), we will go through each entry, explaining its purpose and how to set it up.

!!! example

    ```toml
      [generals]
      n_procs = 8
      job_name = prot-dna
      output_dir = ./outputs
      n_samples = 10
      n_factor = 1.5
      
      [topo-traj]
      topology = ./hmr_cionize_ions_solvent_sc150mM_WRAPPED_RENUM.pdb
      trajectory = ./traj_sample.dcd
      start = 0
      last = -1
      stride = 1
      chunk_size = 100
      
      [interactions]
      selection_1 = protein
      selection_2 = nucleic
      min_prevalence = 0
      resolution = atom
      interactions = all
      annotations = False
      
      [cutoffs]
      dist_cut_CloseContact = 3.0
      dist_cut_Hydrophobic = 4.5
      dist_cut_Ionic = 4.5
      dist_cut_Metalic = 2.8
      
      dist_cut_HA = 3.5
      dist_cut_DA = 3.5
      min_ang_DHA = 130
      max_ang_DHA = 180
      
      dist_cut_XA = 3.5
      dist_cut_XD = 3.5
      min_ang_DXA = 0
      max_ang_DXA = 90
      
      dist_cut_FaceToFace = 5.5
      min_ang_nn_FaceToFace = 0
      max_ang_nn_FaceToFace = 35
      min_ang_nc_FaceToFace =  0
      max_ang_nc_FaceToFace = 30
      
      dist_cut_EdgeToFace = 6.5
      min_ang_nn_EdgeToFace = 50
      max_ang_nn_EdgeToFace = 90
      min_ang_nc_EdgeToFace = 0
      max_ang_nc_EdgeToFace = 30
      intersect_radius_EdgeToFace = 1.5
      
      dist_cut_PiStacking = 5.5
      min_ang_nn_PiStacking = 0
      max_ang_nn_PiStacking = 90
      min_ang_nc_PiStacking = 0
      max_ang_nc_PiStacking = 30
      
      dist_cut_PiCation = 4.5
      min_ang_PiCation = 0
      max_ang_PiCation = 30
      dist_cut_PiAnion = 4.5
      min_ang_PiAnion = 0
      max_ang_PiAnion = 30
    ```

    *Parameters in the configuration file are internally validated, so if you set a parameter that is not recognized or has an invalid value, InterMap will raise an error and inform you about the issue.*

    *The only optional section is `[cutoffs]`, which allows you to customize the cutoffs for each interaction type. If you do not include this section, InterMap will use the default cutoffs shown here.*

### **`[generals]`**

??? note "n_procs"

    **Value:** int >= 1

    **Description:** The number of processes to use for parallel execution.

    **Tip:** Increasing this value can speed up the analysis, especially for large trajectories. However, the performance may not improve significantly beyond a certain point due to the overhead of managing multiple processes. 

??? note "job_name"

    **Value:** str

    **Description:** A name for the job. It will be used as the basename to create output files.

    **Tip:** Choose a descriptive name that reflects the system or analysis you are performing. This will help you identify the output files later.

??? note "output_dir"

    **Value:** valid path

    **Description:** The directory where the output files will be saved.

    **Tip:** Ensure the specified directory is in a writable location. If it does not exist, InterMap will attempt to create it.

??? note "n_samples"

    **Value:** int > 0

    **Description:** The number of samples in your trajectory that InterMap explores before starting the detection stage. This helps infer the maximum number of interactions in your system.

    **Tip:** A good starting point is to set this value to 10. If you encounter memory issues, you should increase it.

??? note "n_factor"

    **Value:** float

    **Description:** A factor to scale the maximum number of interactions detected among the `n_samples` frames. This is used to estimate the memory requirements for the analysis.

    **Tip:** A good starting point is to set this value to  1.5, which provides a reasonable estimate of memory requirements. If you encounter memory issues, you can midly increase it.

!!! info "Memory Allocation"

    The `n_samples` and `n_factor` parameters are used to estimate the memory requirements for the analysis. If you set them too low, InterMap may not be able to handle the whole trajectory, leading to memory errors. If unsure about the values, start with the default (10 and 1.5, respectively) and increase them later upon encountering issues.

### **`[topo-traj]`**

??? note "topology"

    **Value:** valid path

    **Description:** The path to the topology file (e.g., PDB or PSF file) that describes the molecular system.

    **Tip:** Ensure that the topology file is compatible with the trajectory file you will be using. InterMap suports the topology formats supported by MDAnalysis.

??? note "trajectory"

    **Value:** valid path(s)

    **Description:** The path to the trajectory file (e.g., DCD or XTC file) that contains the molecular dynamics simulation data. You can specify multiple trajectories by separating them with commas if they correspond to the same topology.

    **Tip:** Ensure the trajectory file is compatible with the topology file. If you have multiple trajectories, they should all be in the same format and correspond to the same molecular system. InterMap suports the trajectory formats supported by MDAnalysis

??? note "start"

    **Value:** int >= 0

    **Description:** The starting frame of the trajectory to analyze. If set to 0, it will start from the first frame.

    **Tip:** If you want to analyze a specific trajectory portion, you can set this value to the desired starting frame. 

??? note "last"

    **Value:** int >= -1

    **Description:** The last frame of the trajectory to analyze. If set to -1, it will analyze until the last frame of the trajectory.

    **Tip:** If you want to analyze a specific trajectory portion, you can set this value to the desired last frame. 

??? note "stride"

    **Value:** int >= 1

    **Description:** The stride for reading frames from the trajectory. If set to 1, it will read every frame; if set to 2, it will read every second frame, and so on.

    **Tip:** A larger stride can reduce the number of frames processed, which can speed up the analysis but may also miss some interactions. 

??? note "chunk_size"

    **Value:** int >= 1

    **Description:** The number of frames to read from the trajectory. The greater the chunk size, the more memory will be used, which can speed up the analysis. If you encounter memory issues, you can reduce this value.

    **Tip:** A good starting point is to set this value to the number of processors available. If you encounter memory issues, you can reduce it to a lower value. If you have enough memory, you can increase it to speed up the analysis.

!!! info "Supported Molecular Formats"

    The formats of the `topology` and `trajectory` files should be compatible with those supported by MDAnalysis, the underlying library InterMap uses for trajectory processing. The exhaustive list of supported formats can be found in the [MDAnalysis documentation](https://userguide.mdanalysis.org/stable/formats/index.html).

### **`[interactions]`**

??? note "selection_1"

    **Value:** str

    **Description:** The first atomic selection of atoms or residues to analyze. This can be any valid MDAnalysis selection string.

    **Tip:** You can use selections like `protein`, `resid 10-20`, or `name CA` to specify the atoms or residues you want to analyze. If you want to analyze interactions within a single selection, you can set both `selection_1` and `selection_2` to the same value.

??? note "selection_2"

    **Value:** str

    **Description:** The second atomic selection of atoms or residues to analyze. This can also be any valid MDAnalysis selection string.

    **Tip:** You can use selections like `protein`, `resid 10-20`, or `name CA` to specify the atoms or residues you want to analyze. If you want to analyze interactions within a single selection, you can set both `selection_1` and `selection_2` to the same value.

??? note "min_prevalence"

    **Value:** float >= 0

    **Description:** The minimum prevalence of interactions to consider. This threshold filters out interactions
    that occur less frequently than the specified value. If set to 0, all interactions will be considered.

    **Tip:** A good starting point is setting this value to 0, including all interactions. If you want to filter out rare interactions, you can increase this value to a percentage (e.g., 10 for 10% prevalence). The valid range is from 1 to 100.

??? note "resolution"

    **Value:** atom OR residue

    **Description:** The resolution of the interaction fingerprints. The available options are `atom` and `residue`. The
    `atom` resolution provides detailed information at the atomic level, while the `residue` resolution summarizes
    interactions at the residue level.

    **Tip:** If you want detailed information about interactions at the atomic level, set this value to `atom`. If you want a more summarized view at the residue level, set it to `residue`. The choice depends on the level of detail you need for your analysis. Note that the `atom` resolution will produce larger interaction fingerprints, which may require more memory and computation time.

??? note "interactions"

    **Value:** str OR list of interaction OR all

    **Description:** The types of interactions to compute. You can specify a single interaction type or a list of
    interaction types separated by commas. You can set this value to `all` to compute all supported interactions.

    The available interaction types are: `CloseContact`, `VdWContact`, `Hydrophobic`, `Anionic`, `Cationic`, `MetalDonor`, `MetalAcceptor`, `HBAcceptor`, `HBDonor`, `XBAcceptor`, `XBDonor`, `FaceToFace`, `EdgeToFace`, `PiStacking`, `PiCation`, `CationPi`, `PiAnion`, and `AnionPi`.

    **Tip:** To compute all interactions, set this value to `all`. If you are interested in specific interactions, specify them as a comma-separated list (e.g., `CloseContact, Hydrophobic, Cationic`). The choice of interactions depends on the specific analysis you want to perform. The declaration of the interaction types is case-sensitive, so use the correct capitalization.

??? note "annotations"

    **Value:** valid path OR False

    **Description:** The path to a file containing annotations for the selections. This file is prepared by specifying the corresponding selection string for each arbitrary label representing a particular part of your system. InterMap tracks these annotations and can be used later to group interactions by labels, which is helpful for
    analyzing interactions between different parts of a complex system that are not defined by conventional chains,
    segments, etc.

    **Tip:** If you want to use annotations, create a text file with the format shown below and specify its path here.
    You can set this parameter to False if you do not have/need annotations. 

    !!! example
    ```
    core-A-U1 = (resnum 1:67)  or (resnum 89:108) and segid PROA
    hairpin-A-U1 = (resnum 68:88) and segid PROA
    tail-A-U1 = (resnum -24:0) and segid PROA
    ```

!!! info "Selection Syntax"

    The `selection_1` and `selection_2` parameters (as well as the declared `annotations`) can be any valid [MDAnalysis selection string](https://userguide.mdanalysis.org/stable/selections.html). If both selections are equal, InterMap will compute the intra-molecular interactions.

### **`[cutoffs]`**

!!! info "Interactions Cutoffs"

    The `[cutoffs]` section is optional and allows you to define the cutoffs for each interaction type. InterMap will use the default cutoffs defined in the code if this section is not included.

    The cutoffs are specified in Angstroms (for distances) or Degrees (for angles) and are used to determine whether an interaction occurs between a group of atoms. The currently available interactions and cutoffs are listed in the [configuration file example above](basic-usage.md).

??? note "dist_cut_CloseContact"

    **Value:** float

    **Description:** The distance cutoff for the Close Contact interaction type. This is the maximum distance between two atoms to consider them in close contact.

    **Tip:** A good starting point is setting this value to 3.0 Angstroms, a standard threshold for close molecular-system contact.

??? note "dist_cut_Hydrophobic"

    **Value:** float

    **Description:** The distance cutoff for the Hydrophobic interaction type. This is the maximum distance between two atoms to consider them hydrophobic.

    **Tip:** A good starting point is setting this value to 4.5 Angstroms, a standard threshold for hydrophobic interactions in molecular systems.

??? note "dist_cut_Ionic"

    **Value:** float

    **Description:** The distance cutoff for the Cationic / Anionic interaction type. This is the maximum distance between two atoms to consider them under a salt bridge.

    **Tip:** A good starting point is setting this value to 4.5 Angstroms, a standard threshold for ionic interactions in molecular systems.

??? note "dist_cut_Metalic"

    **Value:** float 

    **Description:** The distance cutoff for the Metallic interaction type. This is the maximum distance between two atoms (one of them a metal) to consider them interacting.

    **Tip:** A good starting point is setting this value to 2.8 Angstroms, a standard threshold for metallic interactions in molecular systems.

??? note "dist_cut_HA"

    **Value:** float 

    **Description:** The distance cutoff for the maximum distance between a hydrogen atom and an acceptor atom to consider them in a hydrogen bond.

    **Tip:** A good starting point is setting this value to 3.5 Angstroms, a standard threshold for hydrogen-acceptor distances in hydrogen bonding.

??? note "dist_cut_DA"

    **Value:** float

    **Description:** The distance cutoff for the maximum distance between a donor atom and an acceptor atom to consider them in a hydrogen bond.

    **Tip:** A good starting point is setting this value to 3.5 Angstroms, a standard threshold for donor-acceptor distances in hydrogen bonding.

??? note "min_ang_DHA"

    **Value:** float 

    **Description:** The minimum angle in degrees for the Donor-Hydrogen-Acceptor (DHA) geometry in hydrogen bonding. This angle determines the directionality of the hydrogen bond.

    **Tip:** A good starting point is to set this value to 0 degrees. This allows for flexible hydrogen bond geometries, but you can increase it to restrict the geometry to more linear hydrogen bonds.

??? note "max_ang_DHA"

    **Value:** float 

    **Description:** The maximum angle in degrees for the Donor-Hydrogen-Acceptor (DHA) geometry in hydrogen bonding. This angle determines the directionality of the hydrogen bond.

    **Tip:** A good starting point is to set this value to 180 degrees, allowing for perfectly linear hydrogen bonds.

??? note "dist_cut_XA"

    **Value:** float

    **Description:** The distance cutoff for the X-Acceptor (XA) interaction type in halogen bonding. This is the maximum distance between a halogen atom and an acceptor atom.

    **Tip:** A good starting point is setting this value to 3.5 Angstroms, a standard threshold for halogen-acceptor distances in halogen bonding.

??? note "dist_cut_XD"

    **Value:** float

    **Description:** The distance cutoff for the X-Donor (XD) interaction type in halogen bonding. This is the maximum distance between a halogen atom and a donor atom.

    **Tip:** A good starting point is setting this value to 3.5 Angstroms, a standard threshold for halogen-donor distances in halogen bonding.

??? note "min_ang_DXA"

    **Value:** float

    **Description:** The minimum angle in degrees for the Donor-X-Acceptor (DXA) geometry in halogen bonding. This angle determines the directionality of the halogen bond.

    **Tip:** A good starting point is to set this value to 0 degrees, allowing for flexible halogen bond geometries.

??? note "max_ang_DXA"

    **Value:** float 

    **Description:** The maximum angle in degrees for the Donor-X-Acceptor (DXA) geometry in halogen bonding. This angle determines the directionality of the halogen bond.

    **Tip:** A good starting point is to set this value to 180 degrees, which restricts halogen bonds to more linear geometries.

??? note "dist_cut_FaceToFace"

    **Value:** float 

    **Description:** The distance cutoff for the Face-to-Face π-π interaction type. This is the maximum distance between the centroids of two aromatic rings to consider them in a face-to-face π-π stacking interaction.

    **Tip:** A good starting point is setting this value to 5.5 Angstroms, a standard threshold for face-to-face π-π stacking interactions.

??? note "min_ang_nn_FaceToFace"

    **Value:** float 

    **Description:** The minimum angle in degrees between the normal vectors of two aromatic rings for Face-to-Face π-π interactions. This ensures the rings are oriented adequately for stacking.

    **Tip:** A good starting point is to set this value to 0 degrees, allowing for parallel ring orientations.

??? note "max_ang_nn_FaceToFace"

    **Value:** float 

    **Description:** The maximum angle in degrees between the normal vectors of two aromatic rings for Face-to-Face π-π interactions. This ensures the rings are oriented adequately for stacking.

    **Tip:** A good starting point is to set this value to 35 degrees, which allows for slightly tilted but still effective π-π stacking.

??? note "min_ang_nc_FaceToFace"

    **Value:** float 

    **Description:** The minimum angle in degrees between the normal vector of one ring and the centroid vector for Face-to-Face π-π interactions. This parameter helps define the geometry of the stacking interaction.

    **Tip:** A good starting point is setting this value to 0 degrees, allowing optimal stacking geometries.

??? note "max_ang_nc_FaceToFace"

    **Value:** float 
    **Description:** The maximum angle in degrees between the normal vector of one ring and the centroid vector for Face-to-Face π-π interactions. This parameter helps define the geometry of the stacking interaction.

    **Tip:** A good starting point is setting this value to 30 degrees, which maintains good stacking geometry while allowing flexibility.

??? note "dist_cut_EdgeToFace"

    **Value:** float 

    **Description:** The distance cutoff for the Edge-to-Face π-π interaction type. This is the maximum distance between the centroids of two aromatic rings to consider them in an edge-to-face π-π interaction.

    **Tip:** A good starting point is setting this value to 6.5 Angstroms, a standard threshold for edge-to-face π-π interactions.

??? note "min_ang_nn_EdgeToFace"

    **Value:** float 

    **Description:** The minimum angle in degrees between the normal vectors of two aromatic rings for Edge-to-Face π-π interactions. This ensures the rings are oriented adequately for T-shaped interactions.

    **Tip:** A good starting point is to set this value to 50 degrees, which ensures a proper T-shaped geometry.

??? note "max_ang_nn_EdgeToFace"

    **Value:** float 

    **Description:** The maximum angle in degrees between the normal vectors of two aromatic rings for Edge-to-Face π-π interactions. This ensures the rings are oriented adequately for T-shaped interactions.

    **Tip:** A good starting point is to set this value to 90 degrees, allowing for perpendicular ring orientations typical of edge-to-face interactions.

??? note "min_ang_nc_EdgeToFace"

    **Value:** float 

    **Description:** The minimum angle in degrees between the normal vector of one ring and the centroid vector for Edge-to-Face π-π interactions. This parameter helps define the geometry of the T-shaped interaction.

    **Tip:** A good starting point is setting this value to 0 degrees, allowing optimal edge-to-face geometries.

??? note "max_ang_nc_EdgeToFace"

    **Value:** float 

    **Description:** The maximum angle in degrees between the normal vector of one ring and the centroid vector for Edge-to-Face π-π interactions. This parameter helps define the geometry of the T-shaped interaction.

    **Tip:** A good starting point is setting this value to 30 degrees, which maintains good T-shaped geometry while allowing flexibility.

??? note "intersect_radius_EdgeToFace"

    **Value:** float

    **Description:** The intersection radius for Edge-to-Face π-π interactions. This parameter defines the radius around the ring centroid used to determine if the rings are properly positioned for an edge-to-face interaction.

    **Tip:** A good starting point is to set this value to 1.5 Angstroms, which provides a reasonable intersection zone for edge-to-face interactions.

??? note "dist_cut_PiStacking"

    **Value:** float

    **Description:** The distance cutoff for general π-π stacking interactions. This is the maximum distance between the centroids of two aromatic rings to consider them in a π-π stacking interaction.

    **Tip:** A good starting point is setting this value to 5.5 Angstroms, a standard threshold for π-π stacking interactions.

??? note "min_ang_nn_PiStacking"

    **Value:** float 

    **Description:** The minimum angle in degrees between the normal vectors of two aromatic rings for general π-π stacking interactions. This ensures the rings are oriented adequately for stacking.

    **Tip:** A good starting point is to set this value to 0 degrees, allowing for parallel ring orientations.

??? note "max_ang_nn_PiStacking"

    **Value:** float 

    **Description:** The maximum angle in degrees between the normal vectors of two aromatic rings for general π-π stacking interactions. This ensures the rings are oriented adequately for stacking.

    **Tip:** A good starting point is to set this value to 90 degrees, allowing for both parallel and perpendicular ring orientations.

??? note "min_ang_nc_PiStacking"

    **Value:** float 

    **Description:** The minimum angle in degrees between the normal vector of one ring and the centroid vector for general π-π stacking interactions. This parameter helps define the geometry of the stacking interaction.

    **Tip:** A good starting point is setting this value to 0 degrees, allowing optimal stacking geometries.

??? note "max_ang_nc_PiStacking"

    **Value:** float 

    **Description:** The maximum angle in degrees between the normal vector of one ring and the centroid vector for general π-π stacking interactions. This parameter helps define the geometry of the stacking interaction.

    **Tip:** A good starting point is setting this value to 30 degrees, which maintains good stacking geometry while allowing flexibility.

??? note "dist_cut_PiCation"

    **Value:** float

    **Description:** The distance cutoff for π-cation interactions. This is the maximum distance between the centroid of an aromatic ring and a cation to consider them in a π-cation interaction.

    **Tip:** A good starting point is setting this value to 4.5 Angstroms, a standard threshold for π-cation interactions.

??? note "min_ang_PiCation"

    **Value:** float 

    **Description:** The minimum angle in degrees for π-cation interactions. This angle is measured between the normal vector of the aromatic ring and the vector from the ring centroid to the cation.

    **Tip:** A good starting point is setting this value to 0 degrees, allowing optimal π-cation geometries.

??? note "max_ang_PiCation"

    **Value:** float 

    **Description:** The maximum angle in degrees for π-cation interactions. This angle is measured between the normal vector of the aromatic ring and the vector from the ring centroid to the cation.

    **Tip:** A good starting point is setting this value to 30 degrees, which maintains good π-cation geometry while allowing flexibility.

??? note "dist_cut_PiAnion"

    **Value:** float

    **Description:** The distance cutoff for π-anion interactions. This is the maximum distance between the centroid of an aromatic ring and an anion to consider them in a π-anion interaction.

    **Tip:** A good starting point is setting this value to 4.5 Angstroms, a standard threshold for π-anion interactions.

??? note "min_ang_PiAnion"

    **Value:** float 

    **Description:** The minimum angle in degrees for π-anion interactions. This angle is measured between the normal vector of the aromatic ring and the vector from the ring centroid to the anion.

    **Tip:** A good starting point is to set this value to 0 degrees, allowing for optimal π-anion geometries.

??? note "max_ang_PiAnion"

    **Value:** float 

    **Description:** The maximum angle in degrees for π-anion interactions. This angle is measured between the normal vector of the aromatic ring and the vector from the ring centroid to the anion.

    **Tip:** A good starting point is setting this value to 30 degrees, which maintains good π-anion geometry while allowing flexibility.

## 2. Running InterMap

InterMap can be run in two ways: through the command line or as a Python module. The command line interface is the most
straightforward way to run InterMap, while the Python module allows for integration into larger workflows without
disrupting the automation process.

### 2.1 Command Line Interface

To run InterMap from the command line, you need to have a conda environment set up with InterMap installed. If you have
not done this yet, please refer to the [installation guide](installation.md) for instructions on installing InterMap
using conda.

Once you have activated the InterMap's conda environment, execute the following command in your terminal.

```bash
intermap <path_to_your_config_file>
```

### 2.2 Python Module

To run InterMap as a Python module, you can use the following code snippet in your Python scripts. This code can be found under the `example` folder of the GitHub repository (`example_script.py`).

```python
import bitarray.util as bu
from intermap import runner

# ----| 1. Execute InterMap with the given configuration file
cfg_path = './prot-dna_InterMap.cfg'
bit_dict = runner.execute(cfg_path=cfg_path)
first_element = bit_dict[list(bit_dict.keys())[0]]

# ----| 2. Decode the bitarrays if they are in bytes format
if isinstance(first_element, bytes):
    bit_dict = {k: bu.sc_decode(v) for k, v in bit_dict.items()}
else:
    bit_dict = bit_dict

# ----| 3. Continue with further processing of bit_dict as needed
print(f'Processed {len(bit_dict)} interactions from configuration: {cfg_path}')

```

## 3. Output Files

Upon successful execution, InterMap will generate the following output files in the specified `output_dir`.

- `{job_name}_InterMap.cfg`: A copy of the configuration file used for the analysis, which can be helpful for reference
  or reproducibility.

- `topology.format`: A copy of the topology file used in the job.

- `{job_name}_InterMap.log`: A log file containing information about the execution of InterMap, including any warnings
  Alternatively, errors encountered during the analysis.

- `{job_name}_InterMap.pickle`: A pickle file containing the analysis results. This file can be loaded later for
  further analysis or visualization with [InterVis](intervis.md). Also, you can use the [
  `interconvert`](#4-interconvert) command to convert this file into a CSV file.


## 4. InterConvert

To convert the output `{job_name}_InterMap.pickle` into a CSV file that can be easily read and analyzed, you can use the `interconvert` command: 

```bash
interconvert <path_to_intermap_pickle_output_file>  <path_to_your_new_csv_file>
```
