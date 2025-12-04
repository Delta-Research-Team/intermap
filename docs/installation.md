Once you have git and conda/miniconda installed in your system, run:

```bash
git clone https://github.com/Delta-Research-Team/intermap.git
cd intermap
conda env create -f environment.yml
conda activate intermap
```     

The last command activates the newly created conda environment named `intermap` and you will need to run it each time you want to use InterMap.

## Updating
To ensure at any moment that you have the latest version of InterMap, navigate to the previously cloned repository directory and run:

```bash
git pull origin master
conda env update -f environment.yml
```

[//]: # ()
[//]: # (## 0. Prerequisites)

[//]: # ()
[//]: # (Before installing InterMap, ensure you have the following prerequisites installed on your system:)

[//]: # ()
[//]: # (- **Git**: Required for cloning the repository)

[//]: # (- **Conda or Miniconda**: Required for managing the Python environment)

[//]: # ()
[//]: # (## 1. Clone the Repository)

[//]: # ()
[//]: # (First, clone the InterMap repository from GitHub to your local machine and navigate to the cloned directory. Open your)

[//]: # (terminal and run the following commands:)

[//]: # ()
[//]: # (```bash)

[//]: # (git clone https://github.com/Delta-Research-Team/intermap.git)

[//]: # (cd intermap)

[//]: # (```)

[//]: # ()
[//]: # (## 2. Create the Conda Environment)

[//]: # ()
[//]: # (The InterMap repository includes an `environment.yml` file that contains all the necessary dependencies. Create a new)

[//]: # (conda environment using this file:)

[//]: # ()
[//]: # (```bash)

[//]: # (conda env create -f environment.yml)

[//]: # (```)

[//]: # ()
[//]: # (This command will:)

[//]: # ()
[//]: # (- Read the `environment.yml` file)

[//]: # (- Create a new conda environment with all required dependencies)

[//]: # (- Install the specified Python version and packages)

[//]: # ()
[//]: # (!!! note "Environment Creation Time")

[//]: # ()
[//]: # (    The environment creation process may take several minutes depending on your internet connection and system performance, as conda needs to download and install all the required packages.)

[//]: # ()
[//]: # (## 3. Activate the Environment)

[//]: # ()
[//]: # (Once the environment is created successfully, activate it:)

[//]: # ()
[//]: # (```bash)

[//]: # (conda activate intermap)

[//]: # (```)

[//]: # ()
[//]: # (You'll need to activate this environment every time you want to use InterMap. You can verify that the environment is)

[//]: # (active by checking that `&#40;intermap&#41;` appears at the beginning of your command prompt or just by running)

[//]: # (```conda info --envs```, which will list all your conda environments and highlight the active one.)

[//]: # ()
[//]: # (## 4. Keep InterMap Updated)

[//]: # ()
[//]: # (To update InterMap to the latest version:)

[//]: # ()
[//]: # (```bash)

[//]: # (# Navigate to the repository directory)

[//]: # (cd intermap)

[//]: # ()
[//]: # (# Pull the latest changes)

[//]: # (git pull origin master)

[//]: # ()
[//]: # (# Update the environment if needed)

[//]: # (conda env update -f environment.yml)

[//]: # (```)

[//]: # ()
[//]: # (!!! warning "Environment Updates")

[//]: # ()
[//]: # (    When updating, the `environment.yml` file may have changed. Always run `conda env update -f environment.yml` after pulling updates to ensure you have the latest dependencies.)
