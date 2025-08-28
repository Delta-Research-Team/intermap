**InterMap** can be easily installed using conda. Below are the steps to install on your system using conda.

## Installation Steps

1. **Install conda**: If you haven't already, install [Anaconda](https://www.anaconda.com/products/distribution)
   or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) on your system.
2. **Create a new conda environment**: It's a good practice to create a new environment for each project. You can do
   this by running:
   ```bash
   conda create -n intermap python=3.8
   ```
   Replace `intermap` with your desired environment name.
3. **Activate the environment**: Activate the newly created environment:
   ```bash
   conda activate intermap
   ```
4. **Install InterMap**: You can install InterMap directly from the conda channel:
   ```bash
    conda install -c conda-forge intermap
    ```
5. **Verify the installation**: After installation, you can verify that InterMap is installed correctly by running:
6. ```bash
   intermap --version
   ```
   This should display the version of InterMap installed.
7. **Deactivate the environment**: When you're done using InterMap, you can deactivate the environment by running:
   ```bash
   conda deactivate
   ```
