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
