## Option 1: Via Pip (Recommended)

If you already have a working Python environment and want to use InterMap alongside your other tools, simply install it
via pip:

```bash
pip install intermap
```  

## Option 2: Via Conda (For Developers & Isolated Environments)

If you experience dependency conflicts, or if you plan to contribute to the code and run tests, we recommend creating an
isolated environment using our provided configuration:

```bash
# Clone the repository
git clone https://github.com/Delta-Research-Team/intermap.git
cd intermap

# Create and activate the environment
conda env create -f environment.yml
conda activate intermap

# Install the package in editable mode 
pip install -e .

Note: You will need to run conda activate intermap each time you open a new terminal to use the software.
```

## Running Tests

If you installed the package from source using Option 2, you can run the test suite from the root dir of the project.

```bash
NUMBA_DISABLE_JIT=1 pytest tests/
```  

!!! note "Update file paths before running tests"

Before running the tests, make sure to update the paths to the topology and trajectory files in the configuration file within the examples directory. Otherwise, this may lead to errors during test execution.


!!! warning "Disable Numba's JIT compiler for running the tests"

    You must use the previous command to explicitly disable Numba's JIT compiler when running the tests. Otherwise, tests may fail or stall due to compilation overhead.


## Updating

To ensure you always have the latest version of InterMap:

### If installed via Pip:

```bash
pip install --upgrade intermap

```  

### If installed via Conda (Source):

Navigate to the cloned repository directory and run:

```bash
git pull origin master
conda env update -f environment.yml

```  
