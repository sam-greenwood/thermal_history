# Overview

## Leeds thermal evolution model
Open source Python model for calculating the 2D (radius and time) thermo-chemical evolution of planetary interiors e.g. Earth/Mars/Mercury.
`thermal_history` contains numerical implementations for a planet spearated into 3 regions: A silicate mantle, stably stratified liquid outer core, and an isentropic iron alloy core.
All regions are optional allowing flexible calculations based on the interests of the user. For each region, there are a choice of methods available based on different
published models, whilst the framework of the code makes it easy to add additional methods as they are developed for ultimate flexibility.

Full documentation at: https://sam-greenwood.github.io/thermal_history_docs/_build/html/index.html

## Getting Started

Requirements:
1. Works with Linux/MacOS/Windows.
2. conda package manager from either miniconda (see https://docs.conda.io/en/latest/miniconda.html) or a full anaconda installation.
---
### Installation

#### From the command line (recommended)
Simply clone this repo and install the anaconda environment by running from a unix terminal:

```
git clone https://github.com/sam-greenwood/thermal_history.git
cd thermal_history
conda env create -f environment.yml
conda activate thermal_history
```

#### Installing with the anaconda navigator GUI
Download the code by clicking the green code button, then either download the zip file (then extract the files) or download with the github desktop app.
Use the anaconda navigator app to install a new environment, importing the `environment.yml` file.

#### Importing thermal_history in your script
To be able to import `thermal_history` from anywhere, you can append to the PYTHONPATH system variable by putting `export PYTHONPATH=$PYTHONPATH:"full/path"` into your `~/.bashrc` file (or equivalent depending on your shell), where `full/path` is the path to this repo in your file system. Running `source ~/.bashrc` will then set this variable for the current shell. The command `python -c 'import sys; [print(x) for x in sys.path]; import thermal_history'` will print out all locations where python looks for imports and tries to import thermal_history. Make sure the conda environment that was installed has been activated: `conda activate thermal_history`.

Windows users may not be able to define the PYTHONPATH variable, in which case you can use conda to add the correct path. With the anaconda prompt program open, first install conda-build to the thermal_history environment with `conda activate thermal_history` then `conda install conda-build`. Next, naviate to the top level of this git repo and run `conda develop .`. This will add it to the python path whenever the thermal_history environment is active.

Note installing conda-build may want to downgrade your python version (e.g. from 3.10 to 3.7). You can specify `conda install conda-build python=3.10` for example to retain your later python version.


#### ARM on macOS
The latest Apple computers have switched to ARM processors. At the time of writing, the version of numpy and numba conda defaults to installing are incompatible due to the latest numba package not yet compiled on ARM. Follow instructions for installing the compatible version of numpy for numba after trying to first run the code. (`conda install numpy=X.X.X`)

### Documentation

Full documentation is available at https://sam-greenwood.github.io/thermal_history_docs/_build/html/index.html