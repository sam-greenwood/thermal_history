# Overview

## Leeds thermal evolution model
Open source Python model for calculating the 2D (radius and time) thermo-chemical evolution of planetary interiors e.g. Earth/Mars/Mercury.
`thermal_history` contains numerical implementations for a planet spearated into 3 regions: A silicate mantle, stably stratified liquid outer core, and an isentropic iron alloy core.
All regions are optional allowing flexible calculations based on the interests of the user. For each region, there are a choice of methods available based on different
published models, whilst the framework of the code makes it easy to add additional methods as they are developed for ultimate flexibility.

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
To be able to import `thermal_history` from anywhere, you can append to the PYTHONPATH system variable by putting `export PYTHONPATH=$PYTHONPATH:"full/path/to/directory/containing_thermal_history"` into your `~/.bashrc` file (or `~/.zshrc` or equivalent depending on your shell) and running `source ~/.bashrc`. The active terminal and new terminals will allow python to always be aware of that location for imports (proper setup.py coming sometime soon!).

Otherwise you'l need to add to the system path with
```python
import sys
sys.path.append('/full/path/to/directory/containg_thermal_history')
```
before thermal_history is imported. Windows users beware, you'll need to escape each backslash and watch out for special unicode characters if copy and pasting from the file browser, you may need to type it out manually.

#### ARM on macOS
The latest Apple computers have switched to ARM processors. At the time of writing, the version of numpy and numba conda defaults to installing are incompatible due to the latest numba package not yet compiled on ARM. Follow instructions for installing the compatible version of numpy for numba after trying to first run the code. (`conda install numpy=X.X.X`)

### Documentation

Full documentation is included in the `docs/` folder, with a pre-built html page at `docs/_build/html/README.html`. Open this file in a web browser to see them (they will be hosted on readthedocs soon!).