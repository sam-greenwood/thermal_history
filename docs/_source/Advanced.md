# Adding new methods

Whilst this code ships with certain methods for the 3 regions (core/stable layer/mantle), the user can add alternative ones that they develop. To do so is easy enough but needs to be done in a specific way to interface with the whole python package.

## 1. Overview of the package

To start it makes sense to explain the overall structure of the python package, here is an abbrieviated overview:

```
thermal_history
      |
      |
   ___|____________________________________________________________....
  |                  |          |            |                |
model_classes.py   model.py  __init__.py   core_models   stable_layer_models
                                             |                            |
              _______________________________|                            |       
             |                                                            |
             |                                                            |
   __________|_____....                                                   |
  |           |                                                           |
leeds   simple_test                                                       |
  |           |_________________________________________                  |
  |                                                     |                 |    
  |_____________________________________             ___|________         |
  |          |              |           |           |           |         |
main.py   __init__.py   _readme.md   routines      main.py  __init__.py   |
                                        |                                 |
                                        |                                 |
                            ____________|___________....                  |
                           |                |                             |
                        profiles.py   chemistry.py                        |
                                                                          |
                                                                          |
                                                                          |
                                                            ______________|______...
                                                           |               |
                                                        _init__.py   leeds_thermal
                                                                           |
                                      _____________________________________|
                                     |          |              |           |        
                                  main.py   __init__.py   _readme.md   routines
                                                                         |
                                                                         |
                                                                      ___|___...   
                                                                     |
                                                                diffusion.py
```

The code for any specific numerical model is contained within one of the sub-packages depending on it's type (`core_models`/`stable_layer_models`/`mantle_models`). Let's focus on the `core_models`. Here there are 2 different methods for solving the core thermal evolution: `leeds` and `simple_test`. They are each contained within their own folder (treated as another sub-package due to the presence of `__init__.py`), meaning they are isolated from one another and removing any chance of getting them mixed up. When a model is initialised, for example, with `setup_model(prm, core_method='simple_test')`, the code will look to import `thermal_history.core_models.simple_test.main`. Other files can exist within the folder (e.g. `core_models.leeds.routines` contains functions that are called by code in `main.py`) but `main.py` is the required file to be there.

An additional file is present in the model folders: `_readme.md`. This will be explained in the documentation section below.



## 2. Required format for main.py

There are 3 functions expected to be defined in `main.py`. They are: `setup`, `evolve`, `update`. All take a single argument (the ThermalModel class that is created by `setup_model()` is passed to them) and none of them return anything as they instead set the attributes of ThermalModel. There are also 2 dictionaries, `required_params` and `optional_params`, that must de defined within `main.py`. Let's take a look at the `main.py` in simple_test, starting with the top of the file:


```python

#Description of this model3
Description = 'Simple example core model.'

#import our logger
import logging
logger = logging.getLogger(__name__)

import numpy as np
```

These are all optional but can help the code determine some useful information about the model. `Description` is just that, a brief overview of what this model is, used when the utility function `thermal_history.utils.. `compatibility` is a dictionary that tells the code if this method works when used in conjuction with other methods. As an example, this doesn't work with any of the methods for other regions, so can only be run on it's own. Next we have imported the logger, useful if you want to print some diagnostic information to the log file but not necessary. Numpy is imported at the end since it is needed by our model get the value of pi in `evolve()`.

Next we have the first of the required functions: `setup()`


```python
def setup(model):
    #Sets up any inital conditions

    core = model.core
    prm = model.parameters

    core.T = prm.T #Set initial core temperature
```

This function is called when `thermal_history.model.setup_model()` is run. It should be used to initialise any initial conditions and anything else needed before `evolve()` is called. We can see here that the parameters we load in, are accessible at `model.parameters`.

Next is the primary function, `evolve()`:

```python
def evolve(model):
    #Main model

    core = model.core
    prm = model.parameters

    Q = model.mantle.Q_cmb

    mass = (4/3)*np.pi*prm.rc**3 * prm.density

    dT_dt = Q/(-mass*prm.cp)

    #Variables to be saved need to attributes of model.core
    core.dT_dt = dT_dt
```

This contains the guts of the numerical scheme. Here we gather together necessary variables/parameters to calculate the rate of change of temperature of the core. Note that at the end, we make sure to set `dT_dt` as an attribute of `model.core`. This ensures that a) it is accessible by the `update()` function and b) gets saved each timstep into `model.save_dict`. Note that some attributes will never get added to `save_dict`, these are anything that isn't an int/float/numpy array and anything called profiles. Radial profiles can contain lots of data and can potentially slow down the model significantly. If the user wishes to save the radial profiles at each timestep, they can do so in their prefered way when calling `model.evolve()`.

Next we have the final required function, `update()`:

```python
def update(model):
    #Apply any changes between timesteps. Attributes of model.core will be saved to 

    core = model.core

    #Update core temperature
    core.T += core.dT_dt*model.dt

    #Example use of the logger
    if core.T < 0:
        logger.critical('Core temperature is negative!')
```

This is called after `evolve()` and after attributes have been appended to the save_dict. It is used to update the model to the next timestep. Here we also have used the logger to warn us of a nonsensical result just in case. Note that the timestep, `dt`, is an attribute of `model` not `model.core`.

Next is an optional function which is used to format what is printed to screen during the calculation.

```python
def progress(model):
    #Text to be printed to STDOUT on each iteration

    core = model.core

    text = f'    temperature = {core.T:.2f} K'

    return text
```

This should return a string and can be used to include some key information about how the calcualtion is going. By default just the iteration counter is printed by this can append some extra text alongside any of text from other regions if they are also included.

Finally, there are 2 dictionaries that tell the code what parameters are needed so that anything missing can be flagged during the check and is also used to automatically construct parameter files:

```python
#list required parameters {'name': 'description'}
required_params = {'rc': 'Radius of the core',
                    'density': 'Core density.',
                    'cp': 'Specific heat capacity',
                    'T': 'Initial temperature of the core'}

#list optional parameters that if not specified in model.parameters, will get set to a default value.
#{'name': ('description', default_value)}
optional_params = {'example_variable': ('An example optional variable for demonstration purposes', True)}
```

These could just be set to empty dictionaries which may be useful while developing a model but once a model is finalised (and particularly if someone else is to use it) it is a good idea to make sure these are fully populated.




## 3. Documentation
These instructions are aimed at anyone looking to have their changes to the thermal_history package reflected in the documentation.

Documentation for this repo is handled by another Github repo: https://github.com/sam-greenwood/thermal_history_docs. See the `.github/workflows` to see how the docs are automatically built and published to Github pages.

Whenever changes are made to the main thermal_history Github repo are made, this triggers the documentation repo to copy the latest version of the main repository, rebuild the documentation, and publish it to a Github pages site. Note the use of PAT tokens are required to give one repo permission to trigger a workflow once a workflow on another repo has run. Sphinx is used to build the documentation and its configuration file, templates, and any additional information such as the file containing these instructions can be found in `docs/`. A handy script `docs/create_docs.sh` is used to run the commands to build the docs and can be run locally and contains info on required packages to be installed to run. Sphinx automatically searches through the docstrings throughout the thermal_history package and builds the html files for the site. The `docs/_source/conf.py` file is setup to support the Numpy style docstrings for functions/classes.

Sphinx first generates rst files before building the html files from them. The documentation is setup to automatically generate rst files for each module that exists within thermal_history using the included templates. `docs/_source/index.rst` is not automatically generated and is used to control the appearance of the top level page of the documentation. Markdown files have been used to write the README and other pages due to their readability when compared to rst, with Sphinx supporting them through the use of the Myst parser (https://myst-parser.readthedocs.io/en/latest/index.html) as specified in the `docs/_source/conf.py` file.

As previously mentioned, each model folder (e.g. `thermal_history/core_models/leeds`) may contain an optional markdown file. This markdown file contains more in depth information about that particular model and is used when the documentation is built by having the `__init__.py` file include it into its docstring. This ensures that the markdown is placed onto the main page about the specific model in the API. For example, the docstring of the `__init__.py` file in the `leeds` core model directory reads:

```python
"""
.. include:: ../../thermal_history/core_models/leeds/_readme.md
   :parser: myst_parser.sphinx_
"""
#Path to readme is relative to the docs/_source directory.
#You may need to update the 'docutils' python package to get
#the myst_parser support
```

This tells the documentation tool sphinx to add the markdown to the docstring using the myst parser (which adds markdown support to sphinx). You can see the result of this here: https://sam-greenwood.github.io/thermal_history_docs/_build/html/thermal_history.core_models.leeds.html

## 4. Known Issues

Sub/Super adiabatic oscillation with a stratified layer in the core. As a stratified layer begins to grows, it can do so even after Q_c > Q_a due to the thermal diffusion timescale being in mismatch with the system timestep. The layer can then cool faster than the adiabat and it will be eroded. At this point the core assumes an adiabatic profile which is inconsistent, causing an increase in Qc which causes a layer to form. This process can recur indefinitely whilst the layer is thin and easily eroded. By setting mix_layer = True (False by default) any super-adiabatic sub region of the stable is mixed. The temperature profile will then maintain a more realistic form as the layer begins to mix from the CMB down towards rs, reducing the jump in Tc when the layer finally erodes away. This routine takes the spatial grid, the quantity to conserve (in this case the thermal energy) minus the factor of temperature, and the temperature relative to the adiabat. It then mixes the profile to ensure the gradient is equal to or greater than 0 everywhere whilst conserving the total heat. This oscillation often only occurs at high light element concentrations. This routine has not been thoroughly tested. Use at your own discretion.

Mantle melting. Whilst the logic of the mantle melting appears to be correct, there is a bug whereby a mantle melt volume > 0 is not possible without the fine scale temperature being super-solidus at some point initially. This issue is under investigation.
