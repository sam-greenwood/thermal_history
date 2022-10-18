# Adding new methods

Whilst this code ships with certain methods for the 3 regions (core/stable layer/mantle), the user can add alternative ones that they develop. To do so is easy enough but needs to be done in a specific way to interface with the whole python package.

## 1. Overview of the package

To start it makes sense to explain the overall structure of the python package, here is an abbrieviated overview:

```
thermal_history
|
|   model_classes.py
|   model.py
|   __init__.py
|___core_models
|   |
|   |   __init__.py
|   |___leeds
|   |   |
|   |   |   main.py
|   |   |   __init__.py
|   |   |___routines
|   |       |
|   |       |   profiles.py
|   |       |   chemistry.py
|   |       ....
|   |
|   |___simple_test
|       |
|       |   main.py
|       |   __init__.py
|          
|
|
|___stable_layer_models
|   |
|   |   __init__.py
|   |___leeds_thermal
|       |
|       |   main.py
|       |   __init__.py
|       |___routines
|           |
|           |   diffusion.py
|           ...
...

```

The code for any specific numerical model is contained within one of the sub-packages depending on it's type (`core_models`/`stable_layer_models`/`mantle_models`). Let's focus on the `core_models`. Here there are 2 different methods for solving the core thermal evolution: `leeds` and `simple_test`. They are each contained within their own folder (treated as another sub-package due to the presence of `__init__.py`), meaning they are isolated from one another and removing any chance of getting them mixed up. When a model is initialised, for example, with `setup_model(prm, core_method='simple_test')`, the code will look to import `thermal_history.core_models.simple_test.main`. Other files can exist within the folder (e.g. `core_models.leeds.routines` contains functions that are called by code in `main.py`) but `main.py` must be there.


## 2. Required format for main.py

There are 3 functions expected to be defined in `main.py`. They are: `setup`, `evolve`, `update`. All take a single argument (the ThermalModel class that is created by `setup_model()` is passed to them) and none of them return anything as they instead set the attributes of ThermalModel. There are also 2 dictionaries, `required_params` and `optional_params`, that must de defined within. Let's take a look at the `main.py` in simple_test, starting with the top of the file:


```python

#Description of this model3
Description = 'Simple example core model.'

#List individually confirmed compatibility with methods for other regions (give the names of the methods in the list if any).
compatibility = {'stable_layer': [],
                 'mantle': []}

#import our logger
import logging
logger = logging.getLogger(__name__)

import numpy as np
```

These are all optional but can help the code determine some useful information about the model. `Description` is just that, a brief overview of what this model is. `compatibility` is a dictionary that tells the code if this method works when used in conjuction with other methods. As an example, this doesn't work with any of the methods for other regions, so can only be run on it's own. Next we have imported the logger, useful if you want to print some diagnostic information to the log file but not necessary. Numpy is imported at the end since it is needed by our model get the value of pi in `evolve()`.

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