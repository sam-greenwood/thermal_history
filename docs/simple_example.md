### Running a simple example
For a basic example we will solve for the thermal evolution of an isothermal iron core with radius $r_c$, constant density ($\rho$), and specific heat capacity ($C_p$):

$Q = -\frac{4}{3}r_c^3\rho C_p \frac{\mathrm{d}T}{\mathrm{d}t}$

where $Q$ is the heat conducted away from the core by the mantle (although we are not including a mantle model so $Q$ will be imposed a priori).

For any calculation, we need a file with all of the necessary parameters. A parameters file for this example is included (`./examples/simple_test.py`).
Note that the parameters file is a python file. When we run the model, this file will be imported and anything defined inside will get assigned as a parameter.
Inside the file looks like:

```python
#Constants
ys = 60*60*24*365     #Seconds in a year
ev = 1.602e-19        #Electron volt
kb = 1.3806485e-23    #Boltzmanns constant
G  = 6.67e-11         #Gravitational Constant
Na = 6.022140857e23   #Avogadros Constant
Rg = 8.31446261815324 #Gas constant

core         = True
stable_layer = False
mantle       = False

Tc    = 5000    #Intial core temperature (K)
rho   = 6000    #Density (kg/m^3)
cp    = 800     #Specifc heat (J/K/kg)
rc    = 3480e3  #Core radius (m)
```

At the top are some physical constants that may be useful. These are also hard-coded into the parameters (see `thermal_history.model.Parameters`) but can sometimes be useful for converting variables into SI units here in the parameters file. The first 3 booleans tell the code which of the 3 regions we are interested in modelling, with the other variables defining necessary parameters for the method we'll be using. Since this is a very simple case we only have a few parameters, in practise you'll likley have a lot more! Also note that all the parameters are given in SI units. All of the methods included assume SI units for all variables to remove any confusion.

If you ever lose a parameters file and forget what needs to be defined within it, there is a handy utility function that can automatically generate one for a given set of desired methods:

```python
from thermal_history.utils import create_parameters_file

fname='my_new_parameters_file' #Name for parameters file to be created (.py will be appended automatically if not included)

#Change the keyword arguments from their default 'None' to the string name of the method you'd like for any region.
create_parameters_file(fname, core_method=None, stable_layer_method=None, mantle_method=None) 
```


First thing to do is load in the parameters, do this with

```python
from thermal_history.model import Parameters

prm = Parameters('./examples/simple_test.py')
```

Now with the parameters loaded into a class instance, we can setup the model and give it those parameters:

```python
from thermal_history.model import setup_model

model = setup_model(prm, core_method='simple_test')
```

Here we have told it that we want to use the 'simple_test' method for solving our core evolution. This will print some output to STDOUT as it does an internal check
to make sure all required parameters are contained in `prm` and prints the initial conditions (as defined by the 'simple_test' method).

Now we can iterate the model. The user is left to create their own loop to iterate the model, to allow any custom code (e.g. performing checks) at each timestep.
To advance the model we need to call the evolve() method:

```python
dt = 1e6*prm.ys   #1 Myr time step

for i in range(1000):
    model.evolve(dt)
```

This will produce an error! Remember, we're not solving for the mantle so the model has no idea what the CMB heat flow should be. We can set this manually:


```python
dt = 1e6*prm.ys   #1 Myr time step
model.mantle.Q_cmb = 1e12 #1TW

for i in range(1000):
    model.evolve(dt)
```

This now works! By defining our own loop, we could change the value of model.mantle.Q_cmb each iteration if we wish.

All of our results are saved into `model.save_dict`. This dictionary is structured in an equivalent manner to `model`, with a sub-dictionary for each region being modelled.
In our case we just have a core model, so the dictionary `model.save_dict.core` contains the attributes of `model.core` from each timestep contained within a list.

```python
import matplotlib.pyplot as plt

data = model.save_dict['core']

plt.plot(data['T'])
plt.show()
```

If numerical operations are needed to apply to the output, there is a quick built in fuction to convert all lists in `model.save_dict` into numpy arrays e.g.

```python
data = model.save_dict_to_numpy_array()

#Convert time to Myrs
time = data['core']['time']/(1e6*prm.ys)
```


If we wish to save our output data, either the user can write their own code to perform this action in a way of their preference or the built in `model.write_data(fname)` function can be used. This function writes the save_dict to a binary file using pickle. When saving, the parameters, time of save, and any dictionary called 'profiles' that is an attribute of any region (commonly used by various methods that ship with this code e.g. model.core.profiles) get added as well. Additionally, save_dict_to_numpy_array() is used just before saving.

```python
model.write_data('output.pik')

import pickle

loaded_data = pickle.load(open('output.pik','rb'))
```