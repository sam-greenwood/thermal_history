# Overview

Mantle model based on Knibbe and Van Westrenen (2018). Note that a small correction to the steady state conduction solution for the thermal profile in the lithosphere/crust/regolith has been made after communication with the authors.

# Expected attributes of other regions:

Core must have the following attributes (i.e. model.core.X):
| Attribute  |            Description                                           |
|------------|------------------------------------------------------------------|
|`T_cmb`     |  Temperature at the Core-Mantle Boundary                         |

If not using a core model, `T_cmb` must be set manually, e.g.:

```python
from thermal_history.model import Parameters, setup_model

#Load parameters and setup model
prm = Parameters('/path/to/parameters_file')
model = setup_model(prm, mantle_method='knibbe_westrenen18')

#Evolve model
dt = 1e6*prm.ys
for i in range(1000):
    model.core.T_cmb = 3000 #Set core temperature.
    model.evolve(dt)
    
```

## References:
- Knibbe, J.S. and van Westrenen, W., 2018. The thermal evolution of Mercury's Fe–Si core. Earth and Planetary Science Letters, 482, pp.147-159.

---