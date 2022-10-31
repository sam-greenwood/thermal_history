# Overview:

Simplified implimentation of the mantle model of Thiriet et al. (2019). Stagnant lid thickness is treated as constant in time and is prescribed rather than solved for. Used in the Mars study of Greenwood et al. (2021)

# Expected attributes of other regions:

Core must have the following attributes (i.e. model.core.X):
| Attribute   |            Description                                           |
|-------------|------------------------------------------------------------------|
|`T_cmb`      |  The core temperature at the CMB                                 |

If not using a core model, `T_cmb` must be set manually, e.g.:

```python
from thermal_history.model import Parameters, setup_model

#Load parameters and setup model
prm = Parameters('/path/to/parameters_file')
model = setup_model(prm, mantle_method='driscol_bercovici14')

#Evolve model
dt = 1e6*prm.ys
for i in range(1000):
    model.core.T_cmb = 3000 #Set core temperature.
    model.evolve(dt)
    
```

## References:
- Thiriet, M., Breuer, D., Michaut, C. and Plesa, A.C., 2019. Scaling laws of convection for cooling planets in a stagnant lid regime. Physics of the Earth and Planetary Interiors, 286, pp.138-153.

---