# Overview

Mantle model based on Driscol and Bercovici (2014), with the omission of heat loss due to melting. Evolution of a liquid basal magma ocean is included based on Davies et al. (2020), coupling it to the evolution of the solid mantle, as used in Davies and Greenwood (2022).

# Expected attributes of other regions:

Core must have the following attributes (i.e. model.core.X):
| Attribute  |            Description                                        |
|------------|------------------------------------------------------------------|
|`T_cmb`     |  Temperature at the Core-Mantle Boundary                         |

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
- Davies, C.J., Pozzo, M., Gubbins, D. and Alfè, D., 2020. Transfer of oxygen to Earth's core from a long-lived magma ocean. Earth and Planetary Science Letters, 538, p.116208.
- Davies, C.J. & Greenwood, S. 2022. Dynamics in Earth’s core arising from thermo-chemical interactions with the mantle. In: Core-Mantle coevolution – A multidisciplinary approach. Wiley. (In Press). (Accepted).
- Driscoll, P. and Bercovici, D., 2014. On the thermal and magnetic histories of Earth and Venus: Influences of melting, radioactivity, and conductivity. Physics of the Earth and Planetary Interiors, 236, pp.36-51.

---
