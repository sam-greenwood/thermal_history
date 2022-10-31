# Overview:

Adaptation of the `leeds` core model with the iron snow routines changed to calculate a freezing FeS layer at the CMB as per Rückriemen et al. (2018) (enabled by the parameter `iron_snow=True`).

# Expected attributes of other regions:

Mantle must have the following attributes (i.e. model.mantle.X):
| Attribute     |            Description                                           |
|---------------|------------------------------------------------------------------|
|`Q_cmb`        |  Heat flow across Core-Mantle Boundary                           |

If not using a mantle model, `Q_cmb` must be set manually, e.g.:

```python
from thermal_history.model import Parameters, setup_model

#Load parameters and setup model
prm = Parameters('/path/to/parameters_file')
model = setup_model(prm, core_method='solid_FeS')

#Evolve model
dt = 1e6*prm.ys
for i in range(1000):
    model.mantle.Q_cmb = 1e12 #Set CMB heat flow.
    model.evolve(dt)
    
```


## References:
- Rückriemen, T., Breuer, D. and Spohn, T., 2018. Top-down freezing in a Fe–FeS core and Ganymede’s present-day magnetic field. Icarus, 307, pp.172-196.

---