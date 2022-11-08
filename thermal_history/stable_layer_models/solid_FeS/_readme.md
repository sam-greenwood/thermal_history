# Overview:

Stable layer model for thermal stratification of Greenwood et al. (2021). In this model only the temperature structure of the stratified layer is used to evolve its size over time. This model must be used in conjunction with a core model.

# Expected attributes of other regions:

Core must have the following attributes (i.e. model.core.X):
| Attribute  |  Description                                                     |
|------------|------------------------------------------------------------------|
|`rs`        |  Radius of the stable layer interface                            |
|`Tcen`      |  The core temperature at its center                              |
|`dT_dt`     |  Rate of change of 'Tcen'                                        |
|`profiles`  |  A dictionary of radial core profiles (see leeds core model)     |
|`Qs`        |  Secular cooling to add on the contribution from the stable layer|
|`Es, Ej, Ek`|  Entropy terms to add on the contribution from the stable layer  |

Mantle must have the following attributes (i.e. model.mantle.X):
| Attribute     |            Description                                           |
|---------------|------------------------------------------------------------------|
|`Q_cmb`        |  Heat flow across Core-Mantle Boundary                           |

If not using a mantle model, `Q_cmb` must be set manually, e.g.:

```python
from thermal_history.model import Parameters, setup_model

#Load parameters and setup model
prm = Parameters('/path/to/parameters_file')
model = setup_model(prm, core_method='leeds', stable_layer_method='leeds_thermal')

#Evolve model
dt = 1e6*prm.ys
for i in range(1000):
    model.mantle.Q_cmb = 1e12 #Set CMB heat flow.
    model.evolve(dt)
    
```


## Note some key modelling assumptions:
1. Compositional changes in the core do not directly influence growth of a stratified layer.
2. Once the entire core is thermally stable, only when the CMB heat flow exceeds the adiabatic heat flow at the CMB will the core convect again. This should be reasonable for a constant radial core conductivity profile, not necessarily true otherwise.

## References:
- Greenwood, S., Davies, C. J., & Mound, J. E. (2021). On the evolution of thermally stratified layers at the top of Earth’s core. Physics of the Earth and Planetary Interiors, 106763.


---