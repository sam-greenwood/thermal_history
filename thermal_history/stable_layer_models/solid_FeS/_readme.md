# Overview:

Stable layer model for a solid FeS layer to be used in conjunction with the `solid_FeS` core model. This adds a conduction solution to the solid FeS layer but does not grow thermal stratification any deeper into the core. To grow stratification into the core the `leeds_thermal` stable layer model should be used. This model is based upon the `leeds_thermal` and so imports much of its code directly from that model.

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
model = setup_model(prm, core_method='solid_FeS', stable_layer_method='solid_FeS')

#Evolve model
dt = 1e6*prm.ys
for i in range(1000):
    model.mantle.Q_cmb = 1e12 #Set CMB heat flow.
    model.evolve(dt)
    
```


## Note some key modelling assumptions:
1. Only the solid FeS region is taken to be conducting.



---