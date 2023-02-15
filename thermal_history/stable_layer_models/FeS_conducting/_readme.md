# Overview:

Development model for a liquid FeS layer embedded in a bulk Fe-Si core. The model assumes the FeS layer is always conducting and the bulk core can also thermally stratify beneath it. Note the core model `FeS` must also be used in conjunction with this stable layer model.

When the bulk is fully convecting, the conduction solution for the FeS layer is fixed to the adiabat at it's lower boundary. The heat flow from the bulk into the FeS layer is then calculated from the thermal gradient at the base of the FeS layer, used for calculating the cooling rate of the core in the next iteration:

$Q_\mathrm{FeS} = -k_\mathrm{FeS} 4\pi r_\mathrm{FeS}^2 \nabla T_\mathrm{FeS}^+$

The superscript + indicates that the gradient is taken on the upper side of the boundary between bulk and FeS layer.

When the thermal stratification reaches deeper into the bulk, the boundary condition at the base of the thermal layer and the method for advancing the layer deeper into the core is the same as the `leeds_thermal` method.


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
model = setup_model(prm, core_method='FeS', stable_layer_method='FeS')

#Evolve model
dt = 1e6*prm.ys
for i in range(1000):
    model.mantle.Q_cmb = 1e12 #Set CMB heat flow.
    model.evolve(dt)
    
```


## Note some key modelling assumptions:
1. For the conducting method, the diffusion solution for calculating the temperature profile when a discontinuity in thermal conductivity also assumes that the density, specific heat and thermal conductivity are constant in radius for each region (bulk and FeS layer). A critical warning will be logged is this is not adhered to with the input parameters. This restriction should be relatively easy to relax, just need the discretisation scheme in `.routines.diffusion.diffusion_discont` to account for it as `.routines.diffusion.diffusion_uneven` does.
---