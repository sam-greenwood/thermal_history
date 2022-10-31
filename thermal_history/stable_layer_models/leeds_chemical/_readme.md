# Overview:

Stable layer model for chemical stratification of Davies and Greenwood (2022). This model uses the method of Buffett and Seagle (2010) to grow the layer, with extensions to calculate the energy/entropy balance and couple to the rest of the core. The temperature profile is also found by solving the diffusion solution, with the thermal gradient at `rs` controlling the lower chemical boundary condition. Whilst multiple light elements may be present in the iron alloy, only 1 is treated as controlling the stable layer (the light element represented by the first value in the `conc_l` array). Others are assumed to match their uniform values in the bulk. This model must be used in conjunction with a core model.

# Expected attributes of other regions:

Core must have the following attributes (i.e. model.core.X
| Attribute       |            Description                                           |
|-----------------|------------------------------------------------------------------|
|`rs`             |  Radius of the stable layer interface                            |
|`Tcen`           |  The core temperature at its center                              |
|`dT_dt`          |  Rate of change of 'Tcen'                                        |
|`conc_l`         |  Array of light element mass fractions                           |
|`dc_dt`          |  Rate of change of 'conc_l'                                      |
|`profiles`       |  A dictionary of radial core profiles                            |
|`Qs`             |  Secular cooling to add on the contribution from the stable layer|
|`Es, Ej, Ek, 'Ea`|  Entropy terms to add on the contribution from the stable layer  |

Mantle must have the following attributes (i.e. model.mantle.X):
| Attribute       |            Description                                           |
|-----------------|------------------------------------------------------------------|
|`Q_cmb`          |  Heat flow across Core-Mantle Boundary                           |
|`chemical_bc_cmb`|  Chemical boundary conditions for core stable layer solution     |

If not using a mantle model, `Q_cmb` and `chemical_bc_cmb` must be set manually, e.g.:

```python
from thermal_history.model import Parameters, setup_model

#Load parameters and setup model
prm = Parameters('/path/to/parameters_file')
model = setup_model(prm, core_method='leeds', stable_layer_method='leeds_chemical')

#Evolve model
dt = 1e6*prm.ys
for i in range(1000):
    model.mantle.Q_cmb = 1e12 #Set CMB heat flow.
    #Chemical boundary conditions (1/0, float).
    #1/0 indicates fixed flux or fixed value, then provide the value of the boundary condition.
    model.mantle.chemical_bc_cmb = (1, 0)   #Zero flux BC.
    model.evolve(dt)
    
```

## Note some key modelling assumptions:
1. A minimum chemical gradient is assumed at `rs` of $10^-10 \mathrm{m}^{-1}$ to ensure numerical stability. Lower values than this threshold arrise from heat flows approaching or below the adiabatic heat flow at the CMB. In such a case, thermal stratification can begin to grow and would likely outpace the chemical layer. A chemical layer embedded within a thermal layer may then be the appropriate core structure. There is a stable layer method, `leeds_thermochemical', under development that attempts to account for such a scenario. However, it to date has not been used in any publication.
2. As mentioned above, only 1 light element is assumed to diffuse in the stable layer. The light element represented by the first value in the `model.core.conc_l` array (and subsequently any other arrays relating to light elements e.g. `dc_dt`) is solved for in the stable layer. All others are equal to the homogenous interior values.

## References:

- Buffett, B.A. and Seagle, C.T., 2010. Stratification of the top of the core due to chemical interactions with the mantle. Journal of Geophysical Research: Solid Earth, 115(B4).
- Davies, C.J. & Greenwood, S. 2022. Dynamics in Earth’s core arising from thermo-chemical interactions with the mantle. In: Core-Mantle coevolution – A multidisciplinary approach. Wiley. (In Press). (Accepted).

---