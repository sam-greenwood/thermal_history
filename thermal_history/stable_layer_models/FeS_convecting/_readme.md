# Overview:

Development model for a liquid FeS layer embedded in a bulk Fe-Si core. The model assumes the FeS layer is always convecting and that the bulk core thermally stratify beneath it. Note the core model `FeS` must also be used in conjunction with this stable layer model. The FeS layer consequently always has an adiabatic temperature profile.

When the bulk is fully convecting, the whole core is convecting in a 2 layered system. The interface between the bulk and FeS layer cool at the same rate and so the method is the same as for a homogenous core, just taking into account the material property changes in the FeS in the whole core integrals. The only heat source in the FeS layer is the secular cooling ($Q_\mathrm{s, FeS}$) and so with the cooling rate known, the heat flow at the base of the FeS layer ($Q_\mathrm{FeS}$) is the difference between $Q_\mathrm{cmb}$ and $Q_\mathrm{s, FeS}$.

When $Q_\mathrm{FeS}$ falls below the adiabatic heat flow at that radius, the bulk begins to stratify. The heat flow from the bulk into the FeS layer is then calculated from the thermal gradient at the top of the stratified bulk:

$Q_\mathrm{FeS} = -k^- 4\pi r_\mathrm{FeS}^2 \nabla T_\mathrm{FeS}^-$

The superscript - indicates the value is taken on the lower side of the boundary between bulk and FeS layer.

The boundary condition at the base of the thermal layer and the method for advancing the layer deeper into the core is the same as the `leeds_thermal` method. The upper boundary condition is fixed temperature to the base of the FeS adiabat. The FeS adiabat is cooled by calculating $Q_\mathrm{s, FeS} = Q_\mathrm{cmb} - Q_\mathrm{FeS}$. On the first timestep when the bulk becomes stratified, an assumption is made on the FeS cooling rate to ensure the layer can begin to grow (see assumptions list below).

An additional toggle has been added: model.stable_layer.hold. When this is True, the stable layer size is held to the present constant radius and cannot change. This can be useful on the first iteration when some mantle models need to instantaneously adjust from the initial conditions and can yeild a low heat flow promoting a stable layer in the core. The default value is False. Example of use:

```python
from thermal_history.model import Parameters, setup_model

#Load parameters and setup model
prm = Parameters('/path/to/parameters_file')
model = setup_model(prm, core_method='FeS', stable_layer_method='FeS_convecting')

model.stable_layer.hold = True #Hold growth of stable layer

#Evolve model
dt = 1e6*prm.ys
for i in range(1000):

     #Allow stable layer to grow after 1st iteration
    if i > 1:
        model.stable_layer.hold = False
        
    model.mantle.Q_cmb = 1e12 #Set CMB heat flow.
    model.evolve(dt)
    
```


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

2. When the bulk first becomes stratified, due to the boundary conditions chosen, a layer would not initially start to grow with the method above. Instead, only on the iteration when a layer in the bulk is initialised, the cooling rate of the FeS layer is ignored and taken to be zero for the purpose of setting the boundary conditions. Since the FeS layer and bulk cool in tandem when the core is superadiabatic, when $Q_\mathrm{FeS}$ first falls below the adiabatic value, the upper boundary condition is the same as the adiabatic temperature in the bulk. This leads to a diffusion profile in the thin stratified layer in the bulk that is colder than the adiabat (except at the upper boundary), preventing the layer to form. By ignoring the FeS cooling rate on this first iteration, it ensures a stabilising temperature profile and from then on can continue to grow with the correct FeS cooling rate.
---