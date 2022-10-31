# Overview:

Core thermal evolution model first used in this form in Greenwood et al. (2021). The model is based on ongoing research at the University of Leeds and their collaborators, building upon the work by (but not exclusive to) Gubbins et al. (2004) and Davies (2015).

This model calculates the energy and entropy budgets for the the thermo-chemical evolution of an adiabatic core. Support for stably stratified layers is built in and works with the `leeds_thermal` and `leeds_chemical` stable layer models. The iron snow model of Davies and Pommier (2018) is also built into this core model, enabled by the parameter `iron_snow=True`. The user is assumed to have established the freezing regime of the core a priori allowing them to prescribe a bottom-up or top-down freezing regime.

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
model = setup_model(prm, core_method='leeds')

#Evolve model
dt = 1e6*prm.ys
for i in range(1000):
    model.mantle.Q_cmb = 1e12 #Set CMB heat flow.
    model.evolve(dt)
    
```

## References
- Davies, C. (2015). Cooling history of Earth’s core with high thermal conductivity. Physics of the Earth and Planetary Interiors, 247, pp. 65-79.
- Greenwood, S., Davies, C. J., & Mound, J. E. (2021). On the evolution of thermally stratified layers at the top of Earth’s core. Physics of the Earth and Planetary Interiors, 106763.
- Gubbins, D., Alfe, D., Masters, G., Price, G., & Gillan, M. (2003). Gross thermodynamics of two-component core convection. Geophysical Journal International, 157, pp. 1407-1414

---