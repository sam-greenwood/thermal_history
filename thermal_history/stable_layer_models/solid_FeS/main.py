Description = 'Includes a solid FeS layer as a thermally stratified layer. To be used in conjunction with the solid_FeS core model.'


#List individually confirmed compatibility with other regions
compatibility = {'core': ['greenwood21'],
                 'mantle': ['driscoll_bercovici14']}


from ...core_models.leeds.routines import profiles as prof
from ...core_models.leeds.routines.snow import snow_radius
from ...core_models.solid_FeS.routines.chemistry import melting_temp_R18
from ..leeds_thermal.routines.diffusion import diffusion

import numpy as np

#Import required fuctions from leeds_thermal, no changes needed to those
from ..leeds_thermal.main import setup, evolve, update, progress, required_params
from .. import leeds_thermal as leeds_thermal

# Optional parameters. 'Name': ('Description', default_value)
optional_params = {'layer_thickness': ('Initial size of FeS layer, leave to deault value of 0', 0),
                   'max_resolution': ('Size of radial array for FeS layer. Default: 10', 10)}

#New method for solid FeS layer
def FeS_method(model):
    '''Pure thermal stratification

    Solves the thermal diffusion solution and updates the layer size during one time iteration.
    For stability, it may step in smaller time increments until the total model timestep has been
    reached.

    Parameters
    ----------
    model : ThermalModel
        Main model class

    
    '''

    core = model.core
    sl  = model.stable_layer
    prm = model.parameters

    #Read in values from model
    layer_thickness = sl.layer_thickness
    r_s = prm.r_cmb - layer_thickness
    r_s_original = r_s*1

    r_cmb        = prm.r_cmb
    core_adiabat = prm.core_adiabat_params

    Q_cmb  = model.mantle.Q_cmb

    Tcen   = model.core.Tcen
    dTa_dt = model.core.dT_dt

    T_grad_cmb = Q_cmb/(-model.core.profiles['k'][-1]*4*np.pi*r_cmb**2)

    ADR = T_grad_cmb/model.core.profiles['dTa_dr'][-1]
    sl.ADR = ADR

    idx = core._snow_idx
    #Cool adiabat and find new snow radius
    Tcen = model.core.Tcen + dTa_dt*model.dt
    Ta = prof.adiabat(core.profiles['r'], Tcen, core_adiabat)
    Tm = melting_temp_R18(core.profiles['r'], core.profiles['P'],
                          core.conc_l[0]+core.dc_dt[0]*model.dt,
                          prm.core_melting_params[1:].astype('float64'))[1]

    r_snow = snow_radius(core.profiles['r'][:idx+1], Ta[:idx+1], Tm[:idx+1])
    r_s_new = r_snow


    if r_snow < prm.r_cmb:

        #Boundary conditions
        lb_type = 0
        lb = prof.adiabat(r_snow, Tcen, prm.core_adiabat_params)

        ub_type = 1
        ub = T_grad_cmb

        #New grid
        r = np.linspace(r_snow, prm.r_cmb, prm.max_resolution)
        T = prof.adiabat(r, Tcen, core_adiabat)
        T[r>core.r_snow] = np.interp(r[r>core.r_snow], core.profiles['r'], core.profiles['T'])


        #Regrid core radial profiles onto finer stable layer grid
        rho = np.interp(r, core.profiles['r'], core.profiles['rho'])
        k   = np.interp(r, core.profiles['r'], core.profiles['k'])
        cp = np.interp(r, core.profiles['r'], core.profiles['cp'])

        #Thermal diffusivity
        D_T = k/(rho*cp)

        #Conductivity gradient
        dk_dr = np.gradient(k, r, edge_order=2)

        #If whole core is stratified, lower BC should be zero flux
        if r[0] == 0:
            lb = 0

        #Calculate diffusion solution
        T_new = diffusion(T, r, model.dt, D_T, k, dk_dr, (lb_type,lb),(ub_type,ub))


        #Thermal gradient at Ts
        sl.T_grad_s = (T_new[1]-T_new[0])/(r[1]-r[0])

    else:

        r = np.full(prm.max_resolution, prm.r_cmb)
        T_new = np.full(prm.max_resolution, Ta[-1])
        sl.T_grad_s = T_grad_cmb
        lb, ub = 0, 0


    #Save new profiles. Keep original profiles until update() has been called.
    sl._next_profiles['r'] = r
    sl._next_profiles['T'] = T_new

    sl.ds_dt = (r_s_new - r_s_original)/model.dt

    sl.lb_T = lb
    sl.ub_T = ub

#Overwrite with simple FeS method.
leeds_thermal.main.pure_thermal_method = FeS_method

#New function to determine when stable layer solution is needed.
#No longer need to consider sub-adiabatic heat flows, simply when an FeS layer exists.
def pure_thermal_check(self):

    '''Determine if stable layer solution is needed

    Returns
    -------
    Bool
        Whether stable layer solution is needed this iteration
    '''

    #If FeS layer exists run method.
    if self.core.r_snow < self.parameters.r_cmb:
        return True

    return False

#Overwrite leeds_thermal version
leeds_thermal.routines.functions.pure_thermal_check = pure_thermal_check

