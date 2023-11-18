#Short description for use with thermal_history.utils.available_models()
Description = 'Thermal stratification with embedded FeS layer for Mercury. This is under development!'

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import bisect
from thermal_history.utils.optimised_funcs import trapezoid, linspace

import logging
logger = logging.getLogger(__name__)

import thermal_history as th
from .. import leeds_thermal

from ...core_models.leeds.routines import profiles as prof
from ..leeds_thermal.routines import functions as func
from .routines.diffusion import diffusion_uneven, diffusion_discont, diffusion_discont_variable

#Define main functions based on leeds_thermal
# evolve = leeds_thermal.evolve
update = leeds_thermal.update
progress = leeds_thermal.progress

def setup(model):

    sl = model.stable_layer
    core = model.core
    prm = model.parameters
    #Force initial stable layer to FeS layer size, then call leeds_thermal setup
    prm.layer_thickness = prm.FeS_size
    leeds_thermal.setup(model)

    r_fes = prm.r_cmb - prm.FeS_size

    #Initialise bulk and fes profiles
    sl.profiles['bulk'] = {'r': np.array([]),
                            'T': np.array([])}
    sl.profiles['fes'] = {'r': np.linspace(r_fes, prm.r_cmb, 10)} #This will remain fixed throughout simulation.

    #Adiabatic initial temperature
    sl.profiles['fes']['T'] = prof.adiabat(sl.profiles['fes']['r'], core.Tcen, prm.core_adiabat_params)
    model.core.Q_fes = core.profiles['Qa'][core._fes_idx-1]


    #Reset profiles with new temperature profile
    prof.basic_profiles(model)
    prof.temp_dependent_profiles(model)

    #New boolean user can set to fix layer size constant in time.
    if not hasattr(sl, 'hold'):
        sl.hold = False

#Add on extra required parameters
required_params = leeds_thermal.required_params
required_params['FeS_size']         = 'Size of the FeS layer'
required_params['FeS_density']      = 'Uniform density of liquid FeS'
required_params['FeS_cp']           = 'Uniform specific heat of liquid FeS'
required_params['FeS_conductivity'] = 'Uniform thermal conductivity of liquid FeS'

optional_params = leeds_thermal.optional_params


#!# evolve method to account for FeS layer that is thermally conducting
def evolve(model):
    '''Evolve the stable layer, taking into account the presence of a liquid FeS layer.
    Changes from leeds_thermal.main.evolve are commented by #!#

    Parameters
    ----------
    model : ThermalModel
        Main model class
    '''

    sl = model.stable_layer
    core = model.core
    prm = model.parameters

    sl.Q_cmb = model.mantle.Q_cmb

    fes_idx = core._fes_idx           #!# FeS index for core radial arrays
    r_fes = prm.r_cmb - prm.FeS_size

    T_grad_cmb = sl.Q_cmb/(-core.profiles['k'][fes_idx-1]*4*np.pi*prm.r_cmb**2)
    ADR = T_grad_cmb/core.profiles['dTa_dr'][-1]
    sl.ADR = ADR

    #Check if conditions are sub-adiabatic and a layer should begin to grow
    check  = func.pure_thermal_check(model)
    check = True #!# Always a least an FeS layer.

    #Initialise energies/entropies associated with this state.
    Qs, Es, Ek, Ej = 0, 0, 0, 0

    #initialise BC's
    sl.lb_T, sl.ub_T = 0,0
    sl.T_grad_s = 0

    sl.alpha_D = 0
    sl.baro_grad = 0

    #Run the method.
    r_initial = sl.profiles['r']
    T_initial = sl.profiles['T']

    #*# FeS profiles
    FeS_prof = sl.profiles['fes']
    Bulk_prof = sl.profiles['bulk']

    #Calculate Ek
    #*# FeS layer contribution. Doesn't always have even grid spacing across bulk+FeS layer.
    dT_dr = np.gradient(FeS_prof['T'], FeS_prof['r'], edge_order=2)
    Ek = 4*np.pi*trapezoid(FeS_prof['r'], prm.FeS_conductivity * (dT_dr/FeS_prof['T'])**2 * FeS_prof['r']**2)[-1]

    if core.rs < r_fes:
        #*# Contribution from stratified bulk
        dT_dr = np.gradient(Bulk_prof['T'], Bulk_prof['r'], edge_order=2)
        k     = np.interp(Bulk_prof['r'], core.profiles['r'][:fes_idx], core.profiles['k'][:fes_idx])
        Ek   += 4*np.pi*trapezoid(Bulk_prof['r'], k * (dT_dr/Bulk_prof['T'])**2 * Bulk_prof['r']**2)[-1]

    #Time step model
    conducting_FeS(model)  #!# Changed to our new method

    #New radial profiles
    r_new = sl._next_profiles['r']
    T_new = sl._next_profiles['T']

    #Secular cooling. Compare final profiles to those at the beginning of the time-step.
    if r_new[0] <= r_initial[0]:  #Layer has grown

        #Need to append on adiabat values onto initial temperatures
        T_initial_temp = prof.adiabat(r_new, core.Tcen, prm.core_adiabat_params)
        for i in range(r_new.size):
            if r_new[i] >= r_initial[0]:
                T_initial_temp[i] = np.interp(r_new[i], r_initial, T_initial)
        T_initial = T_initial_temp

        r1, T2, T1 = r_new, T_new, T_initial

    #Set cooling rate at r=0 if no convecting region exists.
        if r_initial[0] == 0:
            core.dT_dt = (T2[0]-T1[0])/model.dt

    elif r_new[0] > r_initial[0]:  #Layer has shrunk

        #Need to append on adiabatic values onto new temperature
        T_new_temp = prof.adiabat(r_initial, core.Tcen + core.dT_dt*model.dt, prm.core_adiabat_params)
        for i in range(r_initial.size):
            if r_initial[i] >= r_new[0]:
                T_new_temp[i] = np.interp(r_initial[i], r_new, T_new)

        T_new = T_new_temp

        r1, T2, T1 = r_initial, T_new, T_initial

    dT_dt = (T2-T1)/model.dt
    
    cp  = np.interp(r1, core.profiles['r'], core.profiles['cp']) #!# Just consider values below FeS
    rho = np.interp(r1, core.profiles['r'], core.profiles['rho']) 
    Qs = -4*np.pi*trapezoid(r1, r1**2*dT_dt*rho*cp)[-1]
    Es = -4*np.pi*trapezoid(r1, r1**2*dT_dt*rho*cp*(1/T1[-1] - 1/T1))[-1]
    sl.dT_dt_s = dT_dt[0]  #Save cooling rate at base of layer

    Ej = Es - Ek

    #Save energies/entropies to model attributes
    sl.Qs, sl.Es, sl.Ek, sl.Ej = Qs, Es, Ek, Ej

    #Add onto core results
    core.Qs += Qs
    core.Ej += Ej
    core.Es += Es
    core.Ek += Ek

    #Make sure profiles have size maximum_resolution, so they can be saved.
    n = sl.profiles['r'].size

    #!# For now we won't worry about profiles at each time step having the same size.
    if n > prm.max_resolution:
        pass
        # r = linspace(sl.profiles['r'][0], sl.profiles['r'][-1], prm.max_resolution)
        # T = np.interp(r, sl.profiles['r'], sl.profiles['T'])

    elif n < prm.max_resolution:
        pass
        # r = np.append(sl.profiles['r'], np.ones(prm.max_resolution-n)*sl.profiles['r'][-1])
        # T = np.append(sl.profiles['T'], np.ones(prm.max_resolution-n)*sl.profiles['T'][-1])

    # else:
    r, T= sl.profiles['r'], sl.profiles['T']

    sl.profiles['r'], sl.profiles['T'] = r, T

    sl.profiles['Ta'] = prof.adiabat(r, core.Tcen, prm.core_adiabat_params)

    sl.T_cmb, sl.T_s= T[-1], T[0]

    core.T_cmb,  core.T_s = sl.T_cmb, sl.T_s


#!# Thermal method with FeS
def conducting_FeS(model):
    '''Simple model taking the temperature of the FeS layer to be uniform.

    Parameters
    ----------
    model : ThermalModel
        Main model class

    '''

    #All deviations from leeds_thermal.pure_thermal_method are commented starting with #!#
    #to highlight them.


    core = model.core
    sl  = model.stable_layer
    prm = model.parameters

    #Read in values from model
    layer_thickness = sl.layer_thickness
    r_s = prm.r_cmb - layer_thickness
    r_s_original = r_s*1

    r_cmb      = prm.r_cmb
    resolution = prm.resolution
    min_res    = prm.min_resolution
    max_res    = prm.max_resolution

    TOL = prm.depth_tolerance                   # Minimum layer size
    init_size = prm.init_size + prm.FeS_size  #!# Stable layer must start below FeS layer 

    core_adiabat      = prm.core_adiabat_params
    rho_poly          = prm.core_liquid_density_params

    E_T = prm.entrainment_T

    Tcen   = model.core.Tcen
    dTa_dt = model.core.dT_dt

    #!# Base of FeS layer
    r_fes = prm.r_cmb - prm.FeS_size
    fes_idx = core._fes_idx

    #Temperature gradient at CMB
    T_grad_cmb = core.Q_cmb/(-model.core.profiles['k'][-1]*4*np.pi*r_cmb**2)

    #!# Temp in the FeS layer. 
    #!# Tracking profiles specific to FeS layer helps make sure Q_fes can be calculated
    r_prof_fes = sl.profiles['fes']['r']
    T_prof_fes = sl.profiles['fes']['T']

    #!# Temp in the FeS layer. 
    #!# Tracking profiles specific to FeS layer helps make sure Q_fes can be calculated
    r_prof_bulk = sl.profiles['bulk']['r']
    T_prof_bulk = sl.profiles['bulk']['T']
        
    core.Q_fes = 4*np.pi*r_fes**2 * prm.FeS_conductivity * -(T_prof_fes[1]-T_prof_fes[0])/(r_prof_fes[1]-r_prof_fes[0])

    #!# ADR defined at r_fes rather than r_cmb
    sl.ADR = core.Q_fes/core.profiles['Qa'][core._fes_idx-1]
    core.ADR_s = sl.ADR #Overwrite old value form core model (FeS core model always gives=1)

    #For stability, reduce the time step size when layer is relatively thin.
    factor = 1
    time_gone = 0
    while time_gone < model.dt:

        if sl.hold and layer_thickness == prm.FeS_size:
            #!# Layer should not grow and just the FeS layer should be calculated
            r, T = r_prof_fes, T_prof_fes

        #!# Only solve for FeS layer unless ADR<1
        elif sl.ADR < 1:

            
            #!# Number of grid points in the bulk this iteration
            n_points = int(np.max([min_res,layer_thickness*resolution]))
            n_points = int(np.min([max_res,n_points]))

            #Initialise layer with minium tolerance thickness
            if layer_thickness < TOL+prm.FeS_size:       #!# stable layer must start below FeS layer
                layer_thickness = init_size
                r_s = r_cmb - layer_thickness
                

                r_prof_bulk = np.linspace(r_s, r_fes, n_points) #!# only upto to grid space of r_fes. FeS profiles will get appended to these profiles.
                T_prof_bulk = prof.adiabat(r_prof_bulk, Tcen, core_adiabat)
                sl.profiles['bulk']['r'] = r_prof_bulk
                sl.profiles['bulk']['T'] = T_prof_bulk

            
            #!# Interpolate stable layer in the bulk onto new grid for this iteration.
            r_prof_bulk = np.linspace(r_s, r_fes, n_points) 
            T_prof_bulk = np.interp(r_prof_bulk, sl.profiles['bulk']['r'], sl.profiles['bulk']['T'])

            #!# Radial profiles over whole conducting region
            r = np.append(r_prof_bulk[:-1], r_prof_fes)
            T = np.append(T_prof_bulk[:-1], T_prof_fes)
            

        else:
            #!# Just the FeS layer is conducting
            r, T = r_prof_fes, T_prof_fes

        #Set timestep
        dt_small = model.dt*factor

        if time_gone + dt_small > model.dt:
            dt_small = model.dt - time_gone

        #Regrid core radial profiles onto finer stable layer grid
        rho   = prof.density(r, rho_poly)
        k     = np.interp(r, core.profiles['r'], core.profiles['k'])
        dk_dr = np.gradient(k, r, edge_order=2)
        cp    = np.interp(r, core.profiles['r'], core.profiles['cp'])
        alpha = np.interp(r, core.profiles['r'], core.profiles['alpha'])

        # Thermal diffusivity
        D_T = k/(rho*cp)

        #Cool adiabat
        Tcen = model.core.Tcen + dTa_dt*(time_gone+dt_small)

        Ta, Ta_grad = prof.adiabat(r, Tcen, core_adiabat), prof.adiabat_grad(r, Tcen, core_adiabat)

        #!# Fixed value if bulk is convecting
        if r[0] == r_fes:
            lb = Ta[0]
            lb_type = 0
        else:
            #Lower boundary condition (using entrainment)
            lb = (1-E_T)*Ta_grad[0]
            # lb = prof.adiabat_grad(r, core.Tcen, core_adiabat)[0]
            lb_type = 1

        #If whole core is stratified, lower BC should be zero flux
        if r[0] == 0:
            lb = 0

        #Upper bc based on Q_cmb
        ub_type, ub = 1, T_grad_cmb

        #Calculate diffusion solution
        #!# Check if standard diffusion solution is needed or new method with discontinuity in k
        if r[0]==r_fes:
            #Just FeS layer. diffusion_uneven() gives same results as diffusion() from leeds_thermal if
            #uniform grid is used. Can support an uneven grid spacing if needed.
            T_new = diffusion_uneven(T, r, dt_small, D_T, k, dk_dr, (lb_type,lb),(ub_type, ub)) #!# Uneven grid.

        else:

            #!# For diffusion solution with a discontinuity, constant material properties in each
            #!# region are assumed for now.

            #Now we allow variable k/D
            #Check values are constant
            # for x in [rho, cp, k]:
            #     if np.max(np.abs(np.diff(x[r<r_fes]))) > 0:
            #         logger.critical('Diffusion solution with discontinuity in k assumes rho, cp, k are constant in radius in bulk')

            #     if np.max(np.abs(np.diff(x[r>r_fes]))) > 0:
            #         logger.critical('Diffusion solution with discontinuity in k assumes rho, cp, k are constant in radius in FeS')

            #!# Thermal diffusivity
            # D_lower = k[0]/(rho[0]*cp[0])
            # D_upper = k[-1]/(rho[0]*cp[0])

            # k_lower = k[0]
            # k_upper = k[-1]

            #!# Diffusion solution for a discontinuity in k.
            # T_new = diffusion_discont(T, r, r_fes, dt_small, D_lower, D_upper, k_lower, k_upper, (lb_type, lb), (ub_type, ub))

            T_new = diffusion_discont_variable(T, r, r_fes, dt_small, D_T, k, dk_dr, (lb_type, lb), (ub_type, ub))


        #!# Only need to change layer thickness if stable layer will grow/receed.
        if sl.ADR < 1 and layer_thickness > prm.FeS_size:
            #Mix layer (experimental function, not used by default)
            if prm.mix_layer:
                T_rel = func.mix_profile(r, cp*rho*r**2, T_new-Ta)
                T_new = T_rel + Ta

            #Determine if layer should retreat
            #!# Just consider areas below FeS layer
            #*# New boolean model.stable_layer.hold allows layer size to be fixed manually.
            if sl.hold:
                r_s_new = r[0]
            else:
                r_s_new = func.retreat_layer(r[r<r_fes], T_new[r<r_fes], np.zeros(r.size)[r<r_fes], Tcen, 0, core_adiabat, (resolution,min_res,max_res), alpha[r<r_fes], 0)

            #while core is stratified
            if r[0] == 0 and sl.ADR < 1:
                r_s_new = 0

            #Grow layer
            elif r_s_new == r[0]:

                #Find radius on adiabat that matches temperature at base of layer
                def f(guess,T,Tcen):
                    return T - prof.adiabat(guess,Tcen,core_adiabat)
                
                if T_new[0] >= Tcen: #Layer has reached center of core
                    r_s_new = 0
                else:
                    #!# Limit search to 0 -> r_fes
                    r_s_new = bisect(f,0,r_fes,args=(T_new[0],Tcen),maxiter=200) 

            #Shrink layer
            else:
                #Give up if layer will still not grow after generous time step increases.
                if r_s_new >= r_fes-TOL:     #!# replaced r_cmb with r_fes

                    #No stratification
                    r_s_new = r_fes        #!# replaced r_cmb with r_fes
                    factor = factor*2      
                    time_gone =- dt_small  #¡# Reset timer back to start

                    if factor > 64 or sl.ADR >= 1:
                        dt_small = model.dt - time_gone
                        factor = 1

            #Ensure stable layer doesn't go below ICB
            ri_new = core.ri + core.dri_dt*model.dt
            if r_s_new < ri_new:
                r_s_new = np.min([ri_new+10e3, prm.r_cmb]) #At least 10km above ICB.

            layer_thickness = r_cmb - r_s_new

            #Regrid domain, keeping temperature relative to the adiabat.
            T_rel = T_new[r<r_fes] - prof.adiabat(r[r<r_fes], Tcen, core_adiabat)

            r_prof_bulk, T_rel = func.change_domain_size(T_rel, r[r<r_fes], r_s_new, (resolution,min_res,max_res))

            T_prof_bulk = T_rel + prof.adiabat(r_prof_bulk, Tcen, core_adiabat)


            r, T = np.append(r_prof_bulk, r_prof_fes), np.append(T_prof_bulk, T_new[r>=r_fes])

        else:
            r_s_new = r_fes
            T = T_new

        time_gone += dt_small
  


    #Save new profiles. Keep original profiles until update() has been called.
    sl._next_profiles['r'] = r
    sl._next_profiles['T'] = T

    #!# Save separate bulk and FeS profiles.
    sl.profiles['bulk']['r'] = r[:-r_prof_fes.size]
    sl.profiles['bulk']['T'] = T[:-r_prof_fes.size]

    sl.profiles['fes']['r'] = r[-r_prof_fes.size:]
    sl.profiles['fes']['T'] = T[-r_prof_fes.size:]

    sl.ds_dt = (r_s_new - r_s_original)/model.dt

    if r[1]-r[0] > 0:
        sl.T_grad_s = (T[1]-T[0])/(r[1]-r[0])
    else:
        sl.T_grad_s = 0

    sl.lb_T = lb
    sl.ub_T = ub




