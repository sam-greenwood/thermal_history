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
from .routines.diffusion import diffusion_uneven, diffusion_discont

#Define main functions based on leeds_thermal
# evolve = leeds_thermal.evolve
update = leeds_thermal.update
progress = leeds_thermal.progress

def setup(model):

    sl = model.stable_layer
    prm = model.parameters
    #Force initial stable layer to FeS layer size, then call leeds_thermal setup
    prm.layer_thickness = prm.FeS_size
    leeds_thermal.setup(model)

    #Set relevant method
    if prm.FeS_method == 'isothermal':
        th.stable_layer_models.FeS.main.evolve = simple_isothermal_evolve

        profiles = sl.profiles
        profiles['T'] = np.full(profiles['r'].size, model.core.T_cmb) #Isothermal initial state

    elif prm.FeS_method == 'conducting':
        th.stable_layer_models.FeS.main.evolve = conducting_evolve

        r_fes = prm.r_cmb - prm.FeS_size

        #Initialise bulk and fes profiles
        sl.profiles['bulk'] = {'r': np.array([]),
                               'T': np.array([])}
        sl.profiles['fes'] = {'r': np.linspace(r_fes, prm.r_cmb, 3)} #This will remain fixed throughout simulation.

        sl.profiles['fes']['T'] = prof.adiabat(sl.profiles['fes']['r'], model.core.Tcen, prm.core_adiabat_params)

    else:
        raise ValueError('Need to correctly specify method with FeS_method parameter')

    #Reset profiles with new temperature profile
    prof.basic_profiles(model)
    prof.temp_dependent_profiles(model)

#Add on extra required parameters
required_params = leeds_thermal.required_params
required_params['FeS_size']         = 'Size of the FeS layer'
required_params['FeS_density']      = 'Uniform density of liquid FeS'
required_params['FeS_cp']           = 'Uniform specific heat of liquid FeS'
required_params['FeS_conductivity'] = 'Uniform thermal conductivity of liquid FeS'
required_params['FeS_method']       = 'String denoting method. Can be: \'isothermal\', \'conducting\''

optional_params = leeds_thermal.optional_params



### Redefine evolve method to account for FeS layer that is isothermal
def simple_isothermal_evolve(model):
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

    #Check if Qcmb would give a layer this time-step.
    if sl.layer_thickness == 0 and not check:

        sl.ds_dt = 0
        sl.dT_dt_s = 0
        sl._next_profiles['r'] = np.ones(prm.max_resolution)*prm.r_cmb
        sl._next_profiles['T'] = prof.adiabat(sl._next_profiles['r'], core.Tcen, prm.core_adiabat_params)

    #Otherwise run the method.
    else:

        r_initial = sl.profiles['r']
        T_initial = sl.profiles['T']
        dr = r_initial[1]-r_initial[0]


        #Calculate Ek
        if dr == 0: #Layer is still of 0 thickness
            Ek == 0
        else:
            dT_dr = np.gradient(T_initial, dr, edge_order=2)
            k = np.interp(r_initial, core.profiles['r'], core.profiles['k'])
            Ek = 4*np.pi*trapezoid(r_initial, k*(dT_dr/T_initial)**2*r_initial**2)[-1]


        #Time step model
        simple_isothermal_FeS(model)  #!# Changed to our new method

        #New radial profiles
        r_new = sl._next_profiles['r']
        T_new = sl._next_profiles['T']

        T_new = T_new[r_new <= r_fes] #!# Just stable layer of the bulk, no FeS
        r_new = r_new[r_new <= r_fes]
        
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
        
        cp  = np.interp(r1, core.profiles['r'][:fes_idx], core.profiles['cp'][:fes_idx]) #!# Just consider values below FeS
        rho = np.interp(r1, core.profiles['r'][:fes_idx], core.profiles['rho'][:fes_idx]) 
        Qs = -4*np.pi*trapezoid(r1, r1**2*dT_dt*rho*cp)[-1]
        Es = -4*np.pi*trapezoid(r1, r1**2*dT_dt*rho*cp*(1/T1[-1] - 1/T1))[-1]
        sl.dT_dt_s = dT_dt[0]  #Save cooling rate at base of layer


        #!# Add on values for FeS cooling
        Qs += -(4/3)*np.pi*(prm.r_cmb**3 - r_fes**3)*prm.FeS_density*prm.FeS_cp * sl.dT_dt_fes
        Es += 0 #Operates @ T_cmb so no thermodynamic efficiency

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

    if n > prm.max_resolution:

        r = linspace(sl.profiles['r'][0], sl.profiles['r'][-1], prm.max_resolution)
        T = np.interp(r, sl.profiles['r'], sl.profiles['T'])

    elif n < prm.max_resolution:

        r = np.append(sl.profiles['r'], np.ones(prm.max_resolution-n)*sl.profiles['r'][-1])
        T = np.append(sl.profiles['T'], np.ones(prm.max_resolution-n)*sl.profiles['T'][-1])

    else:
        r, T= sl.profiles['r'], sl.profiles['T']

    sl.profiles['r'], sl.profiles['T'] = r, T

    sl.profiles['Ta'] = prof.adiabat(r, core.Tcen, prm.core_adiabat_params)

    sl.T_cmb, sl.T_s= T[-1], T[0]

    core.T_cmb,  core.T_s = sl.T_cmb, sl.T_s


#Thermal method with FeS
def simple_isothermal_FeS(model):
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

    #!# Existing radial profiles excluding FeS layer
    r = sl.profiles['r'][sl.profiles['r'] <= r_fes]
    T = sl.profiles['T'][sl.profiles['r'] <= r_fes]

    #!# Evolve temperature of FeS layer
    if layer_thickness > prm.FeS_size:
        k_fes = np.interp(r_fes-1, core.profiles['r'], core.profiles['k']) #Get value on lower side of boundary
        T_grad_fes = (T[-1]-T[-2])/(r[-1]-r[-2])
        core.Q_fes = -4 * np.pi *r_fes**2 * k_fes * T_grad_fes

        dQ = core.Q_cmb - core.Q_fes
    else:
        dQ = core.Q_cmb - core.Q_rs

    #Cooling rate of FeS layer
    Qs_tilde = -prm.FeS_density*prm.FeS_cp * (4/3)*np.pi*(prm.r_cmb**3 - r_fes**3) #* (core.T_cmb/core.Tcen)
    dT_dt_fes = dQ/Qs_tilde
    sl.dT_dt_fes = dT_dt_fes

    #!# ADR defined at r_fes rather than r_cmb
    sl.ADR = core.Q_fes/core.profiles['Qa'][core._fes_idx-1]

    #!#Only run if stable layer needs to calculated.
    #!#Skip the first iteration as mantle model often gives small heat flows before jumping up again. This messes up the thermal evolution.
    if model.it > 1 and (sl.ADR <= 1 or layer_thickness > prm.FeS_size):

        #For stability, reduce the time step size when layer is relatively thin.
        factor = 1
        time_gone = 0
        while time_gone < model.dt:


            #Initialise layer with minium tolerance thickness
            if layer_thickness < TOL+prm.FeS_size:       #!# stable layer must start below FeS layer

                layer_thickness = init_size
                r_s = r_cmb - layer_thickness

                #Number of grid points in this iteration
                n_points = int(np.max([min_res,layer_thickness*resolution]))
                n_points = int(np.min([max_res,n_points]))

                r = np.linspace(r_s, r_fes, n_points) #!# only upto r_fes
                T = prof.adiabat(r, Tcen, core_adiabat)


            #Set timestep
            dt_small = model.dt*factor

            if time_gone + dt_small > model.dt:
                dt_small = model.dt - time_gone

            fes_idx = core._fes_idx #!# We want values below FeS layer
            #Regrid core radial profiles onto finer stable layer grid
            rho = prof.density(r, rho_poly)
            k   = np.interp(r, core.profiles['r'][:fes_idx], core.profiles['k'][:fes_idx])
            dk_dr = np.gradient(k, r, edge_order=2)
            cp = np.interp(r, core.profiles['r'][:fes_idx], core.profiles['cp'][:fes_idx])
            alpha = np.interp(r, core.profiles['r'][:fes_idx], core.profiles['alpha'][:fes_idx])

            #Thermal diffusivity
            D_T = k/(rho*cp)

            #Cool adiabat
            Tcen = model.core.Tcen + dTa_dt*(time_gone+dt_small)

            Ta, Ta_grad = prof.adiabat(r, Tcen, core_adiabat), prof.adiabat_grad(r, Tcen, core_adiabat)

            #Lower boundary condition (using entrainment)
            lb = (1-E_T)*Ta_grad[0]
            # lb = prof.adiabat_grad(r, core.Tcen, core_adiabat)[0]
            lb_type = 1

            #If whole core is stratified, lower BC should be zero flux
            if r[0] == 0:
                lb = 0

            #!# Fixed T upper boundary condition
            ub_type, ub = 0, core.T_cmb + dT_dt_fes*(time_gone+dt_small)

            #Calculate diffusion solution
            T_new = diffusion(T, r, dt_small, D_T, k, dk_dr, (lb_type,lb),(ub_type, ub))

            # plt.plot(r, T_new-Ta, 'o-')
            # plt.show()
            # breakpoint()


            #Mix layer (experimental function, not used by default)
            if prm.mix_layer:
                T_rel = func.mix_profile(r, cp*rho*r**2, T_new-Ta)
                T_new = T_rel + Ta

            #Determine if layer should retreat
            r_s_new = func.retreat_layer(r, T_new, np.zeros(r.size), Tcen, 0, core_adiabat, (resolution,min_res,max_res), alpha, 0)

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
                    breakpoint()
                    #No stratification
                    r_s_new = r_fes        #!# replaced r_cmb with r_fes
                    factor = factor/2      #!# works better with smaller timestep
                    time_gone =- dt_small  #¡# Reset timer back to start

                    if factor > 64 or sl.ADR >= 1:
                        dt_small = model.dt - time_gone
                        factor = 1

            layer_thickness = r_cmb - r_s_new

            #Regrid domain, keeping temperature relative to the adiabat.
            T_rel = T_new - prof.adiabat(r, Tcen, core_adiabat)

            r, T_rel = func.change_domain_size(T_rel, r, r_s_new, (resolution,min_res,max_res))

            T = T_rel + prof.adiabat(r, Tcen, core_adiabat)

            time_gone += dt_small
    else:
        #!# No stable layer growth
        lb, ub = 0, 0
        r_s_new = r_fes
        r = np.full(2, r_fes)
        T = np.full(2, core.T_cmb+dT_dt_fes*model.dt)

    #!# Add on radial profiles of FeS layer
    r_fes = np.linspace(r_fes, r_cmb, 10)[1:]
    T_fes = np.full(r_fes.size, core.T_cmb+dT_dt_fes*model.dt)
    r = np.append(r, r_fes)
    T = np.append(T, T_fes)

    #Save new profiles. Keep original profiles until update() has been called.
    sl._next_profiles['r'] = r
    sl._next_profiles['T'] = T

    # if layer_thickness > prm.FeS_size:
    #     import matplotlib.pyplot as plt
    #     plt.plot(r, T, 'ro--')
    #     # plt.show()
    #     # breakpoint()

    sl.ds_dt = (r_s_new - r_s_original)/model.dt

    if r[1]-r[0] > 0:
        sl.T_grad_s = (T[1]-T[0])/(r[1]-r[0])
    else:
        sl.T_grad_s = 0

    sl.lb_T = lb
    sl.ub_T = ub







#!# evolve method to account for FeS layer that is thermally conducting
def conducting_evolve(model):
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
    dr = r_initial[1]-r_initial[0]

    #Calculate Ek
    if dr == 0: #Layer is still of 0 thickness
        Ek == 0
    else:
        dT_dr = np.gradient(T_initial, dr, edge_order=2)
        k = np.interp(r_initial, core.profiles['r'], core.profiles['k'])
        Ek = 4*np.pi*trapezoid(r_initial, k*(dT_dr/T_initial)**2*r_initial**2)[-1]

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

    #For stability, reduce the time step size when layer is relatively thin.
    factor = 1
    time_gone = 0
    while time_gone < model.dt:

        #!# Only solve for FeS layer unless ADR<1
        if sl.ADR < 1:

            
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
        rho = prof.density(r, rho_poly)
        k   = np.interp(r, core.profiles['r'], core.profiles['k'])
        dk_dr = np.gradient(k, r, edge_order=2)
        cp = np.interp(r, core.profiles['r'], core.profiles['cp'])
        alpha = np.interp(r, core.profiles['r'], core.profiles['alpha'])

        # Thermal diffusivity
        D_T = k/(rho*cp)

        #Cool adiabat
        Tcen = model.core.Tcen + dTa_dt*(time_gone+dt_small)

        Ta, Ta_grad = prof.adiabat(r, Tcen, core_adiabat), prof.adiabat_grad(r, Tcen, core_adiabat)

        #!# Fixed value if bulk is convecting
        if sl.ADR > 1:
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

            #Check values are constant
            for x in [rho, cp, k]:
                if np.max(np.abs(np.diff(x[r<r_fes]))) > 0:
                    logger.critical('Diffusion solution with discontinuity in k assumes rho, cp, k are constant in radius in bulk')

                if np.max(np.abs(np.diff(x[r>r_fes]))) > 0:
                    logger.critical('Diffusion solution with discontinuity in k assumes rho, cp, k are constant in radius in FeS')

            #!# Thermal diffusivity
            D_lower = k[0]/(rho[0]*cp[0])
            D_upper = k[-1]/(rho[0]*cp[0])

            k_lower = k[0]
            k_upper = k[-1]

            #!# Diffusion solution for a discontinuity in k.
            T_new = diffusion_discont(T, r, r_fes, dt_small, D_lower, D_upper, k_lower, k_upper, (lb_type, lb), (ub_type, ub))


        #!# Only need to change layer thickness if stable layer will grow/receed.
        if sl.ADR < 1 and layer_thickness > prm.FeS_size:
            #Mix layer (experimental function, not used by default)
            if prm.mix_layer:
                T_rel = func.mix_profile(r, cp*rho*r**2, T_new-Ta)
                T_new = T_rel + Ta

            #Determine if layer should retreat
            #!# Just consider areas below FeS layer
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






















# def setup(model):
#     '''Setup the initial stable layer

#     Parameters
#     ----------
#     model : ThermalModel
#         Main model class
#     '''

#     prm = model.parameters
#     sl = model.stable_layer
#     core = model.core

#     sl.layer_thickness = prm.FeS_size

#     sl.layer_thickness = prm.layer_thickness

#     if not prm.core:
#         raise ValueError('Needs a core model')

#     #Initialise profiles dictionaries
#     sl.profiles = {}
#     sl._next_profiles = {}

#     n = prm.max_resolution

#     #Create radial profiles
#     sl.profiles['r'] = linspace(core.rs, prm.r_cmb, n)
#     sl.profiles['T'] = np.interp(sl.profiles['r'], core.profiles['r'], core.profiles['T'])

#     # sl.profiles['T'] = prof.adiabat(sl.profiles['r'], core._Tcen_fes, prm.fes_adiabat_params) #FeS adiabat
#     sl.T_grad_s      = prof.adiabat_grad(core.rs, core.Tcen, prm.core_adiabat_params) #Bulk core adiabatic gradient at rs

# def evolve(model):
#     '''Evolve the stable layer

#     Parameters
#     ----------
#     model : ThermalModel
#         Main model class
#     '''

#     sl = model.stable_layer
#     core = model.core
#     prm = model.parameters

#     sl.Q_cmb = model.mantle.Q_cmb

#     r_fes = prm.r_cmb - prm.FeS_size #radius of FeS layer
#     fes_idx = core._fes_idx          #Index for FeS layer in core profiles


#     #Initialise energies/entropies associated with this state.
#     Qs, Es, Ek, Ej = 0, 0, 0, 0

#     #initialise BC's
#     sl.lb_T, sl.ub_T = 0,0
#     sl.T_grad_s = 0

#     sl.alpha_D = 0
#     sl.baro_grad = 0


#     if sl.layer_thickness > prm.FeS_size:
#         T_grad_fes = np.interp(r_fes, sl.profiles['r'], np.gradient(sl.profiles['T'], sl.profiles['r'], edge_order=2))
#     else:
#         T_grad_fes = core.profiles['dTa_dr'][core._fes_idx-1]

#     heat_in = 4*np.pi*r_fes**2 * core.profiles['k'][core._fes_idx-1] * -T_grad_fes

#     #Check if FeS layer is sub-adiabatic
#     ADR_fes = (core.Q_cmb-heat_in)/core.Qa
#     ADR_fes = core.Q_cmb/core.Qa

#     FeS_convecting = True
#     if ADR_fes < 1:
#         FeS_convecting = False

#     FeS_convecting = True #Force convecting FeS layer

#     #If no thermal stratification exists yet, re-set secular cooling of bulk.
#     if sl.layer_thickness == prm.FeS_size:

#         #Cooling rates in core calculation assumes adiabatic heat flow at rs, correct it here to account for
#         #secular cooling of FeS layer as well.

#         r, rho, T, cp = core.profiles['r'], core.profiles['rho'], core.profiles['T'], core.profiles['cp']
#         Qs_tilde, Es_tilde = en.secular_cool(r, rho, T, cp, r.size, Tcmb=T[-1]) #Secular cooling of whole core (inc FeS layer)
#         Ql_tilde = core.Ql/core.dT_dt   #Re-normalise values
#         Qg_tilde = core.Qg/core.dT_dt

#         dT_dt = (core.Q_cmb-core.Qr)/(Qs_tilde+Ql_tilde+Qg_tilde)  #New corrected cooling rate of core

#         dT_dt_fes = dT_dt * core._Tcen_fes/core.Tcen  #FeS cooling rate

#         core.Q_fes = core.Q_cmb - dT_dt*(Qs_tilde - core.Qs/core.dT_dt) #Heat flow at base of FeS layer (Q_cmb minus secular cooling of FeS layer)

#         #Sub-adiabatic. reset Q_fes to Qa at r_fes and layer will start to grow.
#         if core.Q_fes < core.Qa_rs:

#             # Q_fes = core.Qa_rs #Set to adiabatic heat flow.
#             dT_dt = core.dT_dt

#             # Qs_fes = core.Q_cmb-Q_fes

#             # dT_dt_fes = Qs_fes/Qs_tilde


#             # dT_dt = (core.Q_cmb-core.Qr)/(Qs_tilde+Ql_tilde+Qg_tilde)  #New corrected cooling rate of core

#         #New ADR_fes
#         ADR_fes = core.Q_cmb/core.Qa

#         profiles = core.profiles
#         r, T = profiles['r'], profiles['T']
#         idx = core._fes_idx
#         Qs_tilde, Es_tilde = en.secular_cool(profiles['r'][idx:], profiles['rho'][idx:], profiles['T'][idx:], profiles['cp'][idx:], prm.n_profiles-idx) #Just FeS layer integrals

#         Qs_tilde = Qs_tilde * profiles['T'][idx]/core._Tcen_fes #Renormalise from Tcen value to theoretical FeS adiabat at Tcen
#         Es_tilde = Es_tilde * profiles['T'][idx]/core._Tcen_fes

#         sl.dT_dr_fes = profiles['dTa_dr'][idx]

        
            
#         # if model.it == 1 and core.Q_fes < core.Qa_rs:

#         #     sl.profiles['r'] = linspace(0, r_fes, 60)
#         #     sl.profiles['T'] = prof.adiabat(sl.profiles['r'], core.Tcen, prm.core_adiabat_params)
#         #     sl.layer_thickness = prm.r_cmb
#         #     dT_dt = 0


#         Qs += Qs_tilde*dT_dt_fes    #Add secular cooling of just FeS layer to secular cooling of total layer.
#         Es += Es_tilde*dT_dt_fes


#         #Correct core energies/entropies
#         core.Qs = core.Qs*dT_dt/core.dT_dt
#         core.Ql = core.Ql*dT_dt/core.dT_dt
#         core.Qg = core.Qg*dT_dt/core.dT_dt

#         core.dri_dt = core.Cr*core.dT_dt
#         core.dc_dt  = core.Cr*core.Cc*core.dT_dt

#         core.Es = core.Es*dT_dt/core.dT_dt
#         core.El = core.El*dT_dt/core.dT_dt
#         core.Eg = core.Eg*dT_dt/core.dT_dt

#         core.Ej = core.Es + core.Eg + core.El - core.Ek - core.Ea

#         #Correct cooling rates
#         core.dT_dt = dT_dt
#         core._dT_dt_fes = dT_dt_fes

#     elif FeS_convecting:

#         #FeS layer is convecting
        
#         profiles = core.profiles
#         r, T = profiles['r'], profiles['T']
#         idx = core._fes_idx

#         #FeS layer needs mixing as convection starts again. sl.ADR_fes is value from previous timestep.
#         if model.it > 1 and (sl.ADR_fes < 1 and ADR_fes >= 1):

#             total_heat = trapezoid(r[idx:], r[idx:]**2 * T[idx:])[-1] #Conserve total thermal mass

#             core._Tcen_fes = total_heat / trapezoid(r[idx:], r[idx:]**2 * prof.adiabat(r[idx:], 1, prm.fes_adiabat_params))[-1]

#         #Normliased energies for convecting FeS layer
#         Qs_tilde, Es_tilde = en.secular_cool(profiles['r'][idx:], profiles['rho'][idx:], profiles['T'][idx:], profiles['cp'][idx:], prm.n_profiles-idx) #Just FeS layer integrals
#         Qs_tilde *= profiles['T'][idx]/core._Tcen_fes #Renormalise from Tcen value to theoretical FeS adiabat at Tcen

#         # breakpoint()
#         # core.Q_fes = (core.Qs + core.Qg + core.Ql) + sl._Qs_strat

#         # core.Q_fes = core.Q_cmb #TEST

#         # #Heat flow into base of FeS layer
#         # T_grad_fes = np.interp(r_fes, sl.profiles['r'], np.gradient(sl.profiles['T'], sl.profiles['r'], edge_order=2))
#         T_grad_fes = (T[idx-1]-T[idx-2])/(r[idx-1]-r[idx-2])
#         core.Q_fes = 4*np.pi*r_fes**2 * profiles['k'][idx-1] * -T_grad_fes

#         core._dT_dt_fes = (core.Q_cmb - core.Q_fes)/Qs_tilde

#         sl.Qs_fes = Qs_tilde*core._dT_dt_fes

#         # if model.it >= 7574:
#         #     print(core._Tcen_fes, core._dT_dt_fes, prof.adiabat(prm.r_cmb, core._Tcen_fes, prm.fes_adiabat_params))
#         #     breakpoint()

#     else:
#         #Secular cooling of FeS layer will be accounted for in conduction solution.

#         profiles = core.profiles
#         r, T = profiles['r'], profiles['T']
#         idx = core._fes_idx

#         #Set theoretical Tcen to track the temperature at r_fes
#         core._Tcen_fes = T[idx] / prof.adiabat(r[idx], 1, prm.fes_adiabat_params)
        


#     sl.ADR_fes = ADR_fes
#     # sl.ADR_fes = 1.1 #Force adiabatic FeS layer for now


#     #Check if bulk core is adiabatic
#     ADR = core.Q_fes/core.Qa_fes
#     sl.ADR = ADR

#     #Estimate based on Qc rather than Q_fes
#     sl.ADR = core.Q_cmb / (-4*np.pi*prm.r_cmb**2 * profiles['k'][core._fes_idx-1] * prof.adiabat_grad(prm.r_cmb, core.Tcen, prm.core_adiabat_params))

#     #Profiles before diffusion solution.
#     r_initial = sl.profiles['r']
#     T_initial = sl.profiles['T']

#     #Calculate Ek
#     dT_dr = np.gradient(T_initial, r_initial, edge_order=2)
#     k = np.interp(r_initial, core.profiles['r'], core.profiles['k'])
#     Ek = 4*np.pi*trapz(k*(dT_dr/T_initial)**2*r_initial**2, x=r_initial)

#     #If no layer exists and bulk is super-adiabatic, don't bother with diffusion solution.
#     if sl.ADR >= 1 and sl.layer_thickness == prm.FeS_size:
#         sl._next_profiles['r'] = np.linspace(r_fes, prm.r_cmb, 10)
#         sl._next_profiles['T'] = prof.adiabat(sl._next_profiles['r'], core._Tcen_fes + model.dt*core._dT_dt_fes, prm.fes_adiabat_params) # np.full(10, core.Tcen + model.dt*core.dT_dt)
#         sl.ds_dt = 0
#     else:
#         #Time step layer evolution
#         pure_thermal_method(model)

#     #New radial profiles
#     r_new = sl._next_profiles['r']
#     T_new = sl._next_profiles['T']


#     #Secular cooling. Compare final profiles to those at the beginning of the time-step.
#     if r_new[0] <= r_initial[0]:  #Layer has grown

#         #Need to append adiabat values onto initial temperatures
#         T_initial_temp = prof.adiabat(r_new, core.Tcen, prm.core_adiabat_params)

#         for i in range(r_new.size):
#             if r_new[i] >= r_initial[0]:
#                 T_initial_temp[i] = np.interp(r_new[i], r_initial, T_initial)
#         T_initial = T_initial_temp

#         r1, T2, T1 = r_new, T_new, T_initial

#     elif r_new[0] > r_initial[0]:  #Layer has shrunk

#         #Need to append on adiabatic values onto new temperature
#         T_new_temp = prof.adiabat(r_initial, core.Tcen + core.dT_dt*model.dt, prm.core_adiabat_params)

#         for i in range(r_initial.size):
#             if r_initial[i] >= r_new[0]:
#                 T_new_temp[i] = np.interp(r_initial[i], r_new, T_new)

#         T_new = T_new_temp

#         r1, T2, T1 = r_initial, T_new, T_initial

#     dT_dt = (T2-T1)/model.dt


#     cp  = np.interp(r1, core.profiles['r'], core.profiles['cp'])
#     rho = np.interp(r1, core.profiles['r'], core.profiles['rho'])
#     Qs += -4*np.pi*en.integrate(r1, r1**2*dT_dt*rho*cp)
#     Es += -4*np.pi*en.integrate(r1, r1**2*dT_dt*rho*cp*(1/T1[-1] - 1/T1))
#     sl.dT_dt_s = dT_dt[0]  #Save cooling rate at base of layer

#     # if core.ri > 0:
#     #     breakpoint()


#     # breakpoint()

#     Ej = Es - Ek

#     #Save energies/entropies to model attributes
#     sl.Qs, sl.Es, sl.Ek, sl.Ej = Qs, Es, Ek, Ej

#     # if sl.layer_thickness > prm.FeS_size:
#     #     breakpoint()

#     #Save secular cooling of stratified fluid for next iteration for calculating Q_fes.
#     if r_new[0] < r_fes:
#         sl._Qs_strat = -4*np.pi*en.integrate(r1[r1<=r_fes], r1[r1<=r_fes]**2 * dT_dt[r1<=r_fes] * rho[r1<=r_fes] * cp[r1<=r_fes]) #Secular cooling of stratified fluid not in FeS layer

#     #Add onto core results
#     if prm.core:
#         core.Qs += Qs
#         core.Ej += Ej
#         core.Es += Es
#         core.Ek += Ek

#         # if r_new[0] <= r_fes and ADR_fes >=1: #Convecting FeS layer secular cooling

#         #     profiles = core.profiles
#         #     r, T = profiles['r'], profiles['T']

#         #     idx = core._fes_idx

#         #     #Temp gradient at r_fes
#         #     dT_dr_fes = (T[idx-1]-T[idx-2])/(r[idx-1]-r[idx-2])

#         #     #Conductive heat flow from bulk into FeS layer.
#         #     core.Q_fes = -profiles['k'][idx-1]*dT_dr_fes*4*np.pi*r_fes**2

#         #     Qs_fes = core.Q_cmb-core.Q_fes

#         #     Qs_tilde, Es_tilde = en.secular_cool(profiles['r'][idx:], profiles['rho'][idx:], profiles['T'][idx:], profiles['cp'][idx:], prm.n_profiles-idx)
#         #     Qs_tilde = Qs_tilde * profiles['T'][idx]/core._Tcen_fes
#         #     Es_tilde = Es_tilde * profiles['T'][idx]/core._Tcen_fes #Renormalise from first T value to theoretical FeS adiabat at Tcen
#         #     dT_dt_fes = Qs_fes/Qs_tilde

#         #     #Add on energy and entropy for FeS layer
#         #     core.Qs += Qs_fes
#         #     core.Es += Es_tilde*dT_dt_fes
#         #     core._dT_dt_fes = dT_dt_fes

#     #Make sure profiles have size maximum_resolution, so they can be saved.
#     n = sl.profiles['r'].size

#     if n > prm.max_resolution:

#         r = linspace(sl.profiles['r'][0], sl.profiles['r'][-1], prm.max_resolution)
#         T = np.interp(r, sl.profiles['r'], sl.profiles['T'])

#     elif n < prm.max_resolution:

#         r = np.append(sl.profiles['r'], np.ones(prm.max_resolution-n)*sl.profiles['r'][-1])
#         T = np.append(sl.profiles['T'], np.ones(prm.max_resolution-n)*sl.profiles['T'][-1])

#     else:
#         r, T = sl.profiles['r'], sl.profiles['T']


#     sl.profiles['r'], sl.profiles['T'] = r, T

#     sl.profiles['Ta'] = prof.adiabat(r, core.Tcen, prm.core_adiabat_params)

#     sl.T_cmb, sl.T_s= T[-1], T[0]

#     core.T_cmb,  core.T_s = sl.T_cmb, sl.T_s

#     if prm.FeS_size == 0:
#         core._Tcen_fes = core.Tcen
#         core._dT_dt_fes = core._dT_dt_fes



# def update(model):
#     '''Update function

#     Parameters
#     ----------
#     model : ThermalModel
#         Main model class
#     '''

#     sl = model.stable_layer
#     prm = model.parameters

#     #Update profiles
#     sl.profiles = sl._next_profiles.copy()
#     sl.profiles['Ta'] = prof.adiabat(sl.profiles['r'], model.core.Tcen, prm.core_adiabat_params)

#     #Update core profiles if no core model is being used.
#     if not prm.core:
#         prof.basic_profiles(model)

#     #Update layer size
#     sl.layer_thickness -= sl.ds_dt * model.dt
#     if sl.layer_thickness < 0:
#         sl.layer_thickness = 0

#     model.core.rs = prm.r_cmb - sl.layer_thickness
#     model.core.T_cmb = sl.profiles['T'][-1]

#     #If layer covers entire core, Tcen is temp at base of layer.
#     if sl.layer_thickness == prm.r_cmb:
#         model.core.Tcen = sl.profiles['T'][0]
#         model.core.dT_dt = sl.dT_dt_s

#     #Track theoretical Tcen for fes layer if it is part of the conduction solution.
#     if sl.ADR_fes < 1:
#         model.core._Tcen_fes = model.core.T_cmb / prof.adiabat(prm.r_cmb, 1, prm.fes_adiabat_params)






# def progress(model):
#     '''text to display during calculation

#     Parameters
#     ----------
#     model : ThermalModel
#         Main model class

#     Returns
#     -------
#     str
#         text to display to STDOUT
#     '''

#     sl = model.stable_layer

#     v = (sl.layer_thickness/1000, sl.Q_cmb/1e12, sl.ADR, sl.T_cmb)

#     text = f'    layer thickness: {v[0]:.2f} km    Q_cmb: {v[1]:.2f} TW    ADR(rc): {v[2]:.2f}    T_cmb: {v[3]:.2f} K'

#     return text





# def pure_thermal_method(model):
#     '''Pure thermal stratification

#     Solves the thermal diffusion solution and updates the layer size during one time iteration.
#     For stability, it may step in smaller time increments until the total model timestep has been
#     reached.

#     Parameters
#     ----------
#     model : ThermalModel
#         Main model class

#     '''

#     core = model.core
#     sl  = model.stable_layer
#     prm = model.parameters

#     #Read in values from model
#     layer_thickness = sl.layer_thickness
#     r_s = prm.r_cmb - layer_thickness
#     r_s_original = r_s*1

#     r_cmb      = prm.r_cmb
#     resolution = prm.resolution
#     min_res    = prm.min_resolution
#     max_res    = prm.max_resolution

#     TOL = prm.depth_tolerance   #Minimum layer size
#     init_size = prm.init_size + prm.FeS_size

#     core_adiabat      = prm.core_adiabat_params
#     rho_poly          = prm.core_liquid_density_params

#     E_T = prm.entrainment_T

#     #Heat flows
#     Q_cmb  = model.mantle.Q_cmb
#     Q_fes = core.Q_fes

#     Tcen   = model.core.Tcen
#     dTa_dt = model.core.dT_dt

#     #Base of FeS layer
#     r_fes = prm.r_cmb - prm.FeS_size

#     #Existing radial profiles
#     r = sl.profiles['r']
#     T = sl.profiles['T']

#     #Set upper boundary conditions
#     #If FeS layer is sub-adiabatic, solution as normal
#     ub = prof.adiabat_grad(prm.r_cmb, core._Tcen_fes, prm.fes_adiabat_params)*sl.ADR_fes
#     ub = core.Q_cmb/(4*np.pi*prm.r_cmb**2*-core.profiles['k'][-1])
#     ub_type = 1

#     # if sl.ADR_fes >= 1:
#     #     sl.ADR_fes = 0.99

#     if sl.ADR_fes >= 1:
#         r_upper = r_fes

#     r_upper = r_fes #Force convecting FeS layer for now.


#     #For stability, reduce the time step size when layer is relatively thin.
#     factor = 1
#     # if layer_thickness <= init_size*3:
#     #     factor = 0.01

#     time_gone = 0
#     while time_gone < model.dt:

#         coarse = prm.coarse_resolution
#         fine   = prm.fine_resolution

#         #Initialise layer with minium tolerance thickness
#         if layer_thickness < TOL + prm.FeS_size:

#             layer_thickness = init_size

#             #Number of grid points in this iteration
#             n_points = int(np.max([min_res,layer_thickness*resolution]))
#             n_points = int(np.min([max_res,n_points]))

#             r_s = prm.r_cmb - layer_thickness

#             # r_test = np.linspace(r_s, r_upper, n_points)
            

#             r = adaptive_grid(r_s, r_fes, r_upper, coarse, fine, transition_width=200e3)

#             T = prof.adiabat(r, Tcen, core_adiabat)
#             # T_test = prof.adiabat(r_test, Tcen, core_adiabat)

#         elif r_upper == r_fes:
#             #FeS layer is still convecting so just solve for bulk core that is thermally stratified.

#             #Number of grid points in this iteration
#             n_points = int(np.max([min_res,int((r_upper-r_s)/coarse)]))
#             n_points = int(np.min([max_res,n_points]))

#             r_s = prm.r_cmb - layer_thickness

#             # coarse = np.min([10e3, (r_upper-r_s)/10])
#             # fine = np.min([2e3, (r_upper-r_s)/20])
#             r = np.linspace(r_s, r_upper, n_points)
#             # r = adaptive_grid(r_s, r_fes, r_upper, coarse, fine, transition_width=200e3)

#             T = np.interp(r, sl.profiles['r'], sl.profiles['T'])

#         else:
#             #Make sure solution is on even grid
#             n_points = int(np.max([min_res,layer_thickness*resolution]))
#             n_points = int(np.min([max_res,n_points]))

#             # coarse = np.min([10e3, (r_upper-r_s)/10])
#             # fine = np.min([2e3, (r_upper-r_s)/20])

#             r = adaptive_grid(r_s, r_fes, r_upper, coarse, fine, transition_width=200e3)
#             # r = np.linspace(sl.profiles['r'][0], sl.profiles['r'][-1], n_points)
#             T = np.interp(r, sl.profiles['r'], sl.profiles['T'])

#             # breakpoint()


#         dt_small = model.dt*factor

#         if time_gone + dt_small > model.dt:
#             dt_small = model.dt - time_gone

#         #Regrid core radial profiles onto finer stable layer grid
#         rho = prof.density(r, rho_poly)

#         #n = int((r.size-1)*(layer_thickness-prm.FeS_layer_size)/layer_thickness)

#         # FeS_idx = np.where(r >= prm.r_cmb-prm.FeS_layer_size)[0][0]
#         # if FeS_idx >= r.size-2:
#         #     breakpoint()
#         # k, dk_dr = tanh_conductivity(r, core.profiles['k'][0], prm.FeS_conductivity, FeS_idx)


#         if r_upper == prm.r_cmb:
#             k, dk_dr = tanh_conductivity(r, core.profiles['k'][0], prm.core_fes_conductivity, r_fes, transition_width=5000)
#         else:
#             k = core.profiles['k']
#             dk_dr = np.zeros(k.size)

#             k   = np.interp(r, core.profiles['r'], k)
#             dk_dr = np.interp(r, core.profiles['r'], dk_dr)

#             # k_test = np.full(r_test.size, k[0])
#             # dk_dr_test = np.zeros(r_test.size)



#         # if layer_thickness > prm.FeS_layer_size:
#         #     breakpoint()

#         cp    = np.interp(r, core.profiles['r'], core.profiles['cp'])
#         alpha = np.interp(r, core.profiles['r'], core.profiles['alpha'])

#         #Thermal diffusivity
#         D_T = k/(rho*cp)
#         #dk_dr * 1/(rho*cp). Not the same as the diffusivity gradient

#         # D_T_test = np.full(r_test.size, D_T[0])

#         # dD_dr = np.gradient(k, r, edge_order=2)/(rho*cp)
#         dD_dr = dk_dr/(rho*cp)


#         #Cool adiabat
#         Tcen = core.Tcen + dTa_dt*(time_gone+dt_small)
#         Tcen_fes = core._Tcen_fes + core._dT_dt_fes*(time_gone+dt_small)

#         Ta, Ta_grad = prof.adiabat(r, Tcen, core_adiabat), prof.adiabat_grad(r, Tcen, core_adiabat)

#         Ta[r>r_fes] = prof.adiabat(r[r>r_fes], Tcen_fes, prm.fes_adiabat_params)


#         #Lower boundary condition (using entrainment)
#         lb = (1-E_T)*Ta_grad[0]
#         lb_type = 1

#         #Fes layer is convecting so change upper BC.
#         if r_upper == r_fes:
#             #Fixed temperature boundary condition
#             ub      = prof.adiabat(r_fes, Tcen_fes, prm.fes_adiabat_params)
#             ub_type = 0

#         #If whole core is stratified, lower BC should be zero flux
#         if r[0] == 0:
#             lb = 0
   

#         #Calculate diffusion solution
#         T_new = diffusion2(T, r, dt_small, D_T, k, dk_dr, (lb_type,lb),(ub_type,ub))
#         # T_new_test = diffusion(T_test, r_test, dt_small, D_T_test, k_test, dk_dr_test, (lb_type,lb),(ub_type,ub))

#         Qs = -4*np.pi*trapezoid(r, r**2 * rho* cp *(T_new-T)/dt_small)[-1]
#         # print(f'Qs: {Qs: .3e}   Q_fes: {core.Q_fes-core.Qs: .3e}')

#         # print(model.it, t2-t1, np.min(np.diff(r)))
#         # if model.it > 75 and ub_type == 1:

#         #     Q_top = 4*np.pi*r_cmb**2 * k[-1] * -ub
#         #     Q_bottom = 4*np.pi*rs**2 * k[0] * -lb

#         #     Qs = 
        
#         # if core.ri>0:
#         #     import matplotlib.pyplot as plt
#         #     plt.plot(r, T_new, label='T_new')
#         #     plt.plot(r, T, label='T')
#         #     plt.plot(r, Ta, label='Ta')
#         #     plt.legend(loc=0)
#         #     plt.savefig('test.png')
#         #     breakpoint()



#         #Set adiabat to track diffusion solution
#         if r[0] == 0:
#             Tcen = T_new[0]
#             Ta, Ta_grad = prof.adiabat(r, Tcen-100, core_adiabat), prof.adiabat_grad(r, Tcen, core_adiabat)


#         if np.isnan(T_new[0]):
#             breakpoint()


#         #Mix layer (experimental function, not used by default)
#         if prm.mix_layer:
#             T_rel = func.mix_profile(r, cp*rho*r**2, T_new-Ta)
#             T_new = T_rel + Ta

#         if sl.ADR < 1 and r[0]==0:
#             r_s_new = 0
#         elif r[0] == 0 and model.it == 1:
#             r_s_new = 0
#         else:
#             #Determine if layer should retreat
#             r_s_new = func.retreat_layer(r[r<=r_fes], T_new[r<=r_fes], np.zeros(r.size)[r<=r_fes], Ta[r<=r_fes], 0, alpha[r<=r_fes], 0, density_grad_limit=1e-6)

#         # print(f'it: {model.it}   Qs: {core.Qs/1e12: .3f}   Qa: {core.Qa_rs/1e12: .3f}   size: {sl.layer_thickness/1000}    dT_dt_fes:{core._dT_dt_fes: .3e}')

#         #If whole core is stratified, only ADR>=1 can change the layer thickness.
#         if r[0] == 0 and sl.ADR < 1:
#             r_s_new = 0

#         #Grow layer
#         elif r_s_new == r[0]:

#             #Find radius on adiabat that matches temperature at base of layer
#             def f(guess,T,Tcen):
#                 return T - prof.adiabat(guess,Tcen,core_adiabat)

#             if T_new[0] >= Tcen: #Layer has reached center of core
#                 r_s_new = 0
#             else:
#                 r_s_new = bisect(f,0,r_upper,args=(T_new[0],Tcen),maxiter=200)

#         #Shrink layer
#         else:
#             # import matplotlib.pyplot as plt
#             # plt.plot(r, T_new, label='T_new')
#             # plt.plot(r, T, label='T')
#             # plt.plot(r, Ta, label='Ta')
#             # plt.legend(loc=0)
#             # plt.show()
#             # breakpoint()

#             # # plt.savefig('test.png')
#             # breakpoint()

#             if r_s_new >= r_upper-TOL:

#                 #No stratification. Try a larger timestep.
#                 r_s_new = r_upper
#                 factor = factor*2

#                 if factor > 64 or sl.ADR >= 1: #Give up is large time step doesn't help or core is super-adiabatic.
#                     dt_small = model.dt - time_gone
#                     factor = 1

#         #r_s must be include FeS layer.
#         # r_s_new = np.min([r_s_new, prm.r_cmb-prm.FeS_layer_size])
        
#         #If FeS layer is sub-adiabatic, r_s is at least at base of FeS layer.
#         if sl.ADR_fes < 1:
#             r_s_new = np.min([r_fes, r_s_new])


#         layer_thickness = r_upper - r_s_new

#         #Regrid domain, keeping temperature relative to the adiabat.
#         T_rel = T_new - prof.adiabat(r, Tcen, core_adiabat)

#         r, T_rel = func.change_domain_size(T_rel, r, r_s_new, (resolution,min_res,max_res))

#         T = T_rel + prof.adiabat(r, Tcen, core_adiabat)

#         time_gone += dt_small

#     if r_upper < r_cmb:
#         #Append values in the FeS layer
#         r = np.append(r, np.linspace(r_upper, prm.r_cmb,10)[1:])
#         T = np.append(T, prof.adiabat(r[-9:], Tcen_fes, prm.fes_adiabat_params))

#     #No layer this iteration
#     if r_s_new == r_upper:
#         layer_thickness = prm.FeS_size
#         r_s_new = r_fes
#         r, T = r[:9], T[:9]

#     #Save new profiles. Keep original profiles until update() has been called.
#     sl._next_profiles['r'] = r
#     sl._next_profiles['T'] = T

#     sl.ds_dt = (r_s_new - r_s_original)/model.dt

#     if r[1]-r[0] > 0:
#         sl.T_grad_s = (T[1]-T[0])/(r[1]-r[0])
#     else:
#         sl.T_grad_s = 0

#     sl.lb_T = lb
#     sl.ub_T = ub



#Overwrite existing conductivity function with this one to set a conductivity anomaly at the top of the core.
#Fes layer size must be included in conductivity params.

# def conductivity_FeS(r, P, T, conductivity_params):
    
#     if len(conductivity_params)==1:
#         k = np.ones(r.size)*conductivity_params[0]
    
#     else:

#         k = np.ones(r.size)*conductivity_params[0]

#         k[r>conductivity_params[1]] = conductivity_params[2]

#     return k

# prof.conductivity = conductivity_FeS