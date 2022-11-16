#Radial Profiles Functions
import numpy as np
from thermal_history.utils.optimised_funcs import trapezoid, polyval

#Expose the rest of profiles functions to the ones here
from ..leeds.routines.profiles import *

#Import leeds to overwrite functions we define here
from .. import leeds

import logging
logger = logging.getLogger(__name__)


#Redefine basic_profiles to include FeS layer
def basic_profiles(model, setup=False):
    '''Main function setting radial profiles for core properties.

    Sets radius, rho, g, P, cp, alpha_T, psi, Ta, dTa_dr profiles.
    Also will set the adiabat polynomials if required.

    Parameters
    ----------
    model : ThermalModel class
        Main class being used for the calculation
    setup : bool, optional
        Flag to setup the profiles dictionary during model setup, by default False
    '''

    prm = model.parameters
    n = prm.n_profiles

    core = model.core

    if setup:
        profiles = {}
        for key in ['r', 'rho', 'g', 'P', 'cp', 'alpha', 'psi', 'Ta', 'dTa_dr']:
            profiles[key] = np.zeros(n)
    else:
        profiles = core.profiles
        #Get existing pressure at new ri to keep new pressure field with this value.
        #This ensures pressure field is not changing through time at ri, avoiding
        #difference between ri(r+dt) and ri(t)+Ct*dT_dt. This could be properly
        #accounted for as in Gubbins et al. (2003)
        # P_ri = np.interp(core.ri, profiles['r'], profiles['P'])

    #Radial array
    r_cmb = prm.r_cmb
    ri    = core.ri
    r_snow= core.r_snow
    rs    = core.rs
    r_fes = r_cmb-prm.FeS_size #!# Include r_fes

    #Just want the unique radii in order
    radii = np.unique(np.sort([0, ri, r_snow, rs, r_fes, r_cmb]))  #!# Include r_fes

    #Set radial grid
    r = radial_grid(radii, n)
    #Add to profiles after pressure has been calculated

    #Find indexes for the start of each layer
    for i in range(n):
        if r[i] == ri:
            core._ri_idx = i
        if r[i] == rs:
            core._rs_idx = i
        if r[i] == r_snow:
            core._snow_idx = i
        if r[i] == r_fes:       #!# Include r_fes
            core._fes_idx = i

    ri_idx = core._ri_idx
    fes_idx = core._fes_idx     #!# Include r_fes

    #Density
    rho = density(r, prm.core_liquid_density_params)
    if ri > 0:
        rho[:ri_idx] = density(r[:ri_idx], prm.core_solid_density_params)

    rho[fes_idx:] = prm.FeS_density #!# Include r_fes

    profiles['rho'] = rho
    ##########################


    #Gravity
    g = gravity(r, prm.core_liquid_density_params)
    if ri > 0:
        g_ri = g[ri_idx]
        g[:ri_idx] = gravity(r[:ri_idx], prm.core_solid_density_params)
        g[ri_idx:] = g[ri_idx:] + (g[ri_idx-1]-g_ri)*(ri**2/r[ri_idx:]**2)

    profiles['g'] = g
    ##########################


    #Gravitational Potential
    psi = grav_potential(r, prm.core_liquid_density_params)
    if ri > 0:
        psi_ri = psi[ri_idx-1]
        psi[ri_idx:] = grav_potential(r[ri_idx:], prm.core_liquid_density_params) + (psi[ri_idx-1]-psi_ri)*r[ri_idx:]/ri

    psi = psi - psi[-1] #shift reference to psi(cmb)=0
    profiles['psi'] = psi
    ##########################


    #Pressure. Doesn't update with each iteration.
    if setup:
        P = trapezoid(r[::-1], -rho[::-1]*g[::-1])[::-1] + prm.P_cmb  #Integrate top down and swap order back
        profiles['P'] = P  #swap order back
    else:
        P = np.interp(r, profiles['r'], profiles['P']) #Pressure field remains the same to keep dri_dt consisent with intercept of Tm/T
        # P = trapezoid(r[::-1], -rho[::-1]*g[::-1]) + prm.P_cmb  #Integrate top down
        # P = P[::-1] #swap order back
        # profiles['P'] = P

        # P += P_ri-P[ri_idx] #Keep Pressure field constant in time (at least in the viscinity of ri)


    profiles['P'] = P
    profiles['r'] = r

    ##########################


    #Specific heat
    profiles['cp'] = specific_heat(r, P, prm.core_cp_params)
    profiles['cp'][fes_idx:] = prm.FeS_cp                      #!# Include FeS
    ##########################

    #Thermal expansion
    profiles['alpha'] = thermal_expansivity(r, P, prm.core_alpha_T_params)
    ##########################

    #Adiabatic temperature profile
    if setup:

        if prm.core_adiabat_params is None:
            prm.core_adiabat_params = fit_adiabat_polynomials(r, g, profiles['cp'], profiles['alpha'])

        core.Tcen = core.T_cmb/polyval(prm.core_adiabat_params[::-1], prm.r_cmb)

        if not prm.Tcen is None:
            core.T_cmb = prm.T_cen*polyval(prm.core_adiabat_params[::-1], prm.r_cmb)
            core.Tcen = prm.T_cen
            logger.warning('Tcen has been specified, takes priority over T_cmb')


    Ta = adiabat(r, core.Tcen, prm.core_adiabat_params)
    profiles['Ta'] = Ta

    if setup:
        #Set initial temperature profile to adiabat
        profiles['T'] = Ta
    ##########################

    #Adiabatic gradient
    dTa_dr = adiabat_grad(r, core.Tcen, prm.core_adiabat_params)
    profiles['dTa_dr'] = dTa_dr
    ##########################

    core.profiles = profiles

leeds.routines.profiles.basic_profiles = basic_profiles


#Redefine basic_profiles to include FeS layer conductivity
def temp_dependent_profiles(model, setup=False):
    '''Set the temperature dependent profiles.

    Sets the k, Qa, T profiles.

    Parameters
    ----------
    model : ThermalModel class
        Main class being used for the calculation
    setup : bool, optional
        Flag to setup the profiles dictionary during model setup, by default False
    '''

    core = model.core
    sl   = model.stable_layer
    prm  = model.parameters

    n = prm.n_profiles

    profiles = core.profiles

    if setup:
        for key in ['k', 'Qa', 'T']:
            profiles[key] = np.zeros(n)


    r = profiles['r']
    P = profiles['P']

    #Temperature profile
    T = profiles['Ta'].copy()

    #If a stable layer is present, interpolate it's temperature profile onto the upper core.
    if (not setup) and (prm.stable_layer and core.rs < prm.r_cmb):
        T[core._rs_idx:] = np.interp(r[core._rs_idx:], sl.profiles['r'], sl.profiles['T'])

    profiles['T'] = T


    #Conductivity (may depend on temperature depending on user choise of params.)
    k = conductivity(r, P, T, prm.core_conductivity_params)
    core.profiles['k'] = k
    core.profiles['k'][core._fes_idx:] = prm.FeS_conductivity
    ##########################

    #Adiabatic heat flow. Depends on k and hence may depend on T.
    core.profiles['Qa'] = adiabatic_heat_flow(r, k, core.Tcen,  prm.core_adiabat_params)

leeds.routines.profiles.temp_dependent_profiles = temp_dependent_profiles




