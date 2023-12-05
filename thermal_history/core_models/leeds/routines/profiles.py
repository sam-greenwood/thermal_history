#Radial Profiles Functions
import numpy as np
from scipy.optimize import curve_fit
from numba import njit
from thermal_history.utils.optimised_funcs import linspace, trapezoid, polyval

import logging
logger = logging.getLogger(__name__)

G = 6.67e-11 #Gravitational constant

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

    #Just want the unique radii in order
    radii = np.unique(np.sort([0, ri, r_snow, rs, r_cmb]))

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

    ri_idx = core._ri_idx

    #Density
    rho = density(r, prm.core_liquid_density_params)
    if ri > 0:
        rho[:ri_idx] = density(r[:ri_idx], prm.core_solid_density_params)
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
    
    #CD - added in AP conductivity
    if prm.core_conductivity_params[0] == 'AP':
        k = prm.core_conductivity_params[1].astype('float')*core.conc_l[0]**2 + \
        prm.core_conductivity_params[2].astype('float')*core.conc_l[0]    + \
        prm.core_conductivity_params[3].astype('float')*P**3  + \
        prm.core_conductivity_params[4].astype('float')*P**2  + \
        prm.core_conductivity_params[5].astype('float')*P     + \
        prm.core_conductivity_params[6].astype('float')*T     + \
        prm.core_conductivity_params[7].astype('float')
            
    core.profiles['k'] = k
    ##########################

    #Adiabatic heat flow. Depends on k and hence may depend on T.
    core.profiles['Qa'] = adiabatic_heat_flow(r, k, core.Tcen,  prm.core_adiabat_params)


# @njit
def radial_grid(radii, n_points):
    '''Create the radial grid across the core.

    Creates the adaptive radial grid that ensures 2 radii at each interface to allow
    for resolving dicontinuous variables. E.g. 2 points at the ICB to resolve density jump.
    

    Parameters
    ----------
    radii : array
        radii for layer interfaces throughout the core, in order.
    n_points : int 
        Number of radial grid points

    Returns
    -------
    array
        radial grid
    '''

    sizes = np.zeros(radii.size-1, dtype=int)

    #Calculate proportional number of radial points for each region (need at least 2 for lower and upper radii)
    #These must also be int's as they are used for indexing
    sizes[0] = np.max([int(n_points*radii[1]/radii[-1]), 2])

    for i in range(1,sizes.size):
        sizes[i] = np.max([int(n_points*(radii[i+1]-radii[i])/radii[-1]), 2])

    #The rounding to a whole number of grid points above means there may be too many/few
    #grid points. Take off/add 1 grid point from/to each region progressively until
    #there are the correct total of n_points.
    dn = np.sum(sizes) - n_points

    i = 0
    while not dn == 0:
        if dn>0:
            delta = -1 #points must be taken off
        else:
            delta = 1  #points must be added

        if sizes[i] > 2:
            sizes[i] += delta

        dn = np.sum(sizes) - n_points

        if i == len(sizes)-1: #If we're at the end of sizes, restart again until dn=0
            i = 0
        else:
            i+=1  #Move to the next layer

    #Create radius array
    r = np.zeros(n_points)

    #Add linear spacing between each radii
    start = 0
    for i in range(sizes.size):
        end = start+sizes[i]
        r[start:end] = linspace(radii[i], radii[i+1], sizes[i])
        start=end

    return r


def density(r, density_params):
    '''Density in the core

    Parameters
    ----------
    r : float or array
        radius
    density_params : array
        radial polynomial coefficients for density

    Returns
    -------
    float or array
        density at r
    '''

    rho = polyval(density_params[::-1], r)

    return rho


def gravity(r, density_params):
    '''Gravity in the core

    Parameters
    ----------
    r : float or array
        radius
    density_params : array
        radial polynomial coefficients for density

    Returns
    -------
    float or array
        gravity at r
    '''
    poly = np.zeros(density_params.size+1)

    for i in range(density_params.size):
        poly[i+1] = density_params[i]/(i+3)

    g = (4*np.pi*G)*polyval(poly[::-1],r)

    return g


def grav_potential(r, density_params):
    '''Gravitational potential in the core

    Parameters
    ----------
    r : float or array
        radius
    density_params : array
        radial polynomial coefficients for density

    Returns
    -------
    float or array
        gravititational potential at r
    '''

    #gravity polynomials
    poly = density_params/np.arange(3,density_params.size+3)
    poly = np.concatenate(([0],poly))

    #integrate to gravitational potential
    poly = poly/np.arange(1,poly.size+1)
    poly = np.concatenate(([0],poly))

    psi = (4*np.pi*G)*polyval(poly[::-1],r)

    return psi


def specific_heat(r, P, cp_params):
    '''Specific heat capacity as a function of Pressure.

    Parameters
    ----------
    r : array
        radius
    P : array
        Pressure
    cp_params : array
        pressure polynomial coefficients for specific heat

    Returns
    -------
    array
        specific heat capacity at P

    Raises
    ------
    ValueError
        If more than constant cp is given (cp_params is longer than one number), not yet implemented
    '''


    cp = np.zeros(r.size)
    if len(cp_params)==1:
        cp[:] = cp_params[0]
    else:
       raise ValueError('Not yet implimented! Include code here to calculate pressure dependence of cp')
    return cp


def thermal_expansivity(r, P,  alpha_params):
    '''Thermal expansivity as a function of Pressure.

    Parameterisation of thermal expansivity given in Knibbe and Van Westrenen (2018).
    If only a single value is given in the alpha_params array, then alpha will be constant
    with this value.

    Parameters
    ----------
    r : array
        radius
    P : array
        pressure
    alpha_params : array 
        parameters for thermal expansivity in the order: alpha_T at ambient pressure, bulk modulus, derivative of bulk modulus with pressure

    Returns
    -------
    array
        thermal expansivity at P
    '''

    alpha = np.zeros(r.size)
    if len(alpha_params)==1:
        alpha[:] = alpha_params[0]
    else:
        alpha_0 = alpha_params[0]
        KT = alpha_params[1]
        dK_dP = alpha_params[2]
        for i in range(r.size):
            alpha[i] = alpha_0*KT/(KT + P[i]*dK_dP)

    return alpha


def fit_adiabat_polynomials(r, g, cp, alpha):
    '''Calculate radial polynomial coefficients for adiabatic temperature

    Numerically integrates the adiabatic gradient to get the temperature based on
    given gravity, specific heat and thermal expansivity arrays. Then fits 3rd order
    polynomial coefficients normalised to the temperature at r[0]. The linear temperature
    dependence is also forced to be zero such that the adiabatic gradient is zero at r=0.

    Parameters
    ----------
    r : array
        radius
    g : float or array
        gravity
    cp : float or array
        specific heat
    alpha : float or array
        thermal expansivity

    Returns
    -------
    array
        adiabatic temperature polynomial coefficients
    '''

    if np.array(g).size == 1:
        g = np.full(r.size, g)
    if np.array(cp).size == 1:
        cp = np.full(r.size, cp)
    if np.array(alpha).size == 1:
        alpha = np.full(r.size, alpha)


    #Numerically integrate adiabatic gradient for T[0]=1
    T = np.ones(r.size)
    for i in range(1,r.size):
        dr = r[i]-r[i-1]
        try:
            T[i] = T[i-1]*(1 - dr*g[i-1]*alpha[i-1]/cp[i-1])
        except:
            breakpoint()

    #We require zero adiabatic gradient at r=0, so fit radial polynomials
    #such that poly = [1,0,x1,x2]
    def f(r,*args):
        poly = np.insert(args,0,[1,0])
        return polyval(poly[::-1],r)

    #Solve for x1,x2 and then create full poly array
    poly, _ = curve_fit(f, r, T, p0=[0,0])
    poly = np.insert(poly, 0, [1,0])

    return poly


def adiabat(r, Tcen, adiabat_params):
    '''Adiabatic temperature

    Parameters
    ----------
    r : array
        radius
    Tcen : float
        Temperature at center of core
    adiabat_params : array
        adiabat polynomial coefficients

    Returns
    -------
    array
        Adiabatic temperature at r
    '''

    Ta = Tcen*polyval(adiabat_params[::-1], r)

    return Ta


def adiabat_grad(r, Tcen, adiabat_params):
    '''Adiabatic temperature gradient

    Parameters
    ----------
    r : array
        radius
    Tcen : float
        Temperature at center of core
    adiabat_params : array
        adiabat polynomial coefficients

    Returns
    -------
    array
        Adiabatic temperature gradient at r
    '''

    poly = np.zeros(adiabat_params.size-1)
    for i in range(poly.size):
        poly[i] = adiabat_params[i+1] * (i+1)

    Ta = Tcen*polyval(poly[::-1], r)

    return Ta


def conductivity(r, P, T, conductivity_params):
    '''Thermal conductivity

    Calculates the thermal conductivity that, depending on the supplied parameters,
    can be a function of radius, temperature, or pressure. The first element of the
    conductivity parameters should be a string ('r'/'T'/'P') indicating the type of
    polynomial coefficients. If only a single value is give, the thermal conductivity
    is given as constant with that value.

    Parameters
    ----------
    r : array
        radius
    P : arrat
        pressure
    T : array
        temperature
    conductivity_params : array
        Conducitvity polynomial coefficients

    Returns
    -------
    array
        thermal conductivity at r/P/T

    Raises
    ------
    ValueError
        If more than one value is given in conductivity_params but no string is given indicating type of polynomial coefficients.
    '''

    if len(conductivity_params)==1:
        k = np.ones(r.size)*conductivity_params[0]

    elif conductivity_params[0] == 'T':
        k = polyval(conductivity_params[1:].astype('float')[::-1], T)

    elif conductivity_params[0] == 'P':
        k = polyval(conductivity_params[1:].astype('float')[::-1], P)

    elif conductivity_params[0] == 'r':
        k = polyval(conductivity_params[1:].astype('float')[::-1], r)
        
    elif conductivity_params[0] == 'AP':

        LE = np.array([0.1]) #Set to fixed value for now just as an example
        
        "Set outside routine as k depends on P, T and c"
        k = conductivity_params[1].astype('float')*LE[0]**2 + \
            conductivity_params[2].astype('float')*LE[0]    + \
            conductivity_params[3].astype('float')*P**3  + \
            conductivity_params[4].astype('float')*P**2  + \
            conductivity_params[5].astype('float')*P     + \
            conductivity_params[6].astype('float')*T     + \
            conductivity_params[7].astype('float')
        
    else:
        raise ValueError('Unless only one value is given in core_conductivity_params, first value must be a string denoting if polynomials are in r/T/P')

    return k


def adiabatic_heat_flow(r,k,Tcen,adiabat_poly):
    '''Adiabatic heat flow

    Parameters
    ----------
    r : array
        radius
    k : array
        thermal conductivity
    Tcen : float
        Temperature at center of core
    adiabat_poly : array
        adiabat polynomial coefficients

    Returns
    -------
    array
        Adiabatic heat flow at r
    '''

    Qa = 4*np.pi*r**2 * k * -adiabat_grad(r,Tcen,adiabat_poly)

    return Qa


def entropy_melting(P, ent_melt_params):
    '''Change in entropy on fusion

    Parameters
    ----------
    P : float or array
        pressure
    ent_melt_params : array
        pressure polynomial coefficients for change in entropy on fusion

    Returns
    -------
    array
        Change in entropy on fusion
    '''

    ds_fe = polyval(ent_melt_params[::-1],P)

    return ds_fe


def mass(rho_ic_poly, rho_oc_poly, ri, rs, r_cmb):
    '''mass of regions within the core

    Parameters
    ----------
    rho_ic_poly : array
        inner core density polynomial coefficients
    rho_oc_poly : array
        outer core density polynomial coefficient
    ri : float
        inner core radius
    rs : float
        radius for base of stable layer
    r_cmb : float
        CMB radius

    Returns
    -------
    (float, float, float, float)
        Masses for the inner core, convecting outer core, stable layer, total mass
    '''

    #Integrate density polynomials to get mass polynomials
    poly_ic = np.zeros(rho_ic_poly.size + 3)
    poly_oc = np.zeros(rho_oc_poly.size + 3)

    for i in range(3, rho_ic_poly.size+3):
        poly_ic[i] = rho_ic_poly[i-3]/i

    for i in range(3, rho_oc_poly.size+3):
        poly_oc[i] = rho_oc_poly[i-3]/i


    M_ic   = 4*np.pi*polyval(poly_ic[::-1],ri)
    M_conv = 4*np.pi*(polyval(poly_oc[::-1],rs) - polyval(poly_oc[::-1],ri))
    Ms     = 4*np.pi*(polyval(poly_oc[::-1],r_cmb) - polyval(poly_oc[::-1],rs))

    M = M_ic + M_conv + Ms

    return M_ic, M_conv, Ms, M


def mass_correct(solid_params, liquid_params, ri, M, r_cmb):
    '''Correct outer core density polynomials to account for changing LE concentration.

    Parameters
    ----------
    solid_params : array
        inner core density polynomial coefficients
    liquid_params : array
        outer core density polynomial coefficients
    ri : float
        inner core radius
    M : float
        mass of the core
    r_cmb : float
        CMB radius

    Returns
    -------
    array
        corrected outer core density polynomial coefficients
    '''

    M_ic, M_conv, Ms = mass(solid_params, liquid_params, ri, r_cmb, r_cmb)[:3]
    M_oc = M_conv+Ms

    M_oc_target = M - M_ic

    dM = M_oc_target - M_oc

    drho = 3*dM/(4*np.pi*(r_cmb**3-ri**3))

    liquid_params[0] = liquid_params[0] + drho



    return liquid_params


###############################################################################
#Calculate inner core radius
###############################################################################
# import matplotlib.pyplot as plt
#@jit
def ic_radius(r,Ta,Tm):
    '''Find the inner core radius by intersection of melting curve

    Parameters
    ----------
    r : array
        radius
    Ta : array
        temperature
    Tm : array
        melting temperature

    Returns
    -------
    float
        inner core radius

    Raises
    ------
    ValueError
        If ri comes to be NaN, r, T, Tm arrays will be logged.
    '''

    if Ta[0] >= Tm[0]:
        ri = 0
    elif Ta[-1] <= Tm[-1]:
        ri = r[-1]
    else:

        for i in range(1, r.size):
            if Ta[i] > Tm[i]:
                dr = r[i]-r[i-1]
                m1 = (Ta[i]-Ta[i-1])/dr
                m2 = (Tm[i]-Tm[i-1])/dr
                ri = r[i-1] + (Ta[i-1]-Tm[i-1])/(m2-m1)

                break

    if np.isnan(ri):
        logger.critical(f'Inner core radius is NaN\n\n r   T   Tm\n{np.column_stack((r, Ta, Tm))}')
        raise ValueError('Inner core radius is NaN, logging r, Ta, Tm arrays')

    return ri





