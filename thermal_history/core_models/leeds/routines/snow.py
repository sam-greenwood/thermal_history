from scipy.optimize import bisect, brentq
import numpy as np
from numba import njit
from thermal_history.utils.optimised_funcs import trapezoid
from thermal_history.core_models.leeds.routines.chemistry import riv_Tm

import logging
logger = logging.getLogger(__name__)

@njit
def snow_radius(r,T,Tm):
    '''Finds the snow zone radius by intersection of the temperature with the melting curve.

    Parameters
    ----------
    r : array
        Radius
    T : array
        Temperature
    Tm : array
        Melting temperature

    Returns
    -------
    float
        Snow zone radius
    '''
    

    #Snow hasn't started yet
    if np.min(T-Tm) > 0:
        r_snow = r[-1]

    elif Tm[0] >= T[0]:
        r_snow = r[0]
    
    #Find first intersection of Tm and T
    else:
        for i in range(r.size):

            if Tm[i]>=T[i]:

                dr = r[i]-r[i-1]
                if not dr == 0:
                    m1 = (T[i]-T[i-1])/dr
                    m2 = (Tm[i]-Tm[i-1])/dr
                    r_snow = r[i-1] + (T[i-1]-Tm[i-1])/(m2-m1)
                    break

        if r_snow > r[-1]:
            r_snow = r[-1]
        elif r_snow < r[0]:
            r_snow = r[0]

    return r_snow


def check_top_down_freezing(r, T, Tm):
    '''Checks for top-down freezing regime

    Parameters
    ----------
    r : array   
        radius
    T : array
        Temperature
    Tm : array
        Melting temperature
    '''


    top_down = True
    
    if T[0]-Tm[0] >= 0:
        flag = 'liquid'

        #Still fully liquid
        if np.min(T-Tm) > 0:
            #Check if likely to be top_down
            if np.min(T-(T[-1]-Tm[-1]) - Tm) < 0:
                logger.warning('May not be top-down freezing, currently regions of Tm are steeper than T.')

                #Set top_down to False if freezing starts
                if T[0]<Tm[0]:
                    top_down=False
            
    else:
        flag = 'solid'
    
    #Iterate up and check for intermediate snow zones or bottom up
    for i in range(1,r.size):

        dT = T[i]-Tm[i]
        if dT >= 0 and flag == 'solid':
            top_down = False #Gone from region of sub-liquidus to super-liquidus

        elif dT < 0:
            flag = 'solid' #Entered base of snow_zone

    return top_down





def snow_composition(Ta,Tm_fe,initial_conc,snow_index, conc_l, P, melting_params):

    '''
    Calculates the composition of the slurry such that the melting temperature
    is elevated to the adiabatic temperature. The Williams and Nimmo (2004)
    parameterisation for the melting curve is assumed:
    Tm = Tm_fe*(1-conc_l+initial_conc)

    Parameters
    ----------
    Ta : array
        Adiabatic temperature (array)
    Tm_fe : array
        Melting temperature of pure iron (array)
    intial_conc : float
        Initial concentration of liquid
    snow_index : int
        Index from which snow zone starts in arrays
    melting_params: array
        Array of melting temperature parameters

    Returns
    -------
    array
        Radial concentration profile in snow zone
    '''

    n = snow_index

    if melting_params[0] == 'RI':

        params = melting_params[1:].astype('float64')
    
        #Difference between polynomial evalutation at guessed composition and the target temp
        def f(guess,P,Target,params):
            return riv_Tm(P,guess,params) - Target

        conc_l_snow = np.zeros(P[n:].size)

        conc_l_snow[0] = brentq(f, 0, 1, args=(P[n],Ta[n], params)) #Calculate first value at interface

        #Find composition that gives Tm=Ta.
        for i in range(n+1,P.size):
            try:
                #Try optimiser brackets close to the previous value for speed. We generally expect small changes in composition from grid point to grid point
                conc_l_snow[i-n] = brentq(f, conc_l_snow[i-n-1]-0.001, conc_l_snow[i-n-1]+0.001, args=(P[i],Ta[i], params), rtol=0.0001)
            except:
                #If that fails it's likely because f(a) and f(b) are not opposite signs. Try the full range of possible values
                conc_l_snow[i-n] = brentq(f, 0, 1, args=(P[i],Ta[i], params), rtol=0.0001)


    else:

        if not melting_params[0] == 'WN':
            logger.warning('Assuming Williams and Nimmo melting curve (WN) for snow zone.')

        conc_l_snow = (1 + initial_conc)*(1 - Ta[n:]/Tm_fe[n:])

    if np.min(conc_l_snow) < 0:
       logger.warning('Warning conc_l below zero!: {}'.format(np.min(conc_l_snow)))

    return conc_l_snow

@njit
def Cl_calc(phi,L,T,conc_l,dmu_dc,snow_index):
    '''Returns the Cl factor which normalises changes in slurry mass fraction
    to changes in temperature.

    Parameters
    ----------
    phi : array
        solid fraction
    L : array
        Latent heat
    T : array
        Temperature
    conc_l : array
        mass fraction of alloying light element
    dmu_dc : arraty
        change in chemical potential with mass fraction
    snow_index : int
        Index from which snow zone starts in arrays

    Returns
    -------
    array
        Cl factor
    '''


    n = snow_index
    Cl = np.zeros(phi.size)

    Cl[n:] = -L[n:]*(1-phi[n:])/(T[n:]*conc_l[n:]*dmu_dc[n:])

    return Cl

def latent_snow(r, rho, L, Cl, conc_l_profile, Ta, snow_index):
    '''
    Returns the normalised rate of latent heat release from changes in mass fraction throughout
    the slurry.

    Parameters
    ----------
    r : array
        Radius
    rho : array
        Density
    L : array
        Latent heat
    Cl : array
        Cl factor relating slurry mass fraction changes to temperature
    conc_l_profile : array
        Light element concentration profile
    Ta : array
        Adiabatic temperature
    snow_index : int
        Index from which snow zone starts in arrays

    Returns
    -------
    (float, float)
        Normalised energy and entropy release
    '''

    n = snow_index

    if np.array(L).size == 1:
        L = np.ones(r.size)*L

    Ql_tilde = 4*np.pi*trapezoid(r[n:], rho[n:]*L[n:]*(Cl[n:]/conc_l_profile[n:])*(Ta[n:]/Ta[0])*r[n:]**2)[-1]

    El_tilde = 4*np.pi*trapezoid(r[n:], rho[n:]*L[n:]*(Cl[n:]/conc_l_profile[n:])*(Ta[n:]/Ta[0])*(1/Ta[-1] - 1/Ta[n:])*r[n:]**2)[-1]

    return Ql_tilde, El_tilde

def Cp_factor(r, rho, Cl, Ta, M_liquid, snow_index):
    '''
    Calculates factor relating changes in mass fraction of the liquid region to those of the slurry

    Parameters
    ----------
    r : array
        Radius
    rho : array
        Density
    Cl : array
        Cl factor relating slurry mass fraction changes to temperature
    Ta : array
        Adiabatic temperature
    M_liquid : float
        Mass of the liquid region
    snow_index : int
        Index from which snow zone starts in arrays

    Returns
    -------
    float
        Cp factor
    '''

    n = snow_index

    Cp = -4*np.pi*trapezoid(r[n:], rho[n:]*Cl[n:]*(Ta[n:]/Ta[0])*r[n:]**2)[-1]/M_liquid

    return Cp

def gravitational_freezing(r, rho, psi, alpha_c, Cl, Ta, snow_index):
    '''
    Calculates the gravitational energy/entropy associated with reducing density of snow zone.

    Parameters
    ----------
    r : array
        Radius
    rho : array
        Density
    psi : array
        Gravitational potential
    alpha_c : float
        Chemical expansivity
    Cl : array
        Cl factor relating slurry mass fraction changes to temperature
    Ta : array
        Adiabatic temperature
    snow_index : int
        Index from which snow zone starts in arrays

    Returns
    -------
    (float, float)
        Normalised energy/entropy
    '''

    n = snow_index

    Qg_tilde = -4*np.pi*trapezoid(r[n:], rho[n:]*psi[n:]*alpha_c*Cl[n:]*(Ta[n:]/Ta[0])*r[n:]**2)[-1]

    Eg_tilde = Qg_tilde/Ta[-1]

    return Qg_tilde, Eg_tilde

def gravitational_melting(r, rho, psi, alpha_c, Cp, Cc, Cr, Tcmb, snow_index):
    '''
    Calculates the gravitational energy/entropy associated with increasing density of liquid region

    Parameters
    ----------
    r : array
        Radius
    rho : array
        Density
    psi : array
        Gravitational potential
    alpha_c : float
        Chemical expansivity
    Cp : float
        Cp factor relating changes in mass fraction of the liquid region to those of the slurry
    Cc : float
        Cc factor relating changes in snow zone radius to changes in mass fraction of the liquid
    Cr : float
        Cr factor relating changes in snow zone radius to changes in temperature
    Ta : array
        Adiabatic temperature
    snow_index : int
        Index from which snow zone starts in arrays

    Returns
    -------
    (float, float)
        Normalised energy/entropy
    '''

    n = snow_index
    if n == 0:
        Qg_tilde, Eg_tilde = 0, 0
    else:
        Qg_tilde = 4*np.pi*trapezoid(r[:n], rho[:n]*psi[:n]*alpha_c*(Cp+(Cc*Cr))*r[:n]**2)[-1]
        Eg_tilde = Qg_tilde/Tcmb

    return Qg_tilde, Eg_tilde