#Core model functions
from thermal_history.utils.optimised_funcs import trapezoid

import numpy as np
from scipy.integrate import trapz,cumtrapz

Na = 6.022140857e23  #Avogadros Constant
###############################################################################
#Latent Heat
###############################################################################
def latent_heat(ri,rho,Ta,Cr,L,idx):
    '''Latent heat of fusion during IC growth

    Parameters
    ----------
    ri : float
        Inner core radius
    rho : array
        Density
    Ta : array
        Adiabatic temperature
    Cr : float
        Cr factor that normalises IC growth rate to cooling rate
    L : array
        Latent heat
    idx : int
        Index for ICB in radial arrays

    Returns
    -------
    (float, float)
        Normlised latent heat and associated entropy.
    '''


    Ql_tilda = 4*np.pi*ri**2*rho[idx]*L[idx]*Cr


    El_tilda = Ql_tilda*(Ta[idx]-Ta[-1])/(Ta[idx]*Ta[-1])


    return Ql_tilda, El_tilda

###############################################################################
#Gravitational Energy
###############################################################################
def gravitational(r, Ta, rho, psi, Cr, Cc, M_oc, alpha_c, idx, rs_idx, Tcmb=False):
    '''Gravitational energy associated with growth of IC

    Parameters
    ----------
    r : array
        radius
    Ta : array
        temperature
    rho : array
        density
    psi : array
        gravitational potential
    Cr : float
        Cr factor that normalises IC growth rate to cooling rate 
    Cc : array
        Cc factor that normalises enrichment of light element in OC to inner core growth
    M_oc : float
        Mass of outer core over which light elements are distributed
    alpha_c : array
        chemical expansivities for light elements
    idx : int
        Index for ICB in radial arrays
    rs_idx : int
        Index for rs (base of stable layer) in radial arrays
    Tcmb : bool, optional
        temperature of the CMB. If False, then the end of the Ta array is used, by default False

    Returns
    -------
    (float, float)
        Normlised gravitational energy and associated entropy
    '''

    I = 4*np.pi*trapezoid(r[idx:rs_idx+1],rho[idx:rs_idx+1]*psi[idx:rs_idx+1]*r[idx:rs_idx+1]**2)[-1] - M_oc*psi[idx]

    Qg_tilda = I*np.sum(alpha_c*Cc)*Cr

    if Tcmb == False:
        Eg_tilda = Qg_tilda/Ta[-1]
    else:
        Eg_tilda = Qg_tilda/Tcmb

    return Qg_tilda, Eg_tilda

###############################################################################
#Gravitational Energy for Magnesium precipitation from Du et al. 2017
###############################################################################
def gravitational_precip(r, Ta, rho, psi, idx, Cm, alpha_c):
    '''Gravitational energy associated with precipitation

    Parameters
    ----------
    r : array
        radius
    Ta : array
        temperature
    rho : array
        density
    psi : array
        gravitational potential
    Cm : float
        Factor that normlises precipitation to cooling rate
    alpha_c : float
        chemical expansivity of precipitating element
    idx : int
        Index for ICB in radial arrays

    Returns
    -------
    (float, float)
        Normlised gravitational energy and associated entropy
    '''

    I = 4*np.pi*trapezoid(r[idx:],rho[idx:]*psi[idx:]*r[idx:]**2)[-1]# - M_oc*psi[idx] assuming no uneven partitioning of MgO

    Qg_tilda = I*alpha_c*Cm

    Eg_tilda = Qg_tilda/Ta[-1]

    return Qg_tilda, Eg_tilda

###############################################################################
#Secular cooling
###############################################################################
def secular_cool(r,rho,Ta,cp,idx, Tcmb=False):
    '''Energy from change in sensible energy from secular cooling

    Parameters
    ----------
    r : array
        radius
    rho : array
        density
    Ta : array
        temperature
    cp : array
        specific heat capacity
    idx : int
        Index for rs (base of stable layer) in radial arrays
    Tcmb : bool, optional
        temperature of the CMB. If False, then the end of the Ta array is used, by default False

    Returns
    -------
    (float, float)
        Secular cooling energy and associated entropy
    '''

    #Values are normalised to the cooling rate.
    Tcen = Ta[0]

    if Tcmb == False:
        Tcmb = Ta[-1]

    I1 = 4*np.pi*trapezoid(r[:idx+1],rho[:idx+1]*cp[:idx+1]*Ta[:idx+1]*r[:idx+1]**2)[-1]

    Qs = -I1/Tcen

    I2 = 4*np.pi*trapezoid(r[:idx+1], (Ta[:idx+1]/Tcmb - 1) * rho[:idx+1]*cp[:idx+1]*r[:idx+1]**2)[-1]

    Es = -I2/Tcen

    return Qs, Es
###############################################################################
#Radiogenic heating
###############################################################################
def radiogenic_heating(time,r,rho,Ta,Mass,h0,decay_time):
    '''Heat from radiogenic heating (NOT TESTED)

    Parameters
    ----------
    time : float
        time (relative to present day)
    r : array
        radius
    rho : array
        density
    Ta : array
        temperature
    Mass : float
        mass of core
    h0 : float
        Present day heating rate per unit mass
    decay_time : float
        Half life of radiogenic element

    Returns
    -------
    (float, float)
        Radiogenic heating and associated entropy
    '''
    
    h = h0 * 2**(-time/decay_time)
    Qr = Mass*h

    Er = 4*np.pi*trapezoid(r,h*rho*((1/Ta[-1])-(1/Ta))*r**2)[-1]

    return Qr, Er
###############################################################################
#Entropy of conduction
###############################################################################
def cond_entropy(Ta,Ta_grad,k,r,idx):
    '''Entropy of thermal conduction in adiabatic region

    Parameters
    ----------
    Ta : array
        Adiabatic temperature
    Ta_grad : array
        Adiabatic temperature gradient
    k : array
        Thermal conductivity
    r : array
        radius
    idx : int
        Index for rs (base of stable layer) in radial arrays

    Returns
    -------
    float
        Entropy of thermal conduction
    '''

    Ek = 4*np.pi*trapezoid(r[:idx+1],k[:idx+1]*(Ta_grad[:idx+1]/Ta[:idx+1])**2*r[:idx+1]**2)[-1]
    #pdb.set_trace()
    return Ek
###############################################################################
#Entropy of mixing
###############################################################################
def heat_of_reaction(rho_oc,dmu_dt_poly,Cc,Cr,r_oc,mm):
    '''Entropy of heat of reaction at ICB (NOT TESTED IN A WHILE)

    Parameters
    ----------
    rho_oc : array
        Outer core density
    dmu_dt_poly : array
        radial polynomials for change in chemical potential with temperature.
    Cc : array
        Cc factor that normalises enrichment of light element in OC to inner core growth
    Cr : float
        Cr factor that normalises IC growth rate to cooling rate 
    r_oc : array
        outer core radial grid points
    mm : array
        molar masses of iron and light elements

    Returns
    -------
    float
        Entropy from heat of reaction
    '''

    Eh_tilda = 0
    for i in range(Cc.size):

        dmu_dt = np.polyval(dmu_dt_poly[i][::-1],r_oc)*Na*1000/mm[i+1]
        Eh_tilda += -4*np.pi*trapezoid(r_oc,dmu_dt*rho_oc*r_oc**2)[-1]*Cc[i]*Cr

    return Eh_tilda

###############################################################################
#Entropy of mass diffusion
###############################################################################
def mass_diffusion(r, i, alpha_D, T):
    '''Entropy from barodiffusion of light elements

    Parameters
    ----------
    r : array
        radius
    i : float or array
        mass flux
    alpha_D : array
        barodiffusion coefficient
    T : array
        temperature

    Returns
    -------
    float
        Total entropy of barodiffusion of light elements.
    '''

    if np.max(alpha_D) == 0:
        E_alpha = 0  #alpha_D = 0 if no light element is present
    else:
        E_alpha = 4*np.pi*trapezoid(r,r**2*i**2/(alpha_D*T))[-1]



    return E_alpha

###############################################################################
#Mass flux
###############################################################################
def mass_flux(rho, alpha_c, alpha_D, D, c_grad, g):
    '''Mass flux from diffusion down chemical and pressure gradients

    Parameters
    ----------
    rho : float or array
        density
    alpha_c : float
        chemical expansivities for light element
    alpha_D : float
        barodiffusion coefficients for light element
    D : float
        Self diffusion coefficients for light element
    c_grad : float or array
        chemical gradient
    g : float or array
        gravity

    Returns
    -------
    float or array
        mass flux
    '''
    i = -rho*D*c_grad + alpha_c*alpha_D*g
    return i

###############################################################################
#Numerical integration
###############################################################################
def integrate(x,y,cumulative=0):
    '''
    Depreciated in favour of thermal_history.utils.optimised_funcs.trapezoid
    '''
    if cumulative == 0:
        A = trapz(y,x)
    else:
        A = cumtrapz(y,x,initial=0)
    return A
