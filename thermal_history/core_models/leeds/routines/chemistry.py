#Core chemistry functions
import pdb
from numba import njit

from thermal_history.utils.optimised_funcs import polyval

from . import profiles as prof

import numpy as np
from scipy.optimize import bisect

kb = 1.3806485e-23   #Boltzmanns constant
ev = 1.602e-19        #Electron volt
Na = 6.022140857e23   #Avogadros Constant



def melting_curve(model):
    """Calculates the melting curve and sets the mass fraction of light element in the solid iron core (conc_s) according to
    the relevant partitioning model.
    Different methods for the paramterisation of the melting curve can be set by the core_melting_params parameter.
    In this parameter (a List), the first element is a string, followed by the constants (floats) required for the
    parameterisation of conc_s and the melting curve. The accepted strings are currently:

    'AL': This method is also the default if no string is given. Uses the method of Alfe et al. (2002) to determine conc_s
    based on equal chemical potentials at the ICB and the depression of the iron melting point for the alloy. The iron melting temperature is calculated using
    polynomials in pressure, with the coefficients given after 'AL' in core_melting_params.
    E.g. = core_melting_params = ['AL', Tl0, Tl1, Tl2.....]

    'SG': Uses the Simon-Glatzel parameterisation of the melting curve with fractionation of elements given by fixed
    partitioning behaviour given by the partition_coeff parameter. After 'SG', the parameters required for the melting
    curve are the constants Tl0, beta1 and beta2 (see Knibbe and Van Westrenen (2018) for details).
    E.g. core_melting_params = ['SG', Tl0, beta1, beta2]

    'WN': Uses the Williams and Nimmo (2004) parameterisation (also used in Davies and Pommier (2018)) for the melting curve.
    Like 'SG', the partitioning is also described by fixed parition coefficients (=conc_s/conc_l) in the partition_coeff parameter.
    After 'WN', the constants Tl0 and the polynomial coefficients in pressure should be given. E.g. in Davies and Pommier (2018)
    this would be represented as core_melting_params = ['WN', Tl0, 1, Tl1, Tl2]

    Parameters
    ----------
    model : ThermalModel
        ThermalModel class

    Returns
    -------
    (numpy.ndarray, numpy.ndarray, numpy.ndarray)
        Tuple of 3 arrays for the melting point of iron, melting point of the alloy, and derivative of the alloy melting point with pressure.
        Values are calculated on the same radial grid as all core profiles.

    Raises
    ------
    ValueError
        If an unrecognised string in core_melting_params is given.
    """

    prm = model.parameters
    core = model.core
    P = core.profiles['P']
    melting_params = prm.core_melting_params

    #Change in entropy on melting iron
    core.profiles['dS'] = prof.entropy_melting(core.profiles['P'], prm.entropy_melting_params)

    #Decide on which melting curve parameterisation to use based on 2 letter string at
    #beginning of melting_params. If no string is given, assume values are for Alfe parameterisation.
    #Can add more conditional statements to consider more parameterisations in the furture.

    #Alfe 2002 method
    if melting_params[0] == 'AL' or type(melting_params[0].item()) in (float, int):
        if melting_params[0]=='AL':
            params = melting_params[1:].astype('float64') #Don't pass through first string
        else:
            params = melting_params

        Tm_fe = iron_melting(P, params)

        #Calculate fractionation of LE
        if prm.use_partition_coeff:
            #Calc using fixed partition coefficients
            core.conc_s = core.conc_l * prm.partition_coeff
        else:
            #Calc conc_s based on equal chemical potential
            core.conc_s = LE_frac_dep(P, Tm_fe, core._ri_idx, core.conc_l, prm.mm,
                                       prm.dmu, prm.lambda_liq, prm.lambda_sol, prm.entropy_melting_params)

        core.mf_s = mass_conc2mole_frac(core.conc_s, prm.mm)



        #Melting point depression. Assumes constant dTm across the whole core since we use _ri_idx.
        dTm = melt_pt_dep(core.mf_l, core.mf_s, Tm_fe[core._ri_idx], core.profiles['dS'][core._ri_idx])

        dTm_dP = iron_melting_gradient(P, params)

    #Simon Glatzel
    elif melting_params[0] == 'SG':

        #Simple partitioning
        core.conc_s = core.conc_l * prm.partition_coeff

        params = melting_params[1:].astype('float64') #Don't pass through first string
        Tm_fe = simon_glatzel(P, params)
        dTm   = 0
        dTm_dP = simon_glatzel_gradient(P, params)

    #Williams and Nimmo for iron snow. Needs conc_l in core.profiles.
    elif melting_params[0] == 'WN':

        #Simple partitioning
        core.conc_s = core.conc_l * prm.partition_coeff

        params = melting_params[1:].astype('float64') #Don't use first string
        Tl0 = params[0]

        #Set up conc_l profile on first iteration
        if not 'conc_l' in core.profiles.keys():
            core.profiles['conc_l'] = np.full(prm.n_profiles, core.initial_conc_l, dtype=float)

        if len(params)>1:
            Tm_fe = Tl0*(1+core.initial_conc_l)*polyval(params[1:][::-1], P)
        else: #If Tm has no dependence on P
            Tm_fe = np.full(P.size, Tl0*(1+core.initial_conc_l))

        dTm   = -Tm_fe*core.profiles['conc_l']/(1+core.initial_conc_l)

        dTm_dP = Tl0*(1+core.initial_conc_l-core.profiles['conc_l'])*iron_melting_gradient(P, params[1:])

    else:

        raise ValueError(f'Incorrect string denoting which melting parameterisation to use. See this function (core_models.leeds.routines.chemisty.melting_curve) for valid options. Supplied core_melting_params={melting_params}')

    Tm = Tm_fe + dTm

    return Tm_fe, Tm, dTm_dP

def iron_melting(P, melting_params):
    """Simple polynomial representation of iron melting temperature with pressure

    Parameters
    ----------
    P : int, float, numpy.ndarray
        Pressure values to evaluate at.
    melting_params : numpy.ndarray
        polynomial coefficients in increasing order.

    Returns
    -------
    float, np.ndarray (same shape as P)
        Melting temperature of iron
    """
    return polyval(melting_params[::-1], P)

def iron_melting_gradient(P, melting_params):
    '''Melting temperature gradient

    Parameters
    ----------
    P : float or numpy array
        Pressure value to evalulate gradient at
    melting_params : numpy array
        array containing coefficients of pressure polynomials for melting temperature

    Returns
    -------
    float or numpy array (same type as P)
        Melting temperature gradient evaluated at P.
    '''

    if melting_params.size > 1:
        poly = np.zeros(melting_params.size-1)
        for i in range(poly.size):
            poly[i] = melting_params[i+1] * (i+1)
        return polyval(poly[::-1], P)
    else:
        return P*0

def simon_glatzel(P, melting_params):
    '''Melting temperature of iron alloy parameterised with Simon-Glatzel's equation.

    Parameters
    ----------
    P : float or numpy array
        Pressure value to evalulate gradient at
    melting_params : numpy array
        array containing 3 values: the reference temperature, beta1 and beta2.

    Returns
    -------
    float or numpy array
        The melting temperature at P
    '''

    T0, beta1, beta2 = melting_params
    return T0 * (P/(1e9*beta1) + 1)**(1/beta2)

def simon_glatzel_gradient(P, melting_params):
    '''Gradient of the melting temperature of iron alloy parameterised with Simon-Glatzel's equation.

    Parameters
    ----------
    P : float or numpy array
        Pressure value to evalulate gradient at
    melting_params : numpy array
        array containing 3 values: the reference temperature, beta1 and beta2.

    Returns
    -------
    float or numpy array
        The melting temperature gradient at P
    '''
    T0, beta1, beta2 = melting_params
    return (T0/(beta1*beta2)) * (P/(1e9*beta1) + 1)**(1/beta2 - 1)

###############################################################################
#Convert mole fraction to mass concentration by equations in Labrosse (2014)
###############################################################################
def mole_frac2mass_conc(mf, mm):
    '''Converts from mole fraction to mass fraction for core alloying light elements

    Parameters
    ----------
    mf : array
        mole fractions of light elements
    mm : array 
        molar masses of iron followed by light elements

    Returns
    -------
    array
        mass fractions of light elements (same size as mf)
    '''

    mm = np.array(mm)

    denom = np.dot(mf,mm[1:]) + (1 - np.sum(mf))*mm[0]

    conc = np.array(mf)*np.array(mm)[1:]/denom

    return conc

###############################################################################
#Convert mass concentration to mole fraction by equations in Labrosse (2014)
###############################################################################
def mass_conc2mole_frac(conc, mm):
    '''Converts from mass fraction to mole fraction for core alloying light elements

    Parameters
    ----------
    conc : array
        mass fractions of light elements
    mm : array 
        molar masses of iron followed by light elements

    Returns
    -------
    array
        mole fractions of light elements (same size as conc)
    '''

    mm = np.array(mm)

    denom = np.sum(conc/mm[1:]) + (1-np.sum(conc))/mm[0]

    mf = np.array(conc)/(np.array(mm)[1:]*denom)

    return mf

###############################################################################
#Calculate the mole fraction of light element species in the solid
###############################################################################
def solid_conc(mf_liq, T_m, ds_fe, dmu, lambda_liq, lambda_sol):
    '''Calculates the mole fraction of light element in the inner core from chemical equilibria

    Parameters
    ----------
    mf_liq : array
        mole fraction of light elements in liquid core
    T_m : float
        melting temperature of iron at ICB
    ds_fe : float
         entropy of melting of iron at ICB
    dmu : array
        change in chemical potentials between solid and liquid
    lambda_liq : array
        corrections to chemical potentials in liquid
    lambda_sol : array
        correctinos to chemical potentials in solid

    Returns
    -------
    array
        mole fractions in the inner core
    '''

    def f(guess,dmu_x,lambda_liq_x,lambda_sol_x,mf_liq_x,T_m,ds_fe):
        return dmu_x + mf_liq_x*lambda_liq_x - guess*lambda_sol_x - kb*T_m*np.log(guess/mf_liq_x)*(1+(guess-mf_liq_x)/(ds_fe/kb))


    #Bisection method to find chemcical equilibrium
    mf_sol = np.zeros(mf_liq.size)
    for i in range(mf_liq.size):
        if mf_liq[i] == 0:
            mf_sol[i] = 0
        else:
            lower = 1e-10
            upper = 1

            try:
                mf_sol[i] = bisect(f,lower,upper,args=(dmu[i],lambda_liq[i],lambda_sol[i],mf_liq[i],T_m,ds_fe),maxiter=200)
            except:
                breakpoint()
                raise ValueError('Can\'t solve for mole fraction in the solid. Check radial polynomials for melting temperature and entropy of freezing')


    return mf_sol


###############################################################################
#Calculate the melting temperature including the depression due to light elements
#Equation 12 from Alfe et al.(2002)
###############################################################################
def melt_pt_dep(mf_liq, mf_sol, Tm, ds_fe):
    '''Calculates the melting point depression by Alfe et al. (2002)

    Parameters
    ----------
    mf_liq : array
        mole fraction of light elements in the liquid core
    mf_sol : array
        mole fraction of light elements in the solid core
    Tm : float or array
        Melting temperature
    ds_fe : float or array
        entropy change of melting at Tm

    Returns
    -------
    float or array
        Deflection of melting curve due to presence of light elements
    '''


    dTm = kb*(Tm/ds_fe)*np.sum(mf_sol-mf_liq)

    return dTm

###############################################################################
#Calculate the melting temperature depression and fractionation of light elements.
###############################################################################
def LE_frac_dep(P, Tm_fe, ri_idx, conc_l, mm, dmu, lambda_liq, lambda_sol, ent_mel_poly):
    '''Calculates the concentration of light elements in the solid inner core

    Parameters
    ----------
    P : float
        Pressure at the inner core boundary
    Tm_fe : float
        Melting temperture of pure iron at the ICB
    ri_idx : int
        index for ICB in radial arrays
    conc_l : array
        mass fraction of light element in the liquid
    mm : array
        molar masses of iron and the light elements
    dmu : array
        change in chemical potential between solid/liquid
    lambda_liq : array
        correction to chemical potentials from an ideal solution for liquid
    lambda_sol : array
        correction to chemical potentials from an ideal solution for solid
    ent_mel_poly : array
        radial polynomials for entropy of melting

    Returns
    -------
    array
        mass fraction of light element in the solid
    '''


    #fractionation of light elements
    ds_fe = prof.entropy_melting(P, ent_mel_poly)
    mf_l = mass_conc2mole_frac(conc_l, mm)
    mf_s = solid_conc(mf_l, Tm_fe[ri_idx], ds_fe[ri_idx], dmu, lambda_liq, lambda_sol)

    conc_s = mole_frac2mass_conc(mf_s, mm)

    return conc_s


###############################################################################
#Calculate barodiffusion coefficient as defined in Gubbins et al. (2004)
###############################################################################
def baro_coeff(mf_l, mm, T_av, rho_av, lambda_liq, diffusivity_c):
    '''Calculates the barodiffusion coefficient

    Parameters
    ----------
    mf_liq : array
        mole fraction of light elements in the liquid core
    mm : array
        molar masses of iron and the light elements
    T_av : float
        Average temperature to evaluate at
    rho_av : float
        Average density to evaluate at
    lambda_liq : array
        correction to chemical potentials from an ideal solution for liquid for light elements
    diffusivity_c : array
        Self diffusion coefficients for light elements

    Returns
    -------
    array
        barodiffusion coefficients for light elements
    '''

    if not type(mf_l) == np.ndarray:
        mf_l = np.array([mf_l])

    return baro_coeff_fast(np.array(mf_l, dtype='f8'), np.array(mm, dtype='f8'), T_av, rho_av, lambda_liq, diffusivity_c)

@njit
def baro_coeff_fast(mf_l, mm, T_av, rho_av, lambda_liq, diffusivity_c):
    '''Faster jit compiled version of baro_coeff

    Parameters
    ----------
    mf_liq : array
        mole fraction of light elements in the liquid core
    mm : array
        molar masses of iron and the light elements
    T_av : float
        Average temperature to evaluate at
    rho_av : float
        Average density to evaluate at
    lambda_liq : array
        correction to chemical potentials from an ideal solution for liquid for light elements
    diffusivity_c : array
        Self diffusion coefficients for light elements

    Returns
    -------
    array
        barodiffusion coefficients for light elements
    '''

    # A_av = 0
    # for i in range(mm.size-1):
    #     A_av += np.dot(mf_l, mm[1:]) + (1 - np.sum(mf_l))*mm[0]


    A_av = np.dot(mf_l,mm[1:]) + (1 - np.sum(mf_l))*mm[0]

    dmu_dmf = kb*T_av/mf_l + lambda_liq

    dmu_dc  = dmu_dmf * (A_av/mm[1:]) * Na * (1000/mm[1:])

    alpha_D = rho_av*diffusivity_c/dmu_dc

    return alpha_D


###############################################################################
#Set Ta(r=0) for the specified inner core size from the model_class.state()
###############################################################################
def calibrate_temperature(model):
    '''Calculates the correct value for Tcen for a given inner core radius and melting curve

    Sets model.core.T_cmb, as well as the adiabat/adiabatic gradient profiles with correct temperature
    for given inner core radius

    Parameters
    ----------
    model : ThermalModel class
        Main ThermalModel class being used for the calculation.
    '''


    prm = model.parameters
    core = model.core

    #Get melting curve
    Tm = core.profiles['Tm']
    poly =  prm.core_adiabat_params

    Tm_ri = Tm[core._ri_idx]

    #Calculate necessary Tcen for given ri and Tm
    core.Tcen = Tm_ri/prof.adiabat(core.ri, 1, poly)

    #Reset adiabat profiles
    core.profiles['Ta']     = prof.adiabat(     core.profiles['r'], core.Tcen, poly)
    core.profiles['dTa_dr'] = prof.adiabat_grad(core.profiles['r'], core.Tcen, poly)

    #Reset CMB temperature
    core.T_cmb = core.profiles['Ta'][-1]
