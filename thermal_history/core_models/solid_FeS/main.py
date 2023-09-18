#Description of this model3
Description = 'Based on the leeds model with adaptions for Solid FeS layer as in Rückriemen et al. (2018).'

import numpy as np
import copy

from ..leeds.routines import chemistry as chem    #Import this model's specific version of chemistry
from .routines.chemistry import melting_curve
chem.melting_curve = melting_curve

from thermal_history.utils.optimised_funcs import linspace, polyval, trapezoid
from ..leeds.routines import profiles as prof
from ..leeds.routines import energy as en
from ..leeds.routines import snow as snow
from ..leeds.routines import rivoldini_eos as eos

import logging
logger = logging.getLogger(__name__)


#Required functions: setup(), evolution(), update()

#Use functions from leeds core model. Just the snow_evolution, update and set_Q_rs functions are modified in comparison
from thermal_history.core_models.leeds.main import setup, evolve, progress, required_params, optional_params
import thermal_history.core_models.leeds.main as leeds


#Redefine snow evolution for solid FeS then overwrite leeds.snow_evolution.
def snow_evolution(model):
    '''
    Sub-evolution fuction specifically for a freezing solid FeS layer as in Rückriemen et al. (2018). Called from evolution() if required (prm.iron_snow=True).
    Modified version of the iron snow routine. Differences have been highlighted in the comments with a \*.

    Parameters
    ----------
    model : ThermalModel Class
        Main model

    Returns
    -------
    Dict
        Dictionary containing relevant energies/entropies/other variables to be passed back to the main evolution function.
    '''

    prm = model.parameters
    core = model.core

    rho_l = prm.core_liquid_density_params
    rho_s = prm.core_solid_density_params

    # #Initialise snow parameters Cl, concl_l and phi_snow profiles
    Cl = np.zeros(prm.n_profiles)
    phi_snow = np.zeros(prm.n_profiles)
 
    T, Tm, Tm_fe   = core.profiles['T'], core.profiles['Tm'], core.profiles['Tm_fe']
    r, rho, P, psi, g = core.profiles['r'], core.profiles['rho'], core.profiles['P'], core.profiles['psi'], core.profiles['g']

    #Initialise conc_l profile if not already present
    if not 'conc_l' in core.profiles.keys():
        conc_l_profile = np.full(r.size, core.conc_l[0])
    else:
        conc_l_profile = core.profiles['conc_l']

    #Indicies for key radii
    ri_idx   = core._ri_idx
    snow_idx = core._snow_idx
    rs_idx   = core._rs_idx

    if core.r_snow < prm.r_cmb and not core.r_snow==0:
        # assert hasattr(core, 'dT_dt'), 'No cooling rate has been calculated yet on iteration {}'.format(model.it)

        # *** NOT ASSUMED ***
        # #Concentration of slurry to depress Tm to adiabatic temperature
        # conc_l_profile[snow_idx:] = snow.snow_composition(T, Tm_fe, core.initial_conc_l, snow_idx, core.conc_l[0], P, prm.core_melting_params)
        
        if prm.core_melting_params[0] == 'XX': 
            conc_l_profile[snow_idx:] = float(chem.mole_frac2mass_conc([1.0], prm.mm)) #Mass fraction for pure C
        else:
            conc_l_profile[snow_idx:] = float(chem.mole_frac2mass_conc([0.5], prm.mm)) #Mass fraction for pure FeS
        
        core.profiles.update({'conc_l': conc_l_profile})

        # #Move melting temperature accordingly and update latent heat
        Tm_fe, Tm, dTm_dP = chem.melting_curve(model) #More general function

        # *** NOT ASSUMED ***
        # #Melting temperature constrained to temperature profile with iron snow
        # Tm[core._snow_idx:] = core.profiles['T'][core._snow_idx:]

        #Latent heat in J/kg
        if prm.core_latent_heat is None:
            L = core.profiles['dS'] * prm.Na *1000/prm.mm[0] * Tm   #Latent heat based on dS for iron
        else:
            L = np.ones(Tm.size)*prm.core_latent_heat               #Fixed value for latent heat

        M_liquid = prof.mass(rho_s, rho_l, 0, core.r_snow, core.r_snow)[-1] #Mass of liquid region

        #Cc factor, analagous to inner core, relates changing mass fraction in liquid region to movement of the snow zone radius.
        Cc_snow = 4*np.pi*r[snow_idx]**2*rho[snow_idx]*(conc_l_profile[snow_idx]-core.conc_l[0])/M_liquid


        # *** NOT NEEDED ***
        # Cl factor, normalising changes in slurry mass fraction to changes in temperature (DP eq. 22)
        # although it could be calculated by derivative of WN melting temp (since it equals adiabat in snow zone)
        # Cl[snow_idx:] = -1/Tm_fe[snow_idx:]

        # dmu_dc = (prm.kb*T/conc_l_profile)*prm.Na*(1000/prm.mm[1])
        # Cl = snow.Cl_calc(phi_snow, L, T, conc_l_profile, dmu_dc, snow_idx)
        Cl = np.zeros(r.size)

        #Solid fraction in slurry formed in one timestep. Needed to calculate Ql_snow, total latent heat of moving boundary.
        # phi_snow = (Cl/conc_l_profile) * (T/core.Tcen) * core.dT_dt*model.dt
        phi_snow = np.ones(r.size)

        #Normalised snow zone radius velocity to cooling rate, same as done for inner core
        dTm_dr = -rho[snow_idx]*g[snow_idx]*dTm_dP[snow_idx]
        # dTm_dr = (Tm[snow_idx-2]-Tm[snow_idx-1])/(r[snow_idx-2]-r[snow_idx-1])
        dT_dr  = (T[snow_idx-2]-T[snow_idx-1])/  (r[snow_idx-2]-r[snow_idx-1])

        # dTm_dc = -Tm_fe[snow_idx]/(1+core.initial_conc_l) #Assumes WN melting curve
        if prm.core_melting_params[0] == 'RI': #Use Rivoldini data
            dTm_dc = chem.riv_dTm_dc(P[snow_idx], core.conc_l[0], prm.core_melting_params[1:].astype('float64'))
            x_eu   = eos.xeFeS(1e-9*P[snow_idx])

        # *** INCLUDED RUCKRIEMEN PARAMETERISATION ***
        elif prm.core_melting_params[0] == 'RU': #Ruckriemen parameterisation

            params = prm.core_melting_params[1:].astype('float64') #Don't use first string

            T_eutectic = eos.TeFeS(P[snow_idx]/1e9)
            Tm_fes = params[3] + params[4]*P[snow_idx] - params[5]*P[snow_idx]**2
            x_eu   = eos.xeFeS(1e-9*P[snow_idx])
            
            if not all(np.greater_equal(core.profiles["conc_l"],eos.xeFeS(1e-9*P))):        
                model.critical_failure = True
                model.critical_failure_reason = 'Reached Fe-S eutectic!'
                logger.critical(f'it: {model.it}. Reached Fe-S eutectic!')
                       
            dTm_dc = (Tm_fes-T_eutectic)/(0.3647-x_eu)

        elif prm.core_melting_params[0] == 'XX':
            FeC=prm.core_melting_params[1]
            dTm_dc=FeC.dTmdx(core.conc_l[0],P[snow_idx]*1e-9)
            # check if C in convecting core is in the graphite stability field
            if T[0] <= FeC.Txg(P[snow_idx]*1e-9) or core.conc_l[0]<=FeC.xg(P[0]*1e-9) :
                model.critical_failure =True 
                model.critical_failure_reason ='Outside graphite stability zone in r=0!'
                logger.critical(f'it:{model.it}. Outside graphite stability zone in r=0!')
            if not all(np.greater_equal(core.profiles["conc_l"],eos.pdFeC().xg(1e-9*P))):
                model.critical_failure = True
                model.critical_failure_reason = 'Composition outside stability field!'
                logger.critical(f'it: {model.it}. Composition outside stability field!')
                
        else:
            dTm_dc = 0

 
        if prm.use_new_Cr:
            Cr_snow = (1/(dTm_dr - dT_dr + dTm_dc*Cc_snow))*(T[snow_idx]/core.Tcen)
        else:
            Cr_snow = (1/(dTm_dr - dT_dr))*(T[snow_idx]/core.Tcen)


        # *** NO COMPLEX LATENT HEAT ***
        #Normalised latent heat from moving snow zone radius (DP18 Qlb eq 14). Analagous to inner core growth.
        # Ql_snow_tilde, El_snow_tilde = en.latent_heat(core.r_snow, rho*phi_snow[snow_idx], T, Cr_snow, L, snow_idx) #Multiply density by phi for just mass of freezing solid.

        #Latent heat relase and absorbtion upon freezing and melting. Ql_melting is given negative L to indicate heat absorption.
        #(DP18 Qls and Qll)
        # Ql_freezing_tilde, El_freezing_tilde = snow.latent_snow(r, rho, L,            Cl, conc_l_profile, T, snow_idx)
        # Ql_melting_tilde,  El_melting_tilde  = snow.latent_snow(r, rho, -L[snow_idx], Cl, conc_l_profile, T, snow_idx)

        # *** Simple latent heat at advancing boundary ***
        Ql_snow_tilde, El_snow_tilde = en.latent_heat(core.r_snow, rho, T, -Cr_snow, L, snow_idx) #Negative Cr because these fncts are for inner core

        Ql_freezing_tilde, El_freezing_tilde = 0, 0
        Ql_melting_tilde,  El_melting_tilde = 0, 0
        
        # *** NOT NEEDED ***
        #Cp factor, normalises changing mass fraction in liquid region to changes in slurry mass fraction (which in turn uses Cl to relate to cooling)
        #Uses DP18 eqn 20
        # Cp_snow = snow.Cp_factor(r, rho, Cl, T, M_liquid, snow_idx)
        Cp_snow = 0


        # *** Just gravitational energy in the bulk of the core, handled by CC_snow modifying dc_dt in bulk.
        #Gravitational energy release from removal of solid by melting, then increasing density of the liquid region (DP18 eqn 16/17)
        # Qg_freezing_tilde, Eg_freezing_tilde = snow.gravitational_freezing(r, rho, psi, prm.alpha_c[0], Cl, T, snow_idx)
        # Qg_melting_tilde,  Eg_melting_tilde  = snow.gravitational_melting(r, rho,  psi, prm.alpha_c[0], Cp_snow, Cc_snow, Cr_snow, T[-1], snow_idx)

        Qg_freezing_tilde, Eg_freezing_tilde = 0, 0
        Qg_melting_tilde,  Eg_melting_tilde = en.gravitational(r, T, rho, psi, -Cr_snow, Cc_snow, M_liquid, prm.alpha_c, 0, snow_idx, Tcmb=core.T_cmb) #Negative Cr because these fncts are for inner core
        Qg_snow_tilde, Eg_snow_tilde = 0, 0

        # I showed in the notes on the gravitational energy than this term should include the jump in conc_l
        # not just conc_l as in DP18.
        # Qg_snow_tilde = -4*np.pi*r[snow_idx]**2*rho[snow_idx]*psi[snow_idx]*(conc_l_profile[snow_idx]-conc_l_profile[snow_idx-1])*Cr_snow
        # Eg_snow_tilde = Qg_snow_tilde/T[-1]

        #Set for the correct size (total size of number of LE) to avoid inconsistencies down the line.
        Cc_snow = np.ones(core.conc_l.size)*Cc_snow
        Cc_snow[1:] = 0 #Assumes snow_evolution is only considering 1st light element.

        # Ql_snow_tilde, Ql_freezing_tilde, Ql_melting_tilde, Qg_freezing_tilde, Qg_melting_tilde, Qg_snow_tilde = 0,0,0,0,0,0

        # breakpoint()

    else:
        L = prof.entropy_melting(P, prm.entropy_melting_params) * prm.Na *1000/prm.mm[0] * Tm
        Cc_snow, Cr_snow, Cp_snow = np.zeros(core.conc_l.size), 0, 0
        Ql_snow_tilde, Ql_freezing_tilde, Ql_melting_tilde, Qg_freezing_tilde, Qg_melting_tilde, Qg_snow_tilde = 0,0,0,0,0,0
        El_snow_tilde, El_freezing_tilde, El_melting_tilde, Eg_freezing_tilde, Eg_melting_tilde, Eg_snow_tilde = 0,0,0,0,0,0

        if core.r_snow == 0:
            model.critical_failure = True   #Set flag that critical failure has occured.
            model.critical_failure_reason = 'Snow zone covers whole core'
            logger.critical(f'it: {model.it}. Snow zone has encompassed entire core, not defined how to continue!')


    #Normalised energies and entropies associated with the snow zone. Values to be returned to the main evolution function
    snow_dict = dict(zip(['Ql_b', 'Ql_s', 'Ql_l', 'Qg_b', 'Qg_s', 'Qg_l'], [Ql_snow_tilde, Ql_freezing_tilde, Ql_melting_tilde, Qg_snow_tilde, Qg_freezing_tilde, Qg_melting_tilde]))
    snow_dict.update(dict(zip(['El_b', 'El_s', 'El_l', 'Eg_b', 'Eg_s', 'Eg_l'], [El_snow_tilde, El_freezing_tilde, El_melting_tilde, Eg_snow_tilde, Eg_freezing_tilde, Eg_melting_tilde])))

    snow_dict['Cr'] = Cr_snow
    snow_dict['Cc'] = Cc_snow
    snow_dict['Cp'] = Cp_snow

    #Update profiles
    core.profiles.update({'L': L, 'Cl': Cl, 'phi_snow': phi_snow, 'Tm': Tm, 'conc_l': conc_l_profile, 'T': T})

    return snow_dict

leeds.snow_evolution = snow_evolution

def update(model):
    '''
    Updates parameters based on rates of change calculatied in main fuction.

    Parameters
    ----------
    model : ThermalModel Class
        Main model
    '''

    prm = model.parameters
    core= model.core
    profiles = core.profiles

    #Update with RoC
    core.Tcen   += model.dt*core.dT_dt
    #core.Tcmb    =profiles['T'][-1] # not updated in final result, why?
    core.conc_l += model.dt*core.dc_dt
    core.mf_l    = chem.mass_conc2mole_frac(core.conc_l, prm.mm)
    core.mf_s    = chem.mass_conc2mole_frac(core.conc_s, prm.mm)

    #Update conc_l profile
    if 'conc_l' in profiles.keys():
        profiles['conc_l'][:core._snow_idx] = core.conc_l[0]

    #Update temperature/melting curve to find new snow/ICB radius
    profiles['T'] = prof.adiabat(profiles['r'], core.Tcen, prm.core_adiabat_params)
    if prm.stable_layer:
        sl_prof =  model.stable_layer._next_profiles
        profiles['T'][core._rs_idx:] = np.interp(profiles['r'][core._rs_idx:], sl_prof['r'], sl_prof['T'])

    profiles['Tm_fe'], profiles['Tm'], profiles['dTm_dP'] = chem.melting_curve(model)

        
    #New snow/ICB radius
    if prm.iron_snow:
        idx = core._snow_idx
        core.r_snow = snow.snow_radius(profiles['r'][:idx+1], profiles['T'][:idx+1], profiles['Tm'][:idx+1])
        # if core.r_snow < prm.r_cmb:
        #     breakpoint()

    else:
        core.ri = prof.ic_radius(profiles['r'],
                                prof.adiabat(profiles['r'], core.Tcen, prm.core_adiabat_params),
                                profiles['Tm'])

        if core.ri == prm.r_cmb:
            model.critical_failure = True
            logger.critical(f'it: {model.it}. Inner core has covered entire core!! Not defined how to procede')

    #Update profiles on new grid
    prof.basic_profiles(model)
    prof.temp_dependent_profiles(model)
    profiles['Tm_fe'], profiles['Tm'], profiles['dTm_dP'] = chem.melting_curve(model)

    #Check if top-down freezing
    if prm.iron_snow:
        idx = core._snow_idx
        if idx > 1:
            flag = snow.check_top_down_freezing(profiles['r'][:idx], profiles['T'][:idx], profiles['Tm'][:idx])

            if not flag:
                model.critical_failure = True
                model.critical_failure_reason = 'Not top down freezing'
                logger.critical(f'it: {model.it}. Intermediate freezing occuring, not exclusively top-down from CMB!')


leeds.update = update


def set_Q_rs(model):
    '''
    Sets the heat flow at the top of the convecting region (rs)

    Parameters
    ----------
    model : ThermalModel Class
        Main model

    Returns
    -------
    Float
        The heat flow at rs
    '''
    prm  = model.parameters
    core = model.core
    sl   = model.stable_layer
    rs_idx = core._rs_idx
    snow_idx = core._snow_idx
    k = core.profiles['k']


    if prm.stable_layer:
        #If stable layer is thicker than snow layer
        if core.rs < core.r_snow and core.rs < prm.r_cmb:
            Q_rs = -4 * np.pi * core.rs**2 * k[rs_idx] * sl.T_grad_s

        #If snow layer is thicker than stable layer
        elif core.r_snow < core.rs:
            Q_rs = -4 * np.pi * core.r_snow**2 * k[snow_idx] * sl.T_grad_s

        else:
            Q_rs = model.core.Q_cmb

    else:
        Q_rs = model.core.Q_cmb

    return Q_rs

leeds.set_Q_rs = set_Q_rs

