#Description of this model3
Description = 'Core model based on Greenwood et al. (2021) with iron snow parameterisation of Davies and Pommier (2018).'

#List individually confirmed compatibility with other methods
compatibility = {'stable_layer': ['leeds_thermal'],
                 'mantle': ['knibbe_westrenen18']}

import numpy as np
import copy

from thermal_history.utils.optimised_funcs import linspace, polyval, trapezoid
from .routines import profiles as prof
from .routines import energy as en
from .routines import chemistry as chem
from .routines import snow as snow


import logging
logger = logging.getLogger(__name__)


#Required functions: setup(), evolution(), update()

def setup(model):
    '''
    Sets initial conditions from the input parameters and performs any other housekeeping (e.g. ensuring correct temperature for a given ri)

    Parameters
    ----------
    model : ThermalModel Class
        Main model
    '''

    prm = model.parameters
    core = model.core

    #Core density polynomials
    rho_l = prm.core_liquid_density_params
    rho_s = prm.core_solid_density_params

    #Set initial conditions
    core.ri   = copy.deepcopy(prm.ri)
    core.rs  = prm.r_cmb
    if prm.stable_layer:
        core.rs = prm.r_cmb - prm.layer_thickness
    core.r_snow = copy.deepcopy(prm.r_cmb)

    core.T_cmb = copy.deepcopy(prm.T_cmb) #CMB temperature

    #Set core mass/mole fractions of light elements.
    core.conc_l = copy.deepcopy(np.array(prm.conc_l)).astype('float')
    if not prm.mf_l is None:
        logger.warning(' mf_l has been specified in parameters, this will overrule conc_l')
        core.mf_l = copy.deepcopy(np.array(prm.mf_l)).astype('float')
        core.conc_l = chem.mole_frac2mass_conc(core.mf_l, prm.mm)
    else:
        core.mf_l = chem.mass_conc2mole_frac(core.conc_l, prm.mm)


    #Save initial mass to conserve it as composition evolves
    # M_ic, M_conv, Ms, M = prof.mass(rho_s, rho_l, core.ri, core.r_s, prm.r_cmb)
    core.M0 = prof.mass(rho_s, rho_l, core.ri, core.rs, prm.r_cmb)[-1]
    core.dM = 0 #change in mass from initial state

    #Save initial light element concentration (used by snow zone evolution).
    core.initial_conc_l = core.conc_l[0]

    #Set initial radial profiles
    prof.basic_profiles(model, setup=True)

    #Initial uniform composition profile
    core.profiles['conc_l'] = np.full(prm.n_profiles, core.conc_l[0])

    core.Tcen = core.profiles['Ta'][0] #Temperature at r=0

    #Calclulate melting curve
    core.profiles['Tm_fe'], core.profiles['Tm'], core.profiles['dTm_dP'] = chem.melting_curve(model)
    core.dTm = core.profiles['Tm'][core._ri_idx] -  core.profiles['Tm_fe'][core._ri_idx]


    #Set CMB compositional values for tracking
    core.c_cmb = copy.deepcopy(core.conc_l)

    #If an inner core is present, set the temperature of the core to be consistent with melting curve (profiles are reset)
    if core.ri > 0 and not (core.profiles['Tm'][core._ri_idx] == core.profiles['Ta'][core._ri_idx]):
        logger.warning(' T_cmb not consistent with given ri, resetting the temperature of the core.')
        chem.calibrate_temperature(model)


    #Calculate temperature profile and profiles that do/may depend on temperature:
    prof.temp_dependent_profiles(model, setup=True)

    #Track upper temperature of adiabatic region (can deviate from T_cmb if stable layer is present)
    core.T_upper = core.profiles['T'][core._rs_idx]


def evolve(model):
    '''
    Main evolution fuction for development method. Calculates rates of change for all aspects of core evolution.

    Parameters
    ----------
    model : ThermalModel Class
        Main model

    Raises
    ------
    ValueError
        If model.mantle.Q_cmb is not set either by a mantle method or manually by the user.
    ValueError
        if the cooling rate, dT_dt, is NaN, a common symptom of division by zero in one of the energy terms.
    '''

    ###############################################################################
    #Call relevant functions for the core evolution
    ###############################################################################
    prm = model.parameters
    core = model.core
    mantle = model.mantle

    #Core density polynomials
    rho_l = prm.core_liquid_density_params
    rho_s = prm.core_solid_density_params

    #Set 0th order polynomial coefficient for outer core density (changes to conserve mass as ri grows)
    if hasattr(core,'o_rho_0'):
        rho_l[0] = core.o_rho_0


    #Mass in regions
    # M_ic, M_conv, Ms, M = prof.mass(rho_s, rho_l, core.ri, core.r_s, prm.r_cmb)
    M_ic, M_conv, Ms, M = prof.mass(rho_s, rho_l, core.ri, core.rs, prm.r_cmb)



    #Get profiles

    r, Ta, Ta_grad, g, psi, P, rho, k, cp, alpha = core.profiles['r'],\
                                                   core.profiles['Ta'],\
                                                   core.profiles['dTa_dr'],\
                                                   core.profiles['g'],\
                                                   core.profiles['psi'],\
                                                   core.profiles['P'],\
                                                   core.profiles['rho'],\
                                                   core.profiles['k'],\
                                                   core.profiles['cp'],\
                                                   core.profiles['alpha']

    Tm_fe, Tm = core.profiles['Tm_fe'], core.profiles['Tm']

    #Whole core T profile (inc stable layer)
    T = core.profiles['T']

    #Indicies for key radii in radial profiles
    ri_idx   = core._ri_idx
    snow_idx = core._snow_idx
    rs_idx   = core._rs_idx


    #Check CMB heat flow has been specified
    if (not prm.mantle and not hasattr(mantle, 'Q_cmb')) and not prm.set_cooling_rate:
        message = 'CMB heat flow (ThermalModel.mantle.Q_cmb) must be set when no mantle model is included'
        logger.critical(message)
        raise ValueError(message)

    core.Q_cmb = mantle.Q_cmb   #Heat flow at CMB
    core.Q_rs = set_Q_rs(model) #Heat flow at base of stable layer

    #Latent heat in J/kg
    if prm.core_latent_heat is None:
        L = core.profiles['dS'] * prm.Na *1000/prm.mm[0] * Tm   #Latent heat based on dS for iron
    else:
        L = np.ones(Tm.size)*prm.core_latent_heat               #Fixed value for latent heat


    #Snow zone evolution
    Q_snow_tilde, E_snow_tilde = 0, 0
    if prm.iron_snow:
        snow_dict = snow_evolution(model) #Evaluate snow evolution

        #Iterate over return values from snow_dict
        for key, value in snow_dict.items():
            #Add on snow energies and entropies
            if key[0] == 'Q':
                Q_snow_tilde += value
            elif key[0] == 'E':
                E_snow_tilde += value


    #############################
    #Calculate Cc and Cr for inner core
    if core.ri>0 and core.ri<prm.r_cmb:

        #Temperature gradients
        dTm_dr = (Tm[ri_idx+1]-Tm[ri_idx])/(r[ri_idx+1]-r[ri_idx]) #Forward difference Tm gradient
        dTa_dr = Ta_grad[ri_idx]  #Adiabatic gradient

        dTm_dr = core.profiles['dTm_dP'][ri_idx-1]*(-rho[ri_idx-1]*g[ri_idx-1])
        #breakpoint()

        #Change in composition normalised to cooling rate
        Cc = 4*np.pi*r[ri_idx]**2*rho[ri_idx]*(core.conc_l-core.conc_s)/M_conv

        #Change in melting temp with mole fraction. Assumes AL melting temperature parameterisation.
        if core.profiles['dS'][ri_idx] == 0:
            dTm_dmf = 0
        elif type(prm.core_melting_params[0]) == str and not prm.core_melting_params[0] == 'AL':
            logger.warning('Change in melting temp with mole fraction is only implemented with \'AL\' melting curve parameterisation. Defaulting to 0.')
            dTm_dmf = 0
        else:
            dTm_dmf = -Tm_fe[ri_idx]*prm.kb/core.profiles['dS'][ri_idx]

        #Calculate forward diff gradient in mole fraction with mass concentration for each LE
        dmf_dc = np.zeros(core.conc_l.size)
        for i in range(core.conc_l.size):
            _mf = copy.copy(core.mf_l)
            _mf[i] += 0.001
            _mass = chem.mole_frac2mass_conc(_mf, prm.mm)

            dmf_dc[i] = (_mf[i]-core.mf_l[i])/(_mass[i]-core.conc_l[i])
        #############

        #Change in melting temp with mass fraction
        dTm_dc = dTm_dmf * dmf_dc

        #Change in melting temp with inner core growth
        dTm_dri = np.sum(dTm_dc*Cc)
        dTm_dri = 0 #MAKE ZERO UNTIL PROPERLY TESTED

        #ICB velocity normalised to cooling rate
        Cr = (1/(dTm_dr-dTa_dr + dTm_dri))*(Ta[ri_idx]/core.Tcen)


    else:
        Cc = np.zeros(len(core.conc_l))
        Cr = 0

    #############################
    # Error can be caused by negative gradient in adiabat near center. Cr should always be <= 0.
    # Ta gradient is required to be = 0 @ r=0 in prof.basic_profiles(model, setup=True) but incase
    # the user specifies their own polynomials that create this issue, limit Cr <=0. A jump in the
    # ICB from 0 to some larger radius can be created in this case.
    if Cr > 0:
        logger.warning(f'Cr is negative (={Cr}), setting to 0.')
        Cr = 0
    #############################

    #If snow zone is enabled, use snow defined Cc, Cr.
    if prm.iron_snow:
        Cc = snow_dict['Cc']
        Cr = snow_dict['Cr']


    #############################
    #Calculate normalised energies/entropies and cooling of adiabatic region

    #Latent heat
    Ql_tilda, El_tilda = en.latent_heat(core.ri, rho, Ta, Cr, L, ri_idx)

    #Gravtiational energy
    Qg_tilda, Eg_tilda = 0, 0
    if not prm.iron_snow: #Handled in snow evolution
        Qg_tilda, Eg_tilda = en.gravitational(r, Ta, rho, psi, Cr, Cc, M_conv, prm.alpha_c, ri_idx, rs_idx, Tcmb=T[-1])

    #MgO Precipitation
    if Ta[-1] < prm.precip_temp:
        Qg_mgo_tilda, Eg_mgo_tilda = en.gravitational_precip(r, Ta, rho, psi, 0, prm.Cm_mgo, prm.alpha_c_mgo)
    else:
        Qg_mgo_tilda, Eg_mgo_tilda = 0, 0

    #Secular cooling
    Qs_tilda, Es_tilda = en.secular_cool(r, rho, Ta, cp, rs_idx, Tcmb=T[-1])

    #HoR #Not tested.
    # Eh_tilda = en.heat_of_reaction(rho[ri_idx:],prm.dmu_dT_poly,Cc,Cr,r[ri_idx:],mm) #Not yet implimented
    Eh_tilda = 0

    #Radiogenic heating. Not been tested in a long time! ATM assumes no stable layer
    if prm.core_h0 > 0:
        Qr, Er = en.radiogenic_heating((model.time-4.5e9*prm.ys), r, rho, Ta, core.M0, prm.core_h0, prm.half_life)
    else:
        Qr, Er = 0, 0
   


    #Entropy of conduction down the adiabat
    Ek = en.cond_entropy(Ta, Ta_grad, k, r, rs_idx)

    #############################
    #Entropy of barodiffusion (across entire liquid region of core inc stable layer.)
    Ea = 0
    alpha_D = 0
    if prm.include_baro_diffusion:
        conc_l = core.conc_l.copy()
        mf_l = chem.mass_conc2mole_frac(conc_l, prm.mm)
        alpha_D = chem.baro_coeff(mf_l, prm.mm, np.mean(T), np.mean(rho), prm.lambda_liq, prm.diffusivity_c)
        for i in range(core.conc_l.size):
            i_flux = en.mass_flux(rho[ri_idx:], prm.alpha_c[i], alpha_D[i], prm.diffusivity_c[i], 0, g[ri_idx:])
            Ea += en.mass_diffusion(r[ri_idx:], i_flux, alpha_D[i], T[ri_idx:])

    #############################
    core.alpha_D = alpha_D

    if core.rs == 0: #No convecting region
        dT_dt = core.dT_dt #core.dT_dt is required to be set in stable layer method when layer covers entire core.
    else:
        dT_dt = (core.Q_rs-Qr)/(Qs_tilda+Ql_tilda+Qg_tilda+Qg_mgo_tilda+Q_snow_tilde)

        #If cooling rate is enforced, overrule dT_dt and reset Q_cmb.
        if prm.set_cooling_rate and not prm.stable_layer:
            dT_dt = prm.core_dTa_dt
            core.Q_rs = (Qs_tilda+Ql_tilda+Qg_tilda+Qg_mgo_tilda+Q_snow_tilde)*dT_dt + Qr
            core.Q_cmb = core.Q_rs

    if np.isnan(dT_dt):
        message = 'Cooling rate is NaN, possibly resulting from division by zero.'
        logger.critical(message)
        raise ValueError(message)

    Es = Es_tilda*dT_dt
    El = El_tilda*dT_dt
    Eg = Eg_tilda*dT_dt
    Eg_mgo = Eg_mgo_tilda*dT_dt
    Eh = Eh_tilda*dT_dt

    #Calculate ohmic dissipation
    Ej = Es+El+Eg+Eh+Er+Eg_mgo - Ek - Ea + E_snow_tilde*dT_dt

    #Check to see if limits are placed on Ej. Can only be used if not solving for a mantle or stable layer
    #since Ej will determine Q_cmb.

    Ej_flag = False
    if not prm.mantle and not prm.stable_layer and model.dt<0:
        #Specific Value of Ej specified
        if not prm.Ej_fixed is None:
            Ej = prm.Ej_fixed
            Ej_flag = True
        #Lower bound on Ej specified
        elif not prm.Ej_lower_bound is None:
            Ej = prm.Ej_lower_bound
            Ej_flag = True
        #Value of Ej prior to inner core nucleation set to value immedietly after ICN
        elif prm.Ej_fixed_pre_ic and core.ri == 0:
            Ej = core._Ej_last
            Ej_flag = True

    #If needed, recalculate cooling rate and CMB heat flow based on required value for Ej
    if Ej_flag:
        dT_dt = (Ej + Ek + Ea - Er) / (Es_tilda + El_tilda + Eg_tilda + Eg_mgo_tilda + Eh_tilda + E_snow_tilde)
        core.Q_rs = dT_dt*(Qs_tilda + Ql_tilda + Qg_tilda + Qg_mgo_tilda + Q_snow_tilde) + Qr
        core.Q_cmb = core.Q_rs

        Es = Es_tilda*dT_dt
        El = El_tilda*dT_dt
        Eg = Eg_tilda*dT_dt
        Eg_mgo = Eg_mgo_tilda*dT_dt
        Eh = Eh_tilda*dT_dt


    #Correct outer core density profile to keep mass of core constant.
    #CHECK THIS IS CONSISTENT WITH STABLE LAYER

    #New mass entering the core from the mantle
    core.dM += mantle.cmb_mass_flux

    if core.ri<prm.r_cmb:
        o_rho = prof.mass_correct(rho_s, rho_l, core.ri, core.M0+core.dM, prm.r_cmb)
        core.o_rho_0 = o_rho[0]
    else:
        core.o_rho_0 = rho_l[0]


    #Save current value of Ej if prm.Ej_fixed_pre_ic=True
    if prm.Ej_fixed_pre_ic:
        core._Ej_last = Ej


    #Save values to be tracked in model.save_dict
    core.Qa    = core.profiles['Qa'][-1]
    core.Qa_rs = core.profiles['Qa'][rs_idx]

    core.ADR = core.Q_cmb/core.profiles['Qa'][-1]
    if core.profiles['Qa'][rs_idx] > 0: #Avoid division by 0
        core.ADR_s = core.Q_rs/core.profiles['Qa'][rs_idx]
    else:
        core.ADR_s = 0

    core.Ql = Ql_tilda*dT_dt   #FIX LATENT HEAT T not Ta
    core.El = El

    core.Qg = Qg_tilda*dT_dt
    core.Eg = Eg

    core.Qg_mgo = Qg_mgo_tilda*dT_dt
    core.Eg_mgo = Eg_mgo

    core.Qr = Qr
    core.Er = Er

    core.Qs = Qs_tilda*dT_dt
    core.Es = Es

    core.Ek = Ek
    core.Eh = Eh
    core.Ej = Ej
    core.Ea = Ea


    #Snow variables
    if prm.iron_snow:
        for key, value in snow_dict.items():
            if key[0] == 'Q' or key[0] == 'E':
                setattr(core, key, value*dT_dt)

        core.dr_snow_dt = Cr*dT_dt
        core.Cp = snow_dict['Cp']
        Cp = snow_dict['Cp']
        core.L_r_snow = L[snow_idx]
    else:
        Cp = 0

    core.dT_dt  = dT_dt
    core.dri_dt = Cr*dT_dt
    core.dc_dt  = (Cp+ Cc*Cr)*dT_dt

    core.Cc = Cc
    core.Cr = Cr

    core.o_rho_0 = rho_l[0]
    core.mf_l = chem.mass_conc2mole_frac(core.conc_l, prm.mm)
    core.mf_s = chem.mass_conc2mole_frac(core.conc_s, prm.mm)

    #Mass conservation check
    if core.ri>0:
        core.mass_le_solid = 4*np.pi*trapezoid(r[:ri_idx], r[:ri_idx]**2*rho[:ri_idx]*np.sum(core.conc_s))[-1]
    else:
        core.mass_le_solid = 0

    if core.r_snow > 0:
        if core.r_snow == prm.r_cmb:
            radius = r
            density = rho
        else:
            radius = r[:snow_idx]
            density = rho[:snow_idx]
        core.mass_le_liq  = 4*np.pi*trapezoid(radius, radius**2*density*np.sum(core.conc_l))[-1]
    else:
        core.mass_le_liq = 0

    if core.r_snow < prm.r_cmb:
        core.mass_le_snow = 4*np.pi*trapezoid(r[snow_idx:], r[snow_idx:]**2*rho[snow_idx:]*core.profiles['conc_l'][snow_idx:])[-1]
    else:
        core.mass_le_snow = 0

    core.dTm = Tm[ri_idx] - Tm_fe[ri_idx]
    core.T_upper = Ta[rs_idx]
    core.T_cmb = T[-1]
    core.L_ri = L[ri_idx]



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

        core.r_snow = snow.snow_radius(profiles['r'], profiles['T'], profiles['Tm'])

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

    #Melting temperature constrained to temperature profile with iron snow
    if prm.iron_snow:
    
        flag = snow.check_top_down_freezing(profiles['r'], profiles['T'], profiles['Tm'])
        if not flag:
            model.critical_failure = True
            model.critical_failure_reason = 'Not top down freezing'
            logger.critical(f'it: {model.it}. Intermediate freezing occuring, not exclusively top-down from CMB!')

        elif core.r_snow < prm.r_cmb:
            profiles['Tm'][core._snow_idx:] = core.profiles['T'][core._snow_idx:]


#Dictionary of required parameters by this model. 'Name': 'Description'
required_params = {'T_cmb': 'Initial temperature of the CMB. Can specify Tcen (temperature at the center of the core) instead. Float.',
                   'conc_l': 'Initial light element mass fraction of the core. Preserve order consistency with all other light element property lists. List(float).',
                   'core_solid_density_params': 'Inner core density polynomials (radial). List(float)',
                   'core_liquid_density_params': 'Outer core density polynomials (radial). List(float)',
                   'ri':    'Initial inner core radius. Float',
                   'r_cmb': 'CMB radius. Float',
                   'core_alpha_T_params': 'Core thermal expansivity pressure polynomial coefficients. List(Float)',
                   'core_cp_params': 'Core specific heat capacity pressure polynomial coefficients. List(Float)',
                   'core_conductivity_params': 'List, first element is a string followed by the core thermal conductivity polynomial coefficients. The string (either \'r\'/\'T\'/\'P\') indicates if the polynomials are radial, temperature, or pressure coefficients. Can instead provide a single number which will be taken as the constant conductivity value. List(string, Float, Float...) or List(Float).',
                   'core_melting_params': 'List, first element is a string followed by the constants used for the melting temperature parameterisation. The string element indicates the parameterisation to be used, if no string is given, \'AL\' (method of Alfe et al. 2002) is assumed. See melting_curve() in core_models/leeds/routines/chemistry.py for possible options.',
                   'entropy_melting_params': 'Change of entropy on freezing pressure polynomial coefficients. List(Float).',
                   'mm': 'Molar masses (g/mol) of Iron followed by alloying light elements. List(float)',
                   'alpha_c': 'Chemical expansivity of alloying light elements. List(float)',
                   'diffusivity_c': 'Chemical diffusivities of alloying light elements.  List(float)',
                   'use_partition_coeff': 'Boolean dictating if constant partition coefficients should be used (True) or if partitioning based on chemical equilibrium should be used (False)',
                   'core_h0': 'Present day radiogenic heating per unit mass in the core',
                   'half_life': 'Half life of radioactive isotope in the core.',
                   'partition_coeff': 'Partition coefficients (x_i/x_o) for mass fraction for each light element. Defined as the ratio of inner core concentration over outer core. List(float)',
                   'lambda_sol': 'Linear corrections to chemical potentials (solid phase). List(float)',
                   'lambda_liq': 'Linear corrections to chemical potentials (liquid phase). List(float)',
                   'dmu': 'Change in chemical potential between solid and liquid. List(float)',
                   'n_profiles': 'Number of nodes used in radial profiles for core properties. Float',
                   'P_cmb': 'CMB pressure. Float',
                   'precip_temp': 'Temperature below which precipitation of MgO begins. Float',
                   'Cm_mgo': 'mass fraction of MgO in the core',
                   'alpha_c_mgo': 'MgO chemical expansivity'}

#Optional parameters that the user may omit, in which event are set to the specified default values in this dictionary.
#Use the format {'Name': ('Description', default_value)}
optional_params = {'include_baro_diffusion': ('(Default: True) Boolean, if True, barodiffusion will be included for both the entropy budget and also chemical gradients in any chemcially stratified layer.', True),
                  'core_adiabat_params': ('(Default: None) None or List(float). Radial polynomial coefficients for the core adiabat normalised to T(r=0). If left as None, these will be calculated by fitting to a numerically integrated adiabatic gradient using given alpha/cp/g profiles.', None),
                  'set_cooling_rate': ('(Default: False) Boolean, if True, the core cooling rate is set to the parameter core_dTa_dt and the CMB heat flow is overwritten with the required value. Only used if no stable layer is included', False),
                  'core_dTa_dt': ('(Default: 0) Float, the cooling rate of the core if set_cooling_rate is True', 0),
                  'iron_snow': ('(Default False) Boolean, if True, an top down core crystallisation rather than bottom up is assumed', False),
                  'Ej_fixed': ('(Default None) None or Float, if set to a value, Ej will be fixed to this value at all times. Q_cmb will be altered to ensure this condition', None),
                  'Ej_lower_bound': ('(Default None) None or Float, if set to a value, Ej will be limited to this value if it were to fall below it. Q_cmb will be altered to ensure this condition', None),
                  'Ej_fixed_pre_ic': ('(Default False). Boolean, when dt<0 and if true, Ej is fixed prior to ICN to the value immedietly post ICN. Q_cmb will be altered to ensure this condition.', False),
                  'core_latent_heat': ('(Default None). Latent heat is calculated using the change of entropy on melting unless this is set to a fixed number (J/kg). None or Float', None),
                  'mf_l': ('(Default None) Starting mole fraction of light element in the core. If specified it will overrule conc_l.', None),
                  'Tcen': ('(Default None) Initial temperature at the center of the core. If specified it will overrule T_cmb (but may still be overrulled by contraint of T=Tm at ri.)', None),
                  'use_new_Cr': ('(Default False) Use new Cr factor that takes into account depression of Tm with changing LE concentration', False)
}


def progress(model):
    '''
    Return text for this class when progress at each iteration is printed to STDOUT

    Returns
    -------
    String
        String of text to be printed to screen.
    '''  

    core = model.core
    prm = model.parameters

    v = (core.Tcen, core.Q_rs/1e12, core.ADR_s)

    text = f'    Tcen: {float(v[0]):.2f} ˚K    Q(rs): {float(v[1]):.2f} TW    ADR(rs): {float(v[2]):.2f}'

    if model.parameters.iron_snow:
        text += f'    snow depth: {(prm.r_cmb-core.r_snow)/1000:.2f} km'
    else:
        text += f'    ri: {core.ri/1000:.2f} km'  

    return text

#Extra functions used by the main functions 'setup', 'evolve', and 'update'.


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
    k = core.profiles['k']

    if prm.stable_layer and core.rs < prm.r_cmb:
        Q_rs = -4 * np.pi * core.rs**2 * k[rs_idx] * sl.T_grad_s
    else:
        Q_rs = model.core.Q_cmb

    return Q_rs

def snow_evolution(model):
    '''
    Sub-evolution fuction specifically for an iron snow-zone based on Davies and Pommier (2018). Called from evolution() if required (prm.iron_snow=True)

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

    #Initialise snow parameters Cl, concl_l and phi_snow profiles
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
    ri_idx = core._ri_idx
    snow_idx = core._snow_idx
    rs_idx = core._rs_idx

    if core.r_snow < prm.r_cmb and not core.r_snow==0:
        assert hasattr(core, 'dT_dt'), 'No cooling rate has been calculated yet on iteration {}'.format(model.it)

        #Concentration of slurry to depress Tm to adiabatic temperature
        conc_l_profile[snow_idx:] = snow.snow_composition(T, Tm_fe, core.initial_conc_l, snow_idx, core.conc_l[0], P, prm.core_melting_params)
        core.profiles.update({'conc_l': conc_l_profile})

        #Move melting temperature accordingly and update latent heat
        # Tm = Tm_fe*(1 - conc_l_profile/(1+0.106))
        Tm_fe, Tm, dTm_dP = chem.melting_curve(model) #More general function

        #Melting temperature constrained to temperature profile with iron snow
        Tm[core._snow_idx:] = core.profiles['T'][core._snow_idx:]

        #Latent heat in J/kg
        if prm.core_latent_heat is None:
            L = core.profiles['dS'] * prm.Na *1000/prm.mm[0] * Tm   #Latent heat based on dS for iron
        else:
            L = np.ones(Tm.size)*prm.core_latent_heat               #Fixed value for latent heat
        # L = prof.entropy_melting(P, prm.entropy_melting_params) * prm.Na *1000/prm.mm[0] * Tm

        M_liquid = prof.mass(rho_s, rho_l, 0, core.r_snow, core.r_snow)[-1] #Mass of liquid region

        #Cc factor, analagous to inner core, relates changing mass fraction in liquid region to movement of the snow zone radius.
        Cc_snow = 4*np.pi*r[snow_idx]**2*rho[snow_idx]*(conc_l_profile[snow_idx]-core.conc_l[0])/M_liquid

        #Cl factor, normalising changes in slurry mass fraction to changes in temperature (DP eq. 22)
        # although it could be calculated by derivative of WN melting temp (since it equals adiabat in snow zone)
        #Cl[snow_idx:] = -1/Tm_fe[snow_idx:]

        dmu_dc = (prm.kb*T/conc_l_profile)*prm.Na*(1000/prm.mm[1])
        Cl = snow.Cl_calc(phi_snow, L, T, conc_l_profile, dmu_dc, snow_idx)

        #Solid fraction in slurry formed in one timestep. Needed to calculate Ql_snow, total latent heat of moving boundary.
        phi_snow = (Cl/conc_l_profile) * (T/core.Tcen) * core.dT_dt*model.dt

        #Normalised snow zone radius velocity to cooling rate, same as done for inner core
        dTm_dr = -rho[snow_idx]*g[snow_idx]*dTm_dP[snow_idx]
        # dTm_dr = (Tm[snow_idx-2]-Tm[snow_idx-1])/(r[snow_idx-2]-r[snow_idx-1])
        dT_dr  = (T[snow_idx-2]-T[snow_idx-1])/  (r[snow_idx-2]-r[snow_idx-1])

        dTm_dc = -Tm_fe[snow_idx]/(1+core.initial_conc_l) #Assumes WN melting curve
        if prm.core_melting_params[0] == 'RI': #Use Rivoldini data
            dTm_dc = chem.riv_dTm_dc(P[snow_idx], core.conc_l[0], prm.core_melting_params[1:].astype('float64'))


        if prm.use_new_Cr:
            Cr_snow = (1/(dTm_dr + dT_dr + dTm_dc*Cc_snow))*(T[snow_idx]/core.Tcen)  #Tested this enough to determine it is required.
        else:
            Cr_snow = (1/(dTm_dr - dT_dr))*(T[snow_idx]/core.Tcen)


        #Normalised latent heat from moving snow zone radius (DP18 Qlb eq 14). Analagous to inner core growth.
        Ql_snow_tilde, El_snow_tilde = en.latent_heat(core.r_snow, rho*phi_snow[snow_idx], T, Cr_snow, L, snow_idx) #Multiply density by phi for just mass of freezing solid.

        #Latent heat relase and absorbtion upon freezing and melting. Ql_melting is given negative L to indicate heat absorption.
        #(DP18 Qls and Qll)
        Ql_freezing_tilde, El_freezing_tilde = snow.latent_snow(r, rho, L,            Cl, conc_l_profile, T, snow_idx)
        Ql_melting_tilde,  El_melting_tilde  = snow.latent_snow(r, rho, -L[snow_idx], Cl, conc_l_profile, T, snow_idx)


        #Cp factor, normalises changing mass fraction in liquid region to changes in slurry mass fraction (which in turn uses Cl to relate to cooling)
        #Uses DP18 eqn 20
        Cp_snow = snow.Cp_factor(r, rho, Cl, T, M_liquid, snow_idx)

        #Gravitational energy release from removal of solid by melting, then increasing density of the liquid region (DP18 eqn 16/17)
        Qg_freezing_tilde, Eg_freezing_tilde = snow.gravitational_freezing(r, rho, psi, prm.alpha_c[0], Cl, T, snow_idx)
        Qg_melting_tilde,  Eg_melting_tilde  = snow.gravitational_melting(r, rho,  psi, prm.alpha_c[0], Cp_snow, Cc_snow, Cr_snow, T[-1], snow_idx)


        #I showed in the notes on the gravitational energy than this term should include the jump in conc_l
        #not just conc_l as in DP18.
        Qg_snow_tilde = -4*np.pi*r[snow_idx]**2*rho[snow_idx]*psi[snow_idx]*(conc_l_profile[snow_idx]-conc_l_profile[snow_idx-1])*Cr_snow
        Eg_snow_tilde = Qg_snow_tilde/T[-1]

        #Set for the correct size (total size of number of LE) to avoid inconsistencies down the line.
        Cc_snow = np.ones(core.conc_l.size)*Cc_snow

        # Ql_snow_tilde, Ql_freezing_tilde, Ql_melting_tilde, Qg_freezing_tilde, Qg_melting_tilde, Qg_snow_tilde = 0,0,0,0,0,0

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
