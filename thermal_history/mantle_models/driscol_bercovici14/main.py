Description = 'Mantle model for convection with plate tectonics. \
No heat flow modified by volcanism is used here.'

import numpy as np
import logging
logger = logging.getLogger(__name__)

def setup(model):

    #DB14 method sets the viscosity reference of the mantle to fit present day surface heat flow.

    mantle = model.mantle
    prm = model.parameters

    mantle.Tm   = prm.Tm      #average mantle temperature
    r_surf  = prm.r_surf   #upper radius
    r_cmb  = prm.r_cmb    #lower radius

    T_surf   = prm.T_surf                  #Surface temperature
    alpha    = prm.mantle_alpha_T          #volumetric expansion
    g        = prm.g                       #gravity
    kappa    = prm.mantle_diffusivity_T    #thermal diffusivity
    cp       = prm.mantle_cp               #specific heat capacity
    Rac      = prm.Rac                     #critical Rayleigh Number (upper)
    Rg       = prm.Rg                      #Gas constant
    Av       = prm.activation_energy       #Activation energy

    if prm.basal_magma_ocean:
        r_lower = prm.r_bmo
    else:
        r_lower = r_cmb

    D = r_surf - r_lower      #Mantle Thickness
    r_av = (r_surf+r_lower)/2

    Ta_upper = mantle.Tm*np.exp(-(r_surf-r_av)*alpha*g/cp)

    if not prm.Q_surface == None:
        #On first iteration, set viscosity_reference reference value to match Q_surf
        Q_surface = prm.Q_surface
        delta_upper_calibrated = 4 * np.pi * prm.r_surf**2 * prm.mantle_k_upper * ((Ta_upper-T_surf)/Q_surface)
        nu_upper_calibrated = (delta_upper_calibrated/D)**3*(D**3)*alpha*g*(Ta_upper-T_surf)/(kappa*Rac)
        viscosity_reference = nu_upper_calibrated*10/np.exp(Av/(Rg*mantle.Tm))
        prm.viscosity_ref = viscosity_reference

    if prm.basal_magma_ocean:
        mantle.Tl    = prm.bmo_liquidus
        mantle.r_bmo = prm.r_bmo

        if not prm.core:
            raise ValueError('Including a BMO requires a core model be used as well')


def evolve(model):

    #Read in variables from model and parameters file

    mantle = model.mantle
    core = model.core

    prm = model.parameters

    Tm   = mantle.Tm      #average mantle temperature
    time = model.time

    r_surf  = prm.r_surf   #upper radius
    r_cmb  = prm.r_cmb    #lower radius

    k_upper  = prm.mantle_k_upper  #thermal conductivity in upper mantle
    k_lower  = prm.mantle_k_lower  #thermal conductivity in lower mantle

    T_surf    = prm.T_surf                  #Surface temperature
    Qr0      = prm.Qr0                     #present day radiogenic heating
    Trad     = prm.Trad                    #radiogenic heating decay
    mass     = prm.mantle_mass             #mantle mass
    alpha    = prm.mantle_alpha_T          #volumetric expansion
    g        = prm.g                       #gravity
    kappa    = prm.mantle_diffusivity_T    #thermal diffusivity
    cp       = prm.mantle_cp               #specific heat capacity
    Rac      = prm.Rac                     #critical Rayleigh Number (upper)
    Rg       = prm.Rg                      #Gas constant
    Av       = prm.activation_energy       #Activation energy
    f_visco  = prm.f_visco                 #Ratio of lower to upper viscosity

    #Make sure CMB temperature is specified if no core model is included.
    if not prm.core and not hasattr(core, 'T_cmb'):
        message = 'Core CMB temperature (THERMAL_MODEL.core.T_cmb) must be set when no core model is included'
        logger.critical(message)
        raise ValueError(message)

    T_cmb = core.T_cmb

    #Driscoll and Bercovici Method
    if prm.basal_magma_ocean:
        r_lower = mantle.r_bmo
    else:
        r_lower = r_cmb

    r_av = (r_surf+r_lower)/2

    Ta_upper = Tm*np.exp(-(r_surf-r_av)*alpha*g/cp)
    Ta_lower = Tm*np.exp(-(r_lower-r_av)*alpha*g/cp)

    assert Ta_upper >= T_surf, 'Mantle colder than the surface'

    D = r_surf - r_lower   #Mantle Thickness

    viscosity_reference = prm.viscosity_ref

    nu_upper = viscosity_reference*np.exp(Av/(Rg*Tm))/10
    nu_lower = nu_upper*f_visco

    dT_upper = np.abs(Ta_upper-T_surf)
    delta_upper = D*((kappa*nu_upper*Rac)/((D**3)*alpha*g*dT_upper))**(1/3)

    dT_lower = np.abs(T_cmb-Ta_lower)
    delta_lower = D*((kappa*nu_lower*Rac)/((D**3)*alpha*g*dT_lower))**(1/3)

    mantle.dT_upper = Ta_upper-T_surf
    mantle.dT_lower = T_cmb-Ta_lower

    mantle.delta_upper = delta_upper
    mantle.delta_lower = delta_lower

    #Conductive heat flow through boundaries

    Q_surface = 4*np.pi*r_surf**2*k_upper*((Ta_upper-T_surf)/delta_upper)

    #heat produced via radiogenic production
    t = 4.5e9*prm.ys - time
    Qr = Qr0*np.exp(t/Trad)

    mantle.Qr = Qr

    mantle.Q_cmb = 4*np.pi*r_cmb**2*k_lower*(((T_cmb-Ta_lower))/delta_lower)
    Q_lower = mantle.Q_cmb


    #BMO routine

    #Initialise flux BC with zero value. First value indicates fixed value(0) or fixed flux(1), then second number is the value of the BC.
    mantle.chemical_bc_cmb = np.array([1,0]) 
    if prm.basal_magma_ocean:
        mantle.Q_bmo = Q_lower

        #Calc bmo evolution
        bmo_evolution(model)

        #Split off mass and Qr in bmo
        Qr      = mantle.Qr * (1-mantle.mass_bmo/mass)
        mass   -= mantle.mass_bmo
        Q_lower = mantle.Q_bmo
    
            
    #Calculating energy from mantle secular cooling
    Qs = Q_surface - Q_lower - Qr

    #Mantle cooling rate
    dTm_dt = -Qs/(mass*cp)

    #Save values to model class

    mantle.nu_upper = nu_upper
    mantle.nu_lower = nu_upper*10

    mantle.Qr = Qr
    mantle.Qs = Qs
    mantle.Q_surface = Q_surface

    mantle.dT_dt = dTm_dt
    mantle.T_cmb  = T_cmb

    mantle.T_upper = Ta_upper
    mantle.T_lower = Ta_lower

    # breakpoint()


def update(model):

    prm = model.parameters
    dt = model.dt
    mantle = model.mantle

    mantle.Tm += mantle.dT_dt * dt

    if prm.basal_magma_ocean:
        mantle.Tl       += mantle.dTl_dt_bmo * dt
        mantle.r_bmo    += mantle.dr_dt_bmo * dt
        mantle.cFeO_bmo += mantle.dc_dt_bmo * dt

        if mantle.r_bmo < prm.r_cmb:
            mantle.r_bmo = prm.r_cmb


def progress(model):
    '''
    Return text for this class when progress at each iteration is printed to STDOUT

    Returns
    -------
    String
        String of text to be printed to screen.
    '''  

    mantle = model.mantle
    prm = model.parameters

    v = (mantle.Tm, mantle.Q_cmb/1e12)

    text = f'    Tm: {float(v[0]):.2f} ˚K    Q_cmb: {float(v[1]):.2f} TW'

    if prm.basal_magma_ocean:
        text += f'    BMO size: {(mantle.r_bmo-prm.r_cmb)/1000: .2f} km'

    return text




#Required parameters. 'Name': 'Description'
required_params = {'T_surf': 'Surface temperature. Float',
                   'mantle_alpha_T': 'Thermal expansivity. Float',
                   'g': 'Gravity, assumed constant throughout mantle. Float',
                   'mantle_diffusivity_T': 'Thermal diffusivity. Float',
                   'mantle_cp': 'Mantle specific heat capacity. Float',
                   'Rac': 'Critical Rayleigh number. Float',
                   'Rg': 'Gas constant (J/K/mol). Float',
                   'activation_energy': 'Activation energy for Arrheniun viscosity. Float',
                   'r_surf': 'Planet radius. Float',
                   'r_cmb': 'CMB radius. Float',
                   'viscosity_ref': 'Reference viscosity of the mantle. If Q_surface is set, viscosity will be adjusted to give the value of Q_surface on the first iteration.',
                   'Tm': 'Initial bulk mantle temperature. Float',
                   'mantle_k_upper': 'Upper mantle thermal conductivity. Float',
                   'mantle_k_lower': 'Lower mantle thermal conductivity. Float',
                   'f_visco': 'Viscosity ratio between upper and lower mantle.',
                   'Qr0': 'Present day radiogenic heating rate. Float.',
                   'Trad': 'Radiogenic half life. Float.',
                   'mantle_mass': 'Mass of mantle. Float'
                   }

optional_params = {'Q_surface': ('(Default: None). Initial conductive heat flow through upper mantle thermal boundary layer (inc. crust), used to define viscosity_ref of the mantle.', None),
                   'basal_magma_ocean': ('(Default: False) Toggle to include a basal magma ocean', False),
                   'r_bmo': ('(Default: None). Initial radius of the top of the basal magma ocean. Float', None),
                   'bmo_liquidus': ('(Default: None). Initial bmo liquidus temperature . Float', None),
                   'FeO_mf_bmo': ('(Default: 0) Initial mole fraction of FeO in bmo. Float', 0.5),
                   'partition_FeO': ('(Default: 0) Partition coefficient between mantle bmo and core', 0),
                   'tres': ('(Default 0) Residence time for chemical boundary layer. Float', 0),
                   'ds_bmo': ('(Default 0) Change in entropy on freezing for the BMO. Float', 0),
                   'dc_bmo': ('(Default 0) Difference in FeO between mantle and bmo', 0),
                   'rho_bmo': ('(Default 5000) BMO density. Float', 5000),
                    }





#BMO evolution routine

from .routines import bmo_functions as bmo
from thermal_history.core_models.leeds.routines.chemistry import mass_conc2mole_frac
from thermal_history.core_models.leeds.routines.energy import integrate

def bmo_evolution(model):

    mantle = model.mantle
    core   = model.core
    prm    = model.parameters

    r_bmo = mantle.r_bmo
    Tl    = mantle.Tl

    cp      = prm.mantle_cp       #Mantle specific heat
    r_cmb   = prm.r_cmb
    k_lower = prm.mantle_k_lower
    rhom    = prm.rho_bmo
    dc      = prm.dc_bmo

    D_c = prm.diffusivity_c[0]  #Core O diffusivity. Assumes O is first element specified in parameters.

    #Bulk BMO mass fractions
    if model.it == 1:
        cMgO_bmo, cFeO_bmo = bmo.mole2massconc_bmo(prm.FeO_mf_bmo)  #Get value based on initial mole fraction
        mantle.cFeO_bmo = cFeO_bmo

    cFeO_bmo = mantle.cFeO_bmo

    diff_ls_bmo        = 2*np.sqrt(D_c * prm.tres)      # Chemical BL thickness in BMO

    #Get molar FeO on mantle side of CMB
    cbarO_core = mass_conc2mole_frac(core.c_cmb, prm.mm)[0]
    cbarFeO_bmo_cmb = cbarO_core/prm.partition_FeO

    #Get mass of FeO on mantle side of CMB and flux into core
    cMgO_bmo_cmb, cFeO_bmo_cmb = bmo.mole2massconc_bmo(cbarFeO_bmo_cmb)
    iFeO_bmo                   = bmo.cmb_flux_mass(cFeO_bmo_cmb,cFeO_bmo,
                                diff_ls_bmo, rhom, D_c)

    mass_bmo = 4.0 * np.pi * rhom * (r_bmo**3 - r_cmb**3 ) / 3


    if r_bmo > r_cmb:
        Qbmo = 4*np.pi*r_bmo**2*k_lower*((mantle.dT_lower)/mantle.delta_lower)

        #Normalised latent heat factor from freezing
        L = Tl*prm.ds_bmo
        Ql_tilde = bmo.latent_labrosse(r_cmb, r_bmo, L, dc)

        #HoR from FeO flux
        Qh   = 4*np.pi*r_cmb**2 * (1.0 - cFeO_bmo) * iFeO_bmo * L / dc #Eq 15
        Qh = 0 #Not using Qh

        #Normalised concentration and bmo radius rates of change.
        dcdt_d   = 4*np.pi*r_cmb**2 * (1.0 - cFeO_bmo) * iFeO_bmo     / mass_bmo                 #Eq 13
        drdt_fac = (r_cmb**2 * (1.0 - cFeO_bmo) * iFeO_bmo)/(r_bmo**2 * dc * rhom) # Eq 14

        Qh_hor = 0.0 #No core side HoR

        Qrbmo = mantle.Qr * mass_bmo/prm.mantle_mass

        #Calculate normalised core cooling to dTl_dt, assumes only secular cooling is enegy in the core.
        r_core, T_core ,rho_core, cp_core = core.profiles['r'], core.profiles['T'], core.profiles['rho'], core.profiles['cp']

        I = 4*np.pi* integrate(r_core, cp_core * rho_core * T_core * r_core**2)

        Qs_core_tilde = -I/core.T_cmb

        dTl_dt    = (Qbmo - Qrbmo + Qh + Qh_hor) / (-mass_bmo*cp + Ql_tilde + Qs_core_tilde)

        Q_cmb    = Qs_core_tilde*dTl_dt
        Qsbmo    = -mass_bmo * cp * dTl_dt
        QLbmo    = Ql_tilde * dTl_dt

        dr_dt    = bmo.dadt_fac_labrosse(r_cmb, r_bmo, dc)*dTl_dt + drdt_fac
        dc_dt_e  = bmo.dcdt_fac_labrosse(r_cmb, r_bmo, dc)*dr_dt
        dc_dt    = dc_dt_e + dcdt_d

        FeO_flux = 4*np.pi*r_cmb**2 * iFeO_bmo            # Total flux into core
        dFeO_dt = -4*np.pi*r_cmb**2 * iFeO_bmo/mass_bmo   # Rate of change of core FeO content

    else:

        iFeO_bmo = 0
        Qbmo = 0
        Q_cmb = mantle.Q_cmb
        Qsbmo, QLbmo, Qrbmo, Qh = 0,0,0,0
        dTl_dt, dc_dt, dr_dt = 0,0,0
        FeO_flux, dFeO_dt = 0,0


    #Save values
    mantle.iFeO_bmo = iFeO_bmo
    mantle.MgO_cmb_bmo, mantle.FeO_cmb_bmo = cMgO_bmo_cmb, cFeO_bmo_cmb

    mantle.mass_bmo = mass_bmo

    mantle.Q_cmb = Q_cmb
    mantle.Q_bmo = Qbmo
    mantle.Qs_bmo = Qsbmo
    mantle.Ql_bmo = QLbmo
    mantle.Qr_bmo = Qrbmo
    mantle.Qh_bmo = Qh

    mantle.dTl_dt_bmo = dTl_dt
    mantle.dc_dt_bmo = dc_dt
    mantle.dr_dt_bmo = dr_dt

    mantle.iFeO_cmb      = FeO_flux      
    mantle.dFeOdt_cmb    = dFeO_dt

    #Set O flux boundary condition at CMB.
    dc_dr = -iFeO_bmo / (model.core.profiles['rho'][-1] * D_c) * 16/72
    mantle.chemical_bc_cmb = np.array([1, dc_dr])
