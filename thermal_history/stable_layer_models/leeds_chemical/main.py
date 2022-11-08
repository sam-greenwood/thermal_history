Description = 'Stable layer model for chemical stratification. Uses the Buffett and Seagle (2010) method \
    for changing the size of the layer.'


#List individually confirmed compatibility with other regions
compatibility = {'core': ['leeds'],
                 'mantle': ['driscoll_bercovici14']}

from ...core_models.leeds.routines import profiles as prof
from ...core_models.leeds.routines import energy as en
from ...core_models.leeds.routines import chemistry as chem

from .routines import functions as func
from .routines.diffusion import diffusion

from thermal_history.utils.optimised_funcs import linspace, trapezoid, polyval

import numpy as np
from scipy.integrate import trapz


def setup(model):
    '''Setup the initial stable layer

    Parameters
    ----------
    model : ThermalModel
        Main model class
    '''

    prm = model.parameters
    sl = model.stable_layer
    core = model.core

    #Initial size of stable layer
    sl.layer_thickness = prm.layer_thickness
    core.rs = prm.r_cmb - prm.layer_thickness

    if prm.layer_thickness==0:
        prm.primordial_layer = 0

    #Commented out because a core model is now required.
    # #If no core model is being used, assume simple case of fully molten core
    # if not prm.core:
    #     core.T_cmb = prm.T_cmb
    #     core.ri = 0
    #     core.r_snow = prm.r_cmb

    #     #Set core profiles
    #     prof.basic_profiles(model)
    #     prof.temp_dependent_profiles(model)
    assert prm.core, "A core model is required for this method!"


    #Initialise profiles dictionaries
    sl.profiles = {}
    sl._next_profiles = {}

    n = prm.max_resolution

    #Create radial profiles
    sl.profiles['r'] = linspace(prm.r_cmb - prm.layer_thickness, prm.r_cmb, n)

    sl.profiles['T'] = prof.adiabat(sl.profiles['r'], core.Tcen, prm.core_adiabat_params)
    sl.T_grad_s      = prof.adiabat_grad(core.rs, core.Tcen, prm.core_adiabat_params)

    sl.profiles['c'] = np.full(n, core.conc_l[0])
    sl.c_grad_s      = 0


def evolve(model):
    '''Evolve the stable layer

    Parameters
    ----------
    model : ThermalModel
        Main model class
    '''

    sl = model.stable_layer
    core = model.core
    prm = model.parameters

    sl.Q_cmb = model.mantle.Q_cmb

    T_grad_cmb = sl.Q_cmb/(-core.profiles['k'][-1]*4*np.pi*prm.r_cmb**2)
    ADR = T_grad_cmb/core.profiles['dTa_dr'][-1]
    sl.ADR = ADR

    #Check if conditions are sub-adiabatic and a layer should begin to grow
    check  = func.pure_chemical_check(model)

    #Initialise energies/entropies associated with this state.
    Qs, Es, Ek, Ej, Ea, Ea_baro = 0, 0, 0, 0, 0, 0

    #initialise BC's
    sl.lb_T, sl.ub_T = 0,0
    sl.T_grad_s = 0

    sl.alpha_D = 0
    sl.baro_grad = 0

    #Get/set radial profiles  
    r_initial = sl.profiles['r']
    T_initial = sl.profiles['T']
    c_initial = sl.profiles['c']

    sl.profiles['rho']    = np.interp(r_initial, core.profiles['r'], core.profiles['rho'])
    sl.profiles['g']      = np.interp(r_initial, core.profiles['r'], core.profiles['g'])
    sl.profiles['cp']     = np.interp(r_initial, core.profiles['r'], core.profiles['cp'])
    sl.profiles['k']      = np.interp(r_initial, core.profiles['r'], core.profiles['k'])

    if np.max(np.abs(np.diff(c_initial))) == 0:
        sl.profiles['c_grad'] = np.zeros(r_initial.size)
    else:
        sl.profiles['c_grad'] = np.gradient(c_initial, r_initial, edge_order=2)


    #Check if CMB boundary conditions would give a layer this time-step.
    if sl.layer_thickness == 0 and not check:

        sl.ds_dt = 0
        sl.dT_dt_s = 0
        sl._next_profiles['r'] = np.ones(prm.max_resolution)*prm.r_cmb
        sl._next_profiles['T'] = prof.adiabat(sl._next_profiles['r'], core.Tcen, prm.core_adiabat_params)
        sl._next_profiles['c'] = np.full(sl._next_profiles['r'].size, core.conc_l[0])


    #Otherwise run the method.
    else:

        #Secular cooling. Adiabat is assumed so takes same form as convecting region
        Qs_tilde, Es_tilde = en.secular_cool(r_initial, sl.profiles['rho'], T_initial, sl.profiles['cp'], idx=r_initial.size) #Integral over stable layer
        Qs = Qs_tilde*core.dT_dt
        Es = Es_tilde*core.dT_dt


        sl.dT_dt_s = core.dT_dt * T_initial[0]/core.Tcen  #Save cooling rate at base of layer

        #Calculate Ek
        if sl.layer_thickness == 0: #Layer is still of 0 thickness
            Ek == 0
        else:
            dT_dr = prof.adiabat_grad(r_initial, core.Tcen, prm.core_adiabat_params)
            Ek = en.cond_entropy(T_initial, dT_dr, sl.profiles['k'], r_initial, r_initial.size)

        #Calculate Ea
        conc_l_av = core.conc_l.copy()
        conc_l_av[0] = np.mean(c_initial) #Average composition of fluid in layer. Assumes other light elements are same as convecting region.

        mf_l_av = chem.mass_conc2mole_frac(conc_l_av, prm.mm) #Average mole fractions

        #Get alpha_D for just first light element and save to model (is used later in pure_chemical_method)
        sl.alpha_D = chem.baro_coeff(mf_l_av[0], prm.mm[:2], np.mean(T_initial), np.mean(sl.profiles['rho']), prm.lambda_liq[0], prm.diffusivity_c[0])

        #Mass flux due to barodiffusion
        i_flux_baro = en.mass_flux(sl.profiles['rho'], prm.alpha_c[0], sl.alpha_D, prm.diffusivity_c[0], 0, sl.profiles['g'])

        #Mass flux down compositional gradient and barodiffusion
        i_flux = en.mass_flux(sl.profiles['rho'], prm.alpha_c[0], sl.alpha_D, prm.diffusivity_c[0], sl.profiles['c_grad'], sl.profiles['g'])

        #Take off barodiffusive component
        if not prm.include_baro_diffusion:
            i_flux     -= i_flux_baro
            i_flux_baro = np.zeros(r_initial.size)

        #Save mass fluxes at CMB and rs.
        sl.mass_flux_cmb,      sl.mass_flux_s      = i_flux[-1],      i_flux[0]
        sl.mass_flux_cmb_baro, sl.mass_flux_s_baro = i_flux_baro[-1], i_flux_baro[0]

        Ea_baro = en.mass_diffusion(r_initial, i_flux_baro, sl.alpha_D, T_initial)
        Ea      = en.mass_diffusion(r_initial, i_flux,      sl.alpha_D, T_initial)

        #Time step model
        pure_chemical_method(model)

        # #New radial profiles
        # r_new = sl._next_profiles['r']
        # T_new = sl._next_profiles['T']
        # c_new = sl._next_profiles['c']

        #Ohmic dissipation
        Ej = Es - Ek - Ea

    #Save energies/entropies to model attributes
    sl.Qs, sl.Es, sl.Ek, sl.Ej, sl.Ea, sl.Ea_baro = Qs, Es, Ek, Ej, Ea, Ea_baro

    #Add onto core results
    core.Qs += Qs
    core.Ej += Ej
    core.Es += Es
    core.Ek += Ek
    core.Ea += Ea

    #Save T/c values at CMB/rs
    sl.T_cmb, sl.T_s= sl.profiles['T'][-1], sl.profiles['T'][0]
    core.T_cmb = sl.T_cmb #core values were not set with knowledge of stable layer

    sl.c_cmb, sl.c_s= sl.profiles['c'][-1], sl.profiles['c'][0]
    core.c_cmb[0] = sl.c_cmb

    #Make sure profiles have size maximum_resolution, so that each iteration has same array size.
    n = sl.profiles['r'].size
    keys = [x for x in sl.profiles.keys() if not x=='r'] #All profiles except r.
    corrected_profiles = sl.profiles.copy()
    if n > prm.max_resolution:
        #Interpolate down
        corrected_profiles['r'] = linspace(sl.profiles['r'][0], sl.profiles['r'][-1], prm.max_resolution)
        for key in keys:
            corrected_profiles[key] = np.interp(corrected_profiles['r'], sl.profiles['r'], sl.profiles['T'])

    elif n < prm.max_resolution:
        #Append cmb values to get right size.
        corrected_profiles['r'] = np.append(sl.profiles['r'], np.ones(prm.max_resolution-n)*sl.profiles['r'][-1])
        for key in keys:
            corrected_profiles[key] = np.append(sl.profiles[key], np.ones(prm.max_resolution-n)*sl.profiles[key][-1])

    sl.profiles = corrected_profiles #Set these corrected arrays to profiles dictionary

    sl.profiles['Ta'] = prof.adiabat(sl.profiles['r'], core.Tcen, prm.core_adiabat_params)


def update(model):
    '''Update function

    Parameters
    ----------
    model : ThermalModel
        Main model class
    '''

    sl = model.stable_layer
    prm = model.parameters

    #Update profiles
    sl.profiles = sl._next_profiles.copy()
    sl.profiles['Ta'] = prof.adiabat(sl.profiles['r'], model.core.Tcen, prm.core_adiabat_params)


    #Update core profiles if no core model is being used.
    if not prm.core:
        prof.basic_profiles(model)

    #Update layer size
    sl.layer_thickness -= sl.ds_dt * model.dt
    if sl.layer_thickness < 0:
        sl.layer_thickness = 0
        

    model.core.rs = prm.r_cmb - sl.layer_thickness
    model.core.T_cmb = sl.profiles['T'][-1]

    #If layer covers entire core, Tcen is temp at base of layer.
    if sl.layer_thickness == prm.r_cmb:
        model.core.Tcen = sl.profiles['T'][0]


#Required parameters. 'Name': 'Description'
required_params = {'core_liquid_density_params': 'Outer core density polynomials (radial). List(float)',
                   'r_cmb': 'CMB radius. Float',
                   'core_alpha_T_params': 'Core thermal expansivity pressure polynomial coefficients. List(Float)',
                   'core_cp_params': 'Core specific heat capacity pressure polynomial coefficients. List(Float)',
                   'core_conductivity_params': 'List, first element is a string followed by the core thermal conductivity polynomial coefficients. The string (either \'r\'/\'T\'/\'P\') indicates if the polynomials are radial, temperature, or pressire coefficients. List(string, Float, Float...)',
                   'P_cmb': 'CMB pressure. Float',
                   'diffusivity_c': 'Chemical diffusivities of alloying light elements.  List(float)'}

#Optional parameters. 'Name': ('Description', default_value)
optional_params = {'entrainment_c': ('(Default: 0) Float in the range 0 <= x < 1. Entrainment coefficient, modifies the mass flux out the convecting region by the factor (1-x) to represent entrainment of stable fluid at the base of the stratified layer.', 0),
                   'max_resolution': ('(Default: 100) Maximum number of nodes in grid. Int', 100),
                   'min_resolution': ('(Default: 10) Minumum number of nodes in grid. Int', 10),
                   'resolution': ('(Default: 1/1000) Target number of nodes per meter when constructing the grid. Float.', 1/1000),
                   'depth_tolerance': ('(Default: 1000) Layer size below which the layer is assumed full eroded. Note this default value will be halved when considering chemically stable layers fot stability. Float', 1000),
                   'init_size': ('(Default: 5000) Size the stable layer should be initialised if conditions promote layer growth and a layer does not yet exist. Note this default value will be halved when considering chemically stable layers fot stability. Float', 5000),
                   'mix_layer': ('(Default: False) Boolean, if True, uses an experimental function to calculate the T/c profile of a layer that has the top portion mixed.', False),
                   'ignore_subadiabatic': ('(Default: False) Boolean, if True will ignore the requirement of a super-adiabatic CMB heat flow for this method. A minimum mass flux is imposed at rs instead.', False),
                   'layer_thickness': ('(Default: 0). Float, Initial stable layer thickness.',0),}


def progress(model):
    '''text to display during calculation

    Parameters
    ----------
    model : ThermalModel
        Main model class

    Returns
    -------
    str
        text to display to STDOUT
    '''

    sl = model.stable_layer

    v = (sl.layer_thickness/1000, sl.Q_cmb/1e12, sl.ADR)

    text = f'    layer thickness: {v[0]:.2f} km    Q_cmb: {v[1]:.2f} TW    ADR(rc): {v[2]:.2f}'

    return text


#3 Methods for 3 different cases
def pure_chemical_method(model):

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

    TOL = prm.depth_tolerance   #Minimum layer size
    init_size = prm.init_size

    core_adiabat      = prm.core_adiabat_params
    rho_poly          = prm.core_liquid_density_params

    E_c = prm.entrainment_c

    Q_cmb  = model.mantle.Q_cmb

    Tcen   = model.core.Tcen
    dTa_dt = model.core.dT_dt

    #Assumes the light element of interest is always the element
    #corresponding to the first element in lists/arrays for light elements.
    alpha_c = prm.alpha_c[0]
    alpha_D = sl.alpha_D
    D       = prm.diffusivity_c[0]


    dTa_dt = core.dT_dt
    dc_dt  = core.dc_dt[0].copy()

    #Make sure CMB heat flow is super-adiabatic. 
    T_grad_cmb = Q_cmb/(-model.core.profiles['k'][-1]*4*np.pi*r_cmb**2)

    ADR = T_grad_cmb/model.core.profiles['dTa_dr'][-1]
    sl.ADR = ADR

    assert sl.ADR >= 1 or prm.ignore_subadiabatic, 'ADR <= 1! This method requires a super-adiabatic heat flow'

    
    #Chemical profile
    r = sl.profiles['r']
    c = sl.profiles['c']

    TOL = prm.depth_tolerance   #Minimum layer size
    init_size = prm.init_size

    #For stability, reduce the time step size when layer is relatively thin.
    factor = 1
    if layer_thickness <= init_size:
        factor = 0.1

    time_gone = 0
    while time_gone < model.dt:

        #Initialise layer with minium tolerance thickness
        if layer_thickness < TOL:

            layer_thickness = init_size
            r_s = r_cmb - layer_thickness

            #Number of grid points in this iteration
            n_points = int(np.max([min_res,layer_thickness*resolution]))
            n_points = int(np.min([max_res,n_points]))

            r = np.linspace(r_s, r_cmb, n_points)
            c = np.ones(r.size)*core.conc_l[0]


        #####################
        dt_small = model.dt*factor
        #####################
        if time_gone + dt_small > model.dt:
            dt_small = model.dt - time_gone

        #Enrich convecting region and cool adiabat
        conc_l = core.conc_l[0] + dc_dt *time_gone
        Tcen   = core.Tcen      + dTa_dt*time_gone

        #Regrid core radial profiles onto finer stable layer grid
        rho     = prof.density(r, rho_poly)
        T       = prof.adiabat(r, Tcen, core_adiabat)
        alpha_T = np.interp(r, core.profiles['r'], core.profiles['alpha'])
        g       = np.interp(r, core.profiles['r'], core.profiles['g'])
        


        #Lower boundary condition (using entrainment)
        lb = -(alpha_T[0]/alpha_c)*(T_grad_cmb-prof.adiabat_grad(r[0], Tcen, core_adiabat))*(1-E_c)
        lb_type = 1

        #For stability, very small chemical gradients cause v.large jumps in rs.
        if lb < 1e-10:
            lb = 1e-9
            
        #If whole core is stratified, lower BC should be zero flux
        if r[0] == 0:
            lb_type = 1
            lb = 0
        
        #Upper boundary condition
        ub_type, ub = model.mantle.chemical_bc_cmb
        #Either this tuple is calculated in a mantle model (e.g. driscol_bercovici14 bmo evolution).
        #or is set by the user (equivalent to needing to set Q_cmb if no mantle model is used).
        #This specifies also the type of BC (neumann or dirchlet).

        if prm.include_baro_diffusion:

            #Calculate barodiffusive flux in spherical coordinates.
            dg_dr = np.gradient(g,r[1]-r[0],edge_order=2)
            a = alpha_D*alpha_c/rho  
            baro = -(2*g*a/r +a*dg_dr)

        else:
            baro = 0


        #Calculate diffusion solution
        c_new = diffusion(c, r, dt_small, D, 0, 0, (lb_type,lb), (ub_type,ub), constant=baro)
        T_new = np.zeros(r.size)


        #Mix layer
        if prm.mix_layer:
            c_rel = func.mix_profile(r, rho*r**2, c_new-conc_l)
            c_new = c_rel + conc_l

        sl.c_grad_s = (c_new[1]-c_new[0])/(r[1]-r[0])


        try:
            r_s_new = func.retreat_layer(r, T_new, c_new, 0, conc_l, core_adiabat,
                                        (resolution,min_res,max_res), alpha_T, alpha_c, density_grad_limit=-lb*alpha_c)
        except:
            breakpoint()


        if r_s_new == r[0]:

            r_s_new = func.buffett_seagle_10_growth(r, c_new, conc_l, lb)

            if r_cmb - r_s_new > init_size:
                factor = 1

        else:
            #Give up if layer will still not grow after generous time step increases.
            if r_s_new >= r_cmb-TOL:

                #No stratification
                r_s_new = prm.r_cmb
                factor = factor*2

                if factor > 64:
                    dt_small = model.dt - time_gone
                    factor = 1
            

        layer_thickness = r_cmb - r_s_new

        c_rel = c_new - conc_l

        r, c_rel = func.change_domain_size(c_rel, r, r_s_new, (resolution,min_res,max_res))

        c = c_rel + conc_l

        time_gone += dt_small


    #Save new profiles. Keep original profiles until update() has been called.
    sl._next_profiles['r'] = r
    sl._next_profiles['c'] = c
    sl._next_profiles['T'] = prof.adiabat(r, core.Tcen + core.dT_dt*model.dt, core_adiabat)

    sl.ds_dt = (r_s_new - r_s_original)/model.dt

    if r[1]-r[0] > 0:
        Q = model.mantle.Q_cmb - en.secular_cool(sl.profiles['r'], sl.profiles['rho'], sl.profiles['T'], sl.profiles['cp'], sl.profiles['r'].size)[1]*core.dT_dt
        sl.T_grad_s = Q/(-4*np.pi*core.rs**2 * sl.profiles['k'][0])
        sl.c_grad_s = (c[1]-c[0])/(r[1]-r[0])
    else:
        sl.T_grad_s = 0
        sl.c_grad_s = 0

    sl.lb_T = lb
    sl.ub_T = ub


