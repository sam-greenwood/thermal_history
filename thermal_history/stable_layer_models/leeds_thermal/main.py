Description = 'General stable layer model for thermal/chemical/thermo-chemical stratification.'


#List individually confirmed compatibility with other regions
compatibility = {'core': ['greenwood21'],
                 'mantle': ['driscoll_bercovici14']}

from ...core_models.leeds.routines import profiles as prof
from ...core_models.leeds.routines import energy as en
from ...core_models.leeds.routines import chemistry as chem

from .routines import functions as func
from .routines.diffusion import diffusion

from thermal_history.utils.optimised_funcs import linspace, trapezoid, polyval

import numpy as np
from scipy.optimize import bisect
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


    #If no core model is being used, assume simple case of fully molten core
    if not prm.core:
        core.T_cmb = prm.T_cmb
        core.ri = 0
        core.r_snow = prm.r_cmb

        #Set core profiles
        prof.basic_profiles(model)
        prof.temp_dependent_profiles(model)


    #Initialise profiles dictionaries
    sl.profiles = {}
    sl._next_profiles = {}

    n = prm.max_resolution

    #Create radial profiles
    sl.profiles['r'] = linspace(prm.r_cmb - prm.layer_thickness, prm.r_cmb, n)

    sl.profiles['T'] = prof.adiabat(sl.profiles['r'], core.Tcen, prm.core_adiabat_params)
    sl.T_grad_s      = prof.adiabat_grad(core.rs, core.Tcen, prm.core_adiabat_params)



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
    check  = func.pure_thermal_check(model)

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
            Ek = 4*np.pi*trapz(k*(dT_dr/T_initial)**2*r_initial**2, x=r_initial)


        # rho =  prof.density(r_initial, prm.core_density_liquid_params)
        # g   =  np.interp(r_initial, core.profiles['r'], core.profiles['g'])

    
        #Time step model
        pure_thermal_method(model)

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


        cp  = np.interp(r1, core.profiles['r'], core.profiles['cp'])
        rho = prof.density(r1, prm.core_liquid_density_params)
        Qs = -4*np.pi*en.integrate(r1, r1**2*dT_dt*rho*cp)
        Es = -4*np.pi*en.integrate(r1, r1**2*dT_dt*rho*cp*(1/T1[-1] - 1/T1))
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
                   'P_cmb': 'CMB pressure. Float'}

#Optional parameters. 'Name': ('Description', default_value)
optional_params = {'entrainment_T': ('(Default: 0) Float in the range 0 <= x < 1. Entrainment coefficient, modifies the adiabatic heat flow out the convecting region by the factor (1-x) to represent entrainment of stable fluid at the base of the stratified layer.', 0),
                   'max_resolution': ('(Default: 100) Maximum number of nodes in grid. Int', 100),
                   'min_resolution': ('(Default: 10) Minumum number of nodes in grid. Int', 10),
                   'resolution': ('(Default: 1/1000) Target number of nodes per meter when constructing the grid. Float.', 1/1000),
                   'depth_tolerance': ('(Default: 1000) Layer size below which the layer is assumed full eroded. Note this default value will be halved when considering chemically stable layers fot stability. Float', 1000),
                   'init_size': ('(Default: 5000) Size the stable layer should be initialised if conditions promote layer growth and a layer does not yet exist. Note this default value will be halved when considering chemically stable layers fot stability. Float', 5000),
                   'mix_layer': ('(Default: False) Boolean, if True, uses an experimental function to calculate the T/c profile of a layer that has the top portion mixed.', False),
                   'thermal_stratification': ('(Default: True). Boolean, internal flag to tell the code that thermal stratification is being used', True),
                   'chemical_stratification': ('(Default: False). Boolean, internal flag to tell the code that chemical stratification is not being used', False),
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

def pure_thermal_method(model):
    '''Pure thermal stratification

    Solves the thermal diffusion solution and updates the layer size during one time iteration.
    For stability, it may step in smaller time increments until the total model timestep has been
    reached.

    Parameters
    ----------
    model : ThermalModel
        Main model class

    
    '''

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

    E_T = prm.entrainment_T

    Q_cmb  = model.mantle.Q_cmb

    Tcen   = model.core.Tcen
    dTa_dt = model.core.dT_dt

    T_grad_cmb = Q_cmb/(-model.core.profiles['k'][-1]*4*np.pi*r_cmb**2)

    ADR = T_grad_cmb/model.core.profiles['dTa_dr'][-1]
    sl.ADR = ADR


    r = sl.profiles['r']
    T = sl.profiles['T']


    #Set upper boundary conditions
    ub = T_grad_cmb

    #For stability, reduce the time step size when layer is relatively thin.
    factor = 1
    if layer_thickness <= init_size*3:
        factor = 0.01

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
            T = prof.adiabat(r, Tcen, core_adiabat)


        # dt_small = factor*(prm.ys*1e4 * layer_thickness/300) #Time step scaling for stability
        # if dt_small > model.dt:
        #     dt_small = model.dt

        dt_small = model.dt*factor

        if time_gone + dt_small > model.dt:
            dt_small = model.dt - time_gone

        #Regrid core radial profiles onto finer stable layer grid
        rs_idx = model.core._rs_idx
        rho = prof.density(r, rho_poly)
        k   = np.interp(r, core.profiles['r'], core.profiles['k'])
        cp = np.interp(r, core.profiles['r'], core.profiles['cp'])
        alpha = np.interp(r, core.profiles['r'], core.profiles['alpha'])

        #Thermal diffusivity
        D_T = k/(rho*cp)
        #dk_dr * 1/(rho*cp). Not the same as the diffusivity gradient
        dD_dr = np.gradient(k, r, edge_order=2)/(rho*cp)


        #Cool adiabat
        Tcen = model.core.Tcen + dTa_dt*(time_gone+dt_small)

        Ta, Ta_grad = prof.adiabat(r, Tcen, core_adiabat), prof.adiabat_grad(r, Tcen, core_adiabat)

        #Lower boundary condition (using entrainment)
        lb = (1-E_T)*Ta_grad[0]
        lb_type = 1

        #If whole core is stratified, lower BC should be zero flux
        if r[0] == 0:
            lb = 0


        #Calculate diffusion solution
        T_new = diffusion(T, r, dt_small, D_T, k, dD_dr, (lb_type,lb),(1,ub))

        if r[0] == 0:
            Tcen = T_new[0]

        if np.isnan(T_new[0]):
            breakpoint()


        #Mix layer (experimental function, not used by default)
        if prm.mix_layer:
            T_rel = func.mix_profile(r, cp*rho*r**2, T_new-Ta)
            T_new = T_rel + Ta

        #Determine if layer should retreat
        r_s_new = func.retreat_layer(r, T_new, np.zeros(r.size), Tcen, 0, core_adiabat, (resolution,min_res,max_res), alpha, 0)
        if r[0] == 0 and ADR < 1:
            r_s_new = 0


        #Grow layer
        if r_s_new == r[0]:

            #Find radius on adiabat that matches temperature at base of layer
            def f(guess,T,Tcen):
                return T - prof.adiabat(guess,Tcen,core_adiabat)
            
            if T_new[0] >= Tcen: #Layer has reached center of core
                r_s_new = 0
            else:
                r_s_new = bisect(f,0,r_cmb,args=(T_new[0],Tcen),maxiter=200)

        #Shrink layer
        else:
            #Give up if layer will still not grow after generous time step increases.
            if r_s_new >= r_cmb-TOL:

                #No stratification
                r_s_new = prm.r_cmb
                factor = factor*2

                if factor > 64 or ADR >= 1:
                    dt_small = model.dt - time_gone
                    factor = 1

        layer_thickness = r_cmb - r_s_new

        #Regrid domain, keeping temperature relative to the adiabat.
        T_rel = T_new - prof.adiabat(r, Tcen, core_adiabat)

        r, T_rel = func.change_domain_size(T_rel, r, r_s_new, (resolution,min_res,max_res))

        T = T_rel + prof.adiabat(r, Tcen, core_adiabat)

        time_gone += dt_small


    #Save new profiles. Keep original profiles until update() has been called.
    sl._next_profiles['r'] = r
    sl._next_profiles['T'] = T

    sl.ds_dt = (r_s_new - r_s_original)/model.dt

    if r[1]-r[0] > 0:
        sl.T_grad_s = (T[1]-T[0])/(r[1]-r[0])
    else:
        sl.T_grad_s = 0

    sl.lb_T = lb
    sl.ub_T = ub
