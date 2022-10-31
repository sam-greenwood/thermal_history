Description = 'Thermal stratification with embedded FeS layer for Mercury'


#List individually confirmed compatibility with other regions
compatibility = {'core': ['leeds'],
                 'mantle': ['knibbe_westrenen18']}

from builtins import breakpoint
from ...core_models.leeds.routines import profiles as prof
from ...core_models.leeds.routines import energy as en
from ...core_models.leeds.routines import chemistry as chem

from .routines import functions as func
from .routines.diffusion import diffusion2, diffusion
from .routines.fes_functions import tanh_conductivity, adaptive_grid

from thermal_history.utils.optimised_funcs import linspace, trapezoid, polyval

import time
import numpy as np
from scipy.optimize import bisect
from scipy.integrate import trapz
from copy import deepcopy

import matplotlib.pyplot as plt

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

    sl.layer_thickness = prm.FeS_size

    sl.layer_thickness = prm.layer_thickness

    if not prm.core:
        raise ValueError('Needs a core model')

    #Initialise profiles dictionaries
    sl.profiles = {}
    sl._next_profiles = {}

    n = prm.max_resolution

    #Create radial profiles
    sl.profiles['r'] = linspace(core.rs, prm.r_cmb, n)
    sl.profiles['T'] = np.interp(sl.profiles['r'], core.profiles['r'], core.profiles['T'])

    # sl.profiles['T'] = prof.adiabat(sl.profiles['r'], core._Tcen_fes, prm.fes_adiabat_params) #FeS adiabat
    sl.T_grad_s      = prof.adiabat_grad(core.rs, core.Tcen, prm.core_adiabat_params) #Bulk core adiabatic gradient at rs



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

    r_fes = prm.r_cmb - prm.FeS_size #radius of FeS layer
    fes_idx = core._fes_idx          #Index for FeS layer in core profiles



    #Initialise energies/entropies associated with this state.
    Qs, Es, Ek, Ej = 0, 0, 0, 0

    #initialise BC's
    sl.lb_T, sl.ub_T = 0,0
    sl.T_grad_s = 0

    sl.alpha_D = 0
    sl.baro_grad = 0


    if sl.layer_thickness > prm.FeS_size:
        T_grad_fes = np.interp(r_fes, sl.profiles['r'], np.gradient(sl.profiles['T'], sl.profiles['r'], edge_order=2))
    else:
        T_grad_fes = core.profiles['dTa_dr'][core._fes_idx-1]

    heat_in = 4*np.pi*r_fes**2 * core.profiles['k'][core._fes_idx-1] * -T_grad_fes

    #Check if FeS layer is sub-adiabatic
    ADR_fes = (core.Q_cmb-heat_in)/core.Qa
    ADR_fes = core.Q_cmb/core.Qa

    FeS_convecting = True
    if ADR_fes < 1:
        FeS_convecting = False

    FeS_convecting = True #Force convecting FeS layer

    #If no thermal stratification exists yet, re-set secular cooling of bulk.
    if sl.layer_thickness == prm.FeS_size:

        #Cooling rates in core calculation assumes adiabatic heat flow at rs, correct it here to account for
        #secular cooling of FeS layer as well.

        r, rho, T, cp = core.profiles['r'], core.profiles['rho'], core.profiles['T'], core.profiles['cp']
        Qs_tilde, Es_tilde = en.secular_cool(r, rho, T, cp, r.size, Tcmb=T[-1]) #Secular cooling of whole core (inc FeS layer)
        Ql_tilde = core.Ql/core.dT_dt   #Re-normalise values
        Qg_tilde = core.Qg/core.dT_dt

        dT_dt = (core.Q_cmb-core.Qr)/(Qs_tilde+Ql_tilde+Qg_tilde)  #New corrected cooling rate of core

        dT_dt_fes = dT_dt * core._Tcen_fes/core.Tcen  #FeS cooling rate

        core.Q_fes = core.Q_cmb - dT_dt*(Qs_tilde - core.Qs/core.dT_dt) #Heat flow at base of FeS layer (Q_cmb minus secular cooling of FeS layer)

        #Sub-adiabatic. reset Q_fes to Qa at r_fes and layer will start to grow.
        if core.Q_fes < core.Qa_rs:

            # Q_fes = core.Qa_rs #Set to adiabatic heat flow.
            dT_dt = core.dT_dt

            # Qs_fes = core.Q_cmb-Q_fes

            # dT_dt_fes = Qs_fes/Qs_tilde


            # dT_dt = (core.Q_cmb-core.Qr)/(Qs_tilde+Ql_tilde+Qg_tilde)  #New corrected cooling rate of core

        #New ADR_fes
        ADR_fes = core.Q_cmb/core.Qa

        profiles = core.profiles
        r, T = profiles['r'], profiles['T']
        idx = core._fes_idx
        Qs_tilde, Es_tilde = en.secular_cool(profiles['r'][idx:], profiles['rho'][idx:], profiles['T'][idx:], profiles['cp'][idx:], prm.n_profiles-idx) #Just FeS layer integrals

        Qs_tilde = Qs_tilde * profiles['T'][idx]/core._Tcen_fes #Renormalise from Tcen value to theoretical FeS adiabat at Tcen
        Es_tilde = Es_tilde * profiles['T'][idx]/core._Tcen_fes

        sl.dT_dr_fes = profiles['dTa_dr'][idx]

        
            
        # if model.it == 1 and core.Q_fes < core.Qa_rs:

        #     sl.profiles['r'] = linspace(0, r_fes, 60)
        #     sl.profiles['T'] = prof.adiabat(sl.profiles['r'], core.Tcen, prm.core_adiabat_params)
        #     sl.layer_thickness = prm.r_cmb
        #     dT_dt = 0


        Qs += Qs_tilde*dT_dt_fes    #Add secular cooling of just FeS layer to secular cooling of total layer.
        Es += Es_tilde*dT_dt_fes


        #Correct core energies/entropies
        core.Qs = core.Qs*dT_dt/core.dT_dt
        core.Ql = core.Ql*dT_dt/core.dT_dt
        core.Qg = core.Qg*dT_dt/core.dT_dt

        core.dri_dt = core.Cr*core.dT_dt
        core.dc_dt  = core.Cr*core.Cc*core.dT_dt

        core.Es = core.Es*dT_dt/core.dT_dt
        core.El = core.El*dT_dt/core.dT_dt
        core.Eg = core.Eg*dT_dt/core.dT_dt

        core.Ej = core.Es + core.Eg + core.El - core.Ek - core.Ea

        #Correct cooling rates
        core.dT_dt = dT_dt
        core._dT_dt_fes = dT_dt_fes

    elif FeS_convecting:

        #FeS layer is convecting
        
        profiles = core.profiles
        r, T = profiles['r'], profiles['T']
        idx = core._fes_idx

        #FeS layer needs mixing as convection starts again. sl.ADR_fes is value from previous timestep.
        if model.it > 1 and (sl.ADR_fes < 1 and ADR_fes >= 1):

            total_heat = trapezoid(r[idx:], r[idx:]**2 * T[idx:])[-1] #Conserve total thermal mass

            core._Tcen_fes = total_heat / trapezoid(r[idx:], r[idx:]**2 * prof.adiabat(r[idx:], 1, prm.fes_adiabat_params))[-1]

        #Normliased energies for convecting FeS layer
        Qs_tilde, Es_tilde = en.secular_cool(profiles['r'][idx:], profiles['rho'][idx:], profiles['T'][idx:], profiles['cp'][idx:], prm.n_profiles-idx) #Just FeS layer integrals
        Qs_tilde *= profiles['T'][idx]/core._Tcen_fes #Renormalise from Tcen value to theoretical FeS adiabat at Tcen

        # breakpoint()
        # core.Q_fes = (core.Qs + core.Qg + core.Ql) + sl._Qs_strat

        # core.Q_fes = core.Q_cmb #TEST

        # #Heat flow into base of FeS layer
        # T_grad_fes = np.interp(r_fes, sl.profiles['r'], np.gradient(sl.profiles['T'], sl.profiles['r'], edge_order=2))
        T_grad_fes = (T[idx-1]-T[idx-2])/(r[idx-1]-r[idx-2])
        core.Q_fes = 4*np.pi*r_fes**2 * profiles['k'][idx-1] * -T_grad_fes

        core._dT_dt_fes = (core.Q_cmb - core.Q_fes)/Qs_tilde

        sl.Qs_fes = Qs_tilde*core._dT_dt_fes

        # if model.it >= 7574:
        #     print(core._Tcen_fes, core._dT_dt_fes, prof.adiabat(prm.r_cmb, core._Tcen_fes, prm.fes_adiabat_params))
        #     breakpoint()

    else:
        #Secular cooling of FeS layer will be accounted for in conduction solution.

        profiles = core.profiles
        r, T = profiles['r'], profiles['T']
        idx = core._fes_idx

        #Set theoretical Tcen to track the temperature at r_fes
        core._Tcen_fes = T[idx] / prof.adiabat(r[idx], 1, prm.fes_adiabat_params)
        


    sl.ADR_fes = ADR_fes
    # sl.ADR_fes = 1.1 #Force adiabatic FeS layer for now


    #Check if bulk core is adiabatic
    ADR = core.Q_fes/core.Qa_fes
    sl.ADR = ADR

    #Estimate based on Qc rather than Q_fes
    sl.ADR = core.Q_cmb / (-4*np.pi*prm.r_cmb**2 * profiles['k'][core._fes_idx-1] * prof.adiabat_grad(prm.r_cmb, core.Tcen, prm.core_adiabat_params))

    #Profiles before diffusion solution.
    r_initial = sl.profiles['r']
    T_initial = sl.profiles['T']

    #Calculate Ek
    dT_dr = np.gradient(T_initial, r_initial, edge_order=2)
    k = np.interp(r_initial, core.profiles['r'], core.profiles['k'])
    Ek = 4*np.pi*trapz(k*(dT_dr/T_initial)**2*r_initial**2, x=r_initial)

    #If no layer exists and bulk is super-adiabatic, don't bother with diffusion solution.
    if sl.ADR >= 1 and sl.layer_thickness == prm.FeS_size:
        sl._next_profiles['r'] = np.linspace(r_fes, prm.r_cmb, 10)
        sl._next_profiles['T'] = prof.adiabat(sl._next_profiles['r'], core._Tcen_fes + model.dt*core._dT_dt_fes, prm.fes_adiabat_params) # np.full(10, core.Tcen + model.dt*core.dT_dt)
        sl.ds_dt = 0
    else:
        #Time step layer evolution
        pure_thermal_method(model)

    #New radial profiles
    r_new = sl._next_profiles['r']
    T_new = sl._next_profiles['T']


    #Secular cooling. Compare final profiles to those at the beginning of the time-step.
    if r_new[0] <= r_initial[0]:  #Layer has grown

        #Need to append adiabat values onto initial temperatures
        T_initial_temp = prof.adiabat(r_new, core.Tcen, prm.core_adiabat_params)

        for i in range(r_new.size):
            if r_new[i] >= r_initial[0]:
                T_initial_temp[i] = np.interp(r_new[i], r_initial, T_initial)
        T_initial = T_initial_temp

        r1, T2, T1 = r_new, T_new, T_initial

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
    rho = np.interp(r1, core.profiles['r'], core.profiles['rho'])
    Qs += -4*np.pi*en.integrate(r1, r1**2*dT_dt*rho*cp)
    Es += -4*np.pi*en.integrate(r1, r1**2*dT_dt*rho*cp*(1/T1[-1] - 1/T1))
    sl.dT_dt_s = dT_dt[0]  #Save cooling rate at base of layer

    # if core.ri > 0:
    #     breakpoint()


    # breakpoint()

    Ej = Es - Ek

    #Save energies/entropies to model attributes
    sl.Qs, sl.Es, sl.Ek, sl.Ej = Qs, Es, Ek, Ej

    # if sl.layer_thickness > prm.FeS_size:
    #     breakpoint()

    #Save secular cooling of stratified fluid for next iteration for calculating Q_fes.
    if r_new[0] < r_fes:
        sl._Qs_strat = -4*np.pi*en.integrate(r1[r1<=r_fes], r1[r1<=r_fes]**2 * dT_dt[r1<=r_fes] * rho[r1<=r_fes] * cp[r1<=r_fes]) #Secular cooling of stratified fluid not in FeS layer

    #Add onto core results
    if prm.core:
        core.Qs += Qs
        core.Ej += Ej
        core.Es += Es
        core.Ek += Ek

        # if r_new[0] <= r_fes and ADR_fes >=1: #Convecting FeS layer secular cooling

        #     profiles = core.profiles
        #     r, T = profiles['r'], profiles['T']

        #     idx = core._fes_idx

        #     #Temp gradient at r_fes
        #     dT_dr_fes = (T[idx-1]-T[idx-2])/(r[idx-1]-r[idx-2])

        #     #Conductive heat flow from bulk into FeS layer.
        #     core.Q_fes = -profiles['k'][idx-1]*dT_dr_fes*4*np.pi*r_fes**2

        #     Qs_fes = core.Q_cmb-core.Q_fes

        #     Qs_tilde, Es_tilde = en.secular_cool(profiles['r'][idx:], profiles['rho'][idx:], profiles['T'][idx:], profiles['cp'][idx:], prm.n_profiles-idx)
        #     Qs_tilde = Qs_tilde * profiles['T'][idx]/core._Tcen_fes
        #     Es_tilde = Es_tilde * profiles['T'][idx]/core._Tcen_fes #Renormalise from first T value to theoretical FeS adiabat at Tcen
        #     dT_dt_fes = Qs_fes/Qs_tilde

        #     #Add on energy and entropy for FeS layer
        #     core.Qs += Qs_fes
        #     core.Es += Es_tilde*dT_dt_fes
        #     core._dT_dt_fes = dT_dt_fes

    #Make sure profiles have size maximum_resolution, so they can be saved.
    n = sl.profiles['r'].size

    if n > prm.max_resolution:

        r = linspace(sl.profiles['r'][0], sl.profiles['r'][-1], prm.max_resolution)
        T = np.interp(r, sl.profiles['r'], sl.profiles['T'])

    elif n < prm.max_resolution:

        r = np.append(sl.profiles['r'], np.ones(prm.max_resolution-n)*sl.profiles['r'][-1])
        T = np.append(sl.profiles['T'], np.ones(prm.max_resolution-n)*sl.profiles['T'][-1])

    else:
        r, T = sl.profiles['r'], sl.profiles['T']


    sl.profiles['r'], sl.profiles['T'] = r, T

    sl.profiles['Ta'] = prof.adiabat(r, core.Tcen, prm.core_adiabat_params)

    sl.T_cmb, sl.T_s= T[-1], T[0]

    core.T_cmb,  core.T_s = sl.T_cmb, sl.T_s

    if prm.FeS_size == 0:
        core._Tcen_fes = core.Tcen
        core._dT_dt_fes = core._dT_dt_fes

    # sl.dT_dt_fes = (np.interp(r_fes, r, T)-np.interp(r_fes-1, r ,T))


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
        model.core.dT_dt = sl.dT_dt_s

    #Track theoretical Tcen for fes layer if it is part of the conduction solution.
    if sl.ADR_fes < 1:
        model.core._Tcen_fes = model.core.T_cmb / prof.adiabat(prm.r_cmb, 1, prm.fes_adiabat_params)


#Required parameters. 'Name': 'Description'
required_params = {'core_liquid_density_params': 'Outer core density polynomials (radial). List(float)',
                   'r_cmb': 'CMB radius. Float',
                   'core_alpha_T_params': 'Core thermal expansivity pressure polynomial coefficients. List(Float)',
                   'core_cp_params': 'Core specific heat capacity pressure polynomial coefficients. List(Float)',
                   'core_conductivity_params': 'List, first element is a string followed by the core thermal conductivity polynomial coefficients. The string (either \'r\'/\'T\'/\'P\') indicates if the polynomials are radial, temperature, or pressire coefficients. List(string, Float, Float...)',
                   'P_cmb': 'CMB pressure. Float',
                   'FeS_size': 'Size of the FeS layer'}

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

    v = (sl.layer_thickness/1000, sl.Q_cmb/1e12, sl.ADR, sl.T_cmb)

    text = f'    layer thickness: {v[0]:.2f} km    Q_cmb: {v[1]:.2f} TW    ADR(rc): {v[2]:.2f}    T_cmb: {v[3]:.2f} K'

    return text


#Thermal method with FeS

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
    init_size = prm.init_size + prm.FeS_size

    core_adiabat      = prm.core_adiabat_params
    rho_poly          = prm.core_liquid_density_params

    E_T = prm.entrainment_T

    Q_cmb  = model.mantle.Q_cmb
    Q_fes = core.Q_fes

    Tcen   = model.core.Tcen
    dTa_dt = model.core.dT_dt

    r_fes = prm.r_cmb - prm.FeS_size

    r_upper = prm.r_cmb

    r = sl.profiles['r']
    T = sl.profiles['T']

    #Set upper boundary conditions
    #If FeS layer is sub-adiabatic, solution as normal
    ub = prof.adiabat_grad(prm.r_cmb, core._Tcen_fes, prm.fes_adiabat_params)*sl.ADR_fes
    ub = core.Q_cmb/(4*np.pi*prm.r_cmb**2*-core.profiles['k'][-1])
    ub_type = 1

    # if sl.ADR_fes >= 1:
    #     sl.ADR_fes = 0.99

    if sl.ADR_fes >= 1:
        r_upper = r_fes

    r_upper = r_fes #Force convecting FeS layer for now.


    #For stability, reduce the time step size when layer is relatively thin.
    factor = 1
    # if layer_thickness <= init_size*3:
    #     factor = 0.01

    time_gone = 0
    while time_gone < model.dt:

        coarse = prm.coarse_resolution
        fine   = prm.fine_resolution

        #Initialise layer with minium tolerance thickness
        if layer_thickness < TOL + prm.FeS_size:

            layer_thickness = init_size

            #Number of grid points in this iteration
            n_points = int(np.max([min_res,layer_thickness*resolution]))
            n_points = int(np.min([max_res,n_points]))

            r_s = prm.r_cmb - layer_thickness

            # r_test = np.linspace(r_s, r_upper, n_points)
            

            r = adaptive_grid(r_s, r_fes, r_upper, coarse, fine, transition_width=200e3)

            T = prof.adiabat(r, Tcen, core_adiabat)
            # T_test = prof.adiabat(r_test, Tcen, core_adiabat)

        elif r_upper == r_fes:
            #FeS layer is still convecting so just solve for bulk core that is thermally stratified.

            #Number of grid points in this iteration
            n_points = int(np.max([min_res,int((r_upper-r_s)/coarse)]))
            n_points = int(np.min([max_res,n_points]))

            r_s = prm.r_cmb - layer_thickness

            # coarse = np.min([10e3, (r_upper-r_s)/10])
            # fine = np.min([2e3, (r_upper-r_s)/20])
            r = np.linspace(r_s, r_upper, n_points)
            # r = adaptive_grid(r_s, r_fes, r_upper, coarse, fine, transition_width=200e3)

            T = np.interp(r, sl.profiles['r'], sl.profiles['T'])

        else:
            #Make sure solution is on even grid
            n_points = int(np.max([min_res,layer_thickness*resolution]))
            n_points = int(np.min([max_res,n_points]))

            # coarse = np.min([10e3, (r_upper-r_s)/10])
            # fine = np.min([2e3, (r_upper-r_s)/20])

            r = adaptive_grid(r_s, r_fes, r_upper, coarse, fine, transition_width=200e3)
            # r = np.linspace(sl.profiles['r'][0], sl.profiles['r'][-1], n_points)
            T = np.interp(r, sl.profiles['r'], sl.profiles['T'])

            # breakpoint()


        dt_small = model.dt*factor

        if time_gone + dt_small > model.dt:
            dt_small = model.dt - time_gone

        #Regrid core radial profiles onto finer stable layer grid
        rho = prof.density(r, rho_poly)

        #n = int((r.size-1)*(layer_thickness-prm.FeS_layer_size)/layer_thickness)

        # FeS_idx = np.where(r >= prm.r_cmb-prm.FeS_layer_size)[0][0]
        # if FeS_idx >= r.size-2:
        #     breakpoint()
        # k, dk_dr = tanh_conductivity(r, core.profiles['k'][0], prm.FeS_conductivity, FeS_idx)


        if r_upper == prm.r_cmb:
            k, dk_dr = tanh_conductivity(r, core.profiles['k'][0], prm.core_fes_conductivity, r_fes, transition_width=5000)
        else:
            k = core.profiles['k']
            dk_dr = np.zeros(k.size)

            k   = np.interp(r, core.profiles['r'], k)
            dk_dr = np.interp(r, core.profiles['r'], dk_dr)

            # k_test = np.full(r_test.size, k[0])
            # dk_dr_test = np.zeros(r_test.size)



        # if layer_thickness > prm.FeS_layer_size:
        #     breakpoint()

        cp    = np.interp(r, core.profiles['r'], core.profiles['cp'])
        alpha = np.interp(r, core.profiles['r'], core.profiles['alpha'])

        #Thermal diffusivity
        D_T = k/(rho*cp)
        #dk_dr * 1/(rho*cp). Not the same as the diffusivity gradient

        # D_T_test = np.full(r_test.size, D_T[0])

        # dD_dr = np.gradient(k, r, edge_order=2)/(rho*cp)
        dD_dr = dk_dr/(rho*cp)


        #Cool adiabat
        Tcen = core.Tcen + dTa_dt*(time_gone+dt_small)
        Tcen_fes = core._Tcen_fes + core._dT_dt_fes*(time_gone+dt_small)

        Ta, Ta_grad = prof.adiabat(r, Tcen, core_adiabat), prof.adiabat_grad(r, Tcen, core_adiabat)

        Ta[r>r_fes] = prof.adiabat(r[r>r_fes], Tcen_fes, prm.fes_adiabat_params)


        #Lower boundary condition (using entrainment)
        lb = (1-E_T)*Ta_grad[0]
        lb_type = 1

        #Fes layer is convecting so change upper BC.
        if r_upper == r_fes:
            #Fixed temperature boundary condition
            ub      = prof.adiabat(r_fes, Tcen_fes, prm.fes_adiabat_params)
            ub_type = 0

        #If whole core is stratified, lower BC should be zero flux
        if r[0] == 0:
            lb = 0
   

        #Calculate diffusion solution
        T_new = diffusion2(T, r, dt_small, D_T, k, dk_dr, (lb_type,lb),(ub_type,ub))
        # T_new_test = diffusion(T_test, r_test, dt_small, D_T_test, k_test, dk_dr_test, (lb_type,lb),(ub_type,ub))

        Qs = -4*np.pi*trapezoid(r, r**2 * rho* cp *(T_new-T)/dt_small)[-1]
        # print(f'Qs: {Qs: .3e}   Q_fes: {core.Q_fes-core.Qs: .3e}')

        # print(model.it, t2-t1, np.min(np.diff(r)))
        # if model.it > 75 and ub_type == 1:

        #     Q_top = 4*np.pi*r_cmb**2 * k[-1] * -ub
        #     Q_bottom = 4*np.pi*rs**2 * k[0] * -lb

        #     Qs = 
        
        # if core.ri>0:
        #     import matplotlib.pyplot as plt
        #     plt.plot(r, T_new, label='T_new')
        #     plt.plot(r, T, label='T')
        #     plt.plot(r, Ta, label='Ta')
        #     plt.legend(loc=0)
        #     plt.savefig('test.png')
        #     breakpoint()



        #Set adiabat to track diffusion solution
        if r[0] == 0:
            Tcen = T_new[0]
            Ta, Ta_grad = prof.adiabat(r, Tcen-100, core_adiabat), prof.adiabat_grad(r, Tcen, core_adiabat)


        if np.isnan(T_new[0]):
            breakpoint()


        #Mix layer (experimental function, not used by default)
        if prm.mix_layer:
            T_rel = func.mix_profile(r, cp*rho*r**2, T_new-Ta)
            T_new = T_rel + Ta

        if sl.ADR < 1 and r[0]==0:
            r_s_new = 0
        elif r[0] == 0 and model.it == 1:
            r_s_new = 0
        else:
            #Determine if layer should retreat
            r_s_new = func.retreat_layer(r[r<=r_fes], T_new[r<=r_fes], np.zeros(r.size)[r<=r_fes], Ta[r<=r_fes], 0, alpha[r<=r_fes], 0, density_grad_limit=1e-6)

        # print(f'it: {model.it}   Qs: {core.Qs/1e12: .3f}   Qa: {core.Qa_rs/1e12: .3f}   size: {sl.layer_thickness/1000}    dT_dt_fes:{core._dT_dt_fes: .3e}')

        #If whole core is stratified, only ADR>=1 can change the layer thickness.
        if r[0] == 0 and sl.ADR < 1:
            r_s_new = 0

        #Grow layer
        elif r_s_new == r[0]:

            #Find radius on adiabat that matches temperature at base of layer
            def f(guess,T,Tcen):
                return T - prof.adiabat(guess,Tcen,core_adiabat)

            if T_new[0] >= Tcen: #Layer has reached center of core
                r_s_new = 0
            else:
                r_s_new = bisect(f,0,r_upper,args=(T_new[0],Tcen),maxiter=200)

        #Shrink layer
        else:
            # import matplotlib.pyplot as plt
            # plt.plot(r, T_new, label='T_new')
            # plt.plot(r, T, label='T')
            # plt.plot(r, Ta, label='Ta')
            # plt.legend(loc=0)
            # plt.show()
            # breakpoint()

            # # plt.savefig('test.png')
            # breakpoint()

            if r_s_new >= r_upper-TOL:

                #No stratification. Try a larger timestep.
                r_s_new = r_upper
                factor = factor*2

                if factor > 64 or sl.ADR >= 1: #Give up is large time step doesn't help or core is super-adiabatic.
                    dt_small = model.dt - time_gone
                    factor = 1

        #r_s must be include FeS layer.
        # r_s_new = np.min([r_s_new, prm.r_cmb-prm.FeS_layer_size])
        
        #If FeS layer is sub-adiabatic, r_s is at least at base of FeS layer.
        if sl.ADR_fes < 1:
            r_s_new = np.min([r_fes, r_s_new])


        layer_thickness = r_upper - r_s_new

        #Regrid domain, keeping temperature relative to the adiabat.
        T_rel = T_new - prof.adiabat(r, Tcen, core_adiabat)

        r, T_rel = func.change_domain_size(T_rel, r, r_s_new, (resolution,min_res,max_res))

        T = T_rel + prof.adiabat(r, Tcen, core_adiabat)

        time_gone += dt_small

    if r_upper < r_cmb:
        #Append values in the FeS layer
        r = np.append(r, np.linspace(r_upper, prm.r_cmb,10)[1:])
        T = np.append(T, prof.adiabat(r[-9:], Tcen_fes, prm.fes_adiabat_params))

    #No layer this iteration
    if r_s_new == r_upper:
        layer_thickness = prm.FeS_size
        r_s_new = r_fes
        r, T = r[:9], T[:9]

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



#Overwrite existing conductivity function with this one to set a conductivity anomaly at the top of the core.
#Fes layer size must be included in conductivity params.

def conductivity_FeS(r, P, T, conductivity_params):
    
    if len(conductivity_params)==1:
        k = np.ones(r.size)*conductivity_params[0]
    
    else:

        k = np.ones(r.size)*conductivity_params[0]

        k[r>conductivity_params[1]] = conductivity_params[2]

    return k

prof.conductivity = conductivity_FeS