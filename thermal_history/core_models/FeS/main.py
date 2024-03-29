Description = 'Based on leeds model but with support for a liquid FeS layer. This is under development!'

import numpy as np
from ...utils.optimised_funcs import polyval, trapezoid

#Import redefined profiles functions. Need to at least import it to run the code
#and refefine functions from the leeds model.
from . import profiles

#Define main functions based on leeds
from .. import leeds
update = leeds.update
progress = leeds.progress

def setup(model):

    prm = model.parameters

    #Does not work with iron snow
    assert not prm.iron_snow, 'Method does not work with iron snow.'

    #Call leeds setup
    leeds.setup(model)

    #Set stable layer radius to be bottom of FeS layer
    prm.layer_thickness                = prm.FeS_size
    model.core.rs                      = prm.r_cmb - prm.FeS_size
    model.stable_layer.layer_thickness = prm.FeS_size

#Add on extra required parameters
required_params = leeds.required_params
required_params['FeS_size']    = 'Size of the FeS layer'
required_params['FeS_density'] = 'Uniform density of liquid FeS'
required_params['FeS_cp']      = 'Uniform specific heat of liquid FeS'

optional_params = leeds.optional_params


def evolve(model):

    core = model.core
    sl = model.stable_layer
    prm = model.parameters

    #For the conducting method, initial thermal structure depends on Q_cmb
    #so must be set here rather than in setup. This resets the bulk temperature
    #to be consistent with T_cmb and Q_cmb.
    if model.it == 1 and prm.stable_layer and prm.sl_method_name == 'FeS_conducting':

        #Initial temp in FeS layer is steady state profile fit to T_cmb and Q_cmb

        T_grad_cmb = model.mantle.Q_cmb/(-model.core.profiles['k'][-1]*4*np.pi*prm.r_cmb**2)

        r_prof_fes = sl.profiles['fes']['r']
        r_fes = r_prof_fes[0]        

        A = prm.r_cmb**2 * T_grad_cmb
        B = core.T_cmb + (A/prm.r_cmb)

        #Steady state temperature
        sl.profiles['fes']['T'] = -A/sl.profiles['fes']['r'] + B

        #Reset bulk core temperture
        core.Tcen = sl.profiles['fes']['T'][0]/polyval(prm.core_adiabat_params[::-1], r_fes)

        #Reset profiles with new temperature profile
        leeds.routines.profiles.basic_profiles(model)
        leeds.routines.profiles.temp_dependent_profiles(model)

        #Conductive heat flow from bulk to FeS layer
        core.Q_fes = -4*np.pi*A * core.profiles['k'][core._fes_idx]


    #Run standard leeds core model evolution
    leeds.evolve(model)

    #Account for the fact that it treats the FeS layer as a stable layer
    #and so expects an adiabatic heat flow at rs. Not necesasrily the case when no
    #thermal stratification exists deeper than the FeS layer. Need to adjust
    #cooling rate in this scenario.

    if sl.layer_thickness == prm.FeS_size:

        r_fes = prm.r_cmb - prm.FeS_size
        
        #As core cools, temp at r_fes is continuous between bulk and FeS layer.

        rho_fes = prm.FeS_density
        cp_fes  = prm.FeS_cp

        #Normalised secular cooling for an isothermal FeS layer (normalised to Tcen).
        #Not used if stable layer solution is included.
        Qs_tilde = -rho_fes*cp_fes * (4/3)*np.pi*(prm.r_cmb**3 - r_fes**3) * (core.T_cmb/core.Tcen)

        if prm.stable_layer:
            r_prof_fes, T_prof_fes = (model.stable_layer.profiles['fes'][x] for x in ['r', 'T'])
            Qs_tilde = -4*np.pi*trapezoid(r_prof_fes, prm.FeS_density * prm.FeS_cp * (T_prof_fes/T_prof_fes[0]) * r_prof_fes**2)[-1] * (T_prof_fes[0]/core.Tcen) #Normalised to Tcen

        #Normlised heatflow at r_fes before we modify it.
        Q_tilde = (core.Q_rs-core.Qr)/core.dT_dt

        Q_rs = core.Q_cmb * (1 - Qs_tilde/(Q_tilde + Qs_tilde))

        #If a stable layer is being used, then heat flow is always at least adiabatic at rs
        #Or that defined by conductive profile.
        if prm.stable_layer:
            Q_rs = np.max([core.Q_fes, core.Qa_rs])

        #Correct cooling rate for combined cooling of bulk and FeS layer
        dT_dt = Q_rs/Q_tilde

        #reset energies/entropies with new cooling rate
        for x in ['Qs', 'Ql', 'Qg', 'Es', 'El', 'Eg', 'Eg_mgo', 'dri_dt', 'dc_dt']:
            setattr(core, x, getattr(core, x)*(dT_dt/core.dT_dt))
        core.dT_dt = dT_dt

        #Reset Ej. Contribution from cooling FeS layer is added on when stable layer evolve method is called.
        core.Ej = core.Es + core.El + core.Eg + core.Eh + core.Er + core.Eg_mgo - core.Ek - core.Ea

        #Heat flow extracted from bulk
        core.Q_rs  = Q_rs #core.Q_cmb - Qs_tilde * core.dT_dt #* (core.T_cmb/core.Tcen)

        core.ADR_s = core.Q_rs/core.Qa_rs #Reset ADiabat Ratio.

        core.Q_fes = core.Q_cmb * (1 - Qs_tilde/(Q_tilde + Qs_tilde)) #heat flow at r_fes

        #If no stable layer solution add on secular cooling of FeS layer
        if not prm.stable_layer:
            core.Qs += Qs_tilde * core.dT_dt* (core.T_cmb/core.Tcen) #Secular cooling of FeS layer.
            #No need to add Ek and Es as both are zero for an isothermal FeS layer (which is assumed if no stable layer method).


        

def set_Q_rs(model):
    '''
    Sets the heat flow at the top of the convecting region (rs).
    Modified to account for FeS layer.

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

    r_fes = prm.r_cmb - prm.FeS_size

    #Adiabatic heat flow.
    Q_rs = profiles.adiabatic_heat_flow(core.rs, k[rs_idx-1], core.Tcen, prm.core_adiabat_params)

    if prm.stable_layer:
        Q_rs = core.Q_fes

    return Q_rs
leeds.set_Q_rs = set_Q_rs