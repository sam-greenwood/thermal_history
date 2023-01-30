Description = 'Based on leeds model but with support for a liquid FeS layer. This is under development!'

import numpy as np
from ...utils.optimised_funcs import polyval, trapezoid

#Import redefined profiles functions. Need to at least import it to run the code
#and refefine functions from the leeds model.
from . import profiles
from ..leeds.routines import energy as en

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

    #*# Changes associated with adapting this for a convecting FeS layer are marked with #*#

    #*# removed setup for conducting region, not needed. Assumes Ta is the same in bulk and FeS region.

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

        #*# Integrate across FeS layer to get normailsed secular cooling
        profiles = core.profiles
        idx      = core._fes_idx
        r, Ta, Ta_grad, k, cp, = (profiles[key] for key in ['r','Ta','Ta_grad','k','cp'])
        
        #*# Normalised secular cooling for FeS layer (normalised to Tcen).
        Qs_tilde_core, Es_tilde_core = en.secular_cool(r, rho_fes, Ta, cp, r.size)
        Qs_tilde_bulk, Es_tilde_bulk = en.secular_cool(r, rho_fes, Ta, cp, idx)
        Qs_tilde_fes, Es_tilde_fes = Qs_tilde_core-Qs_tilde_bulk,  Es_tilde_core-Es_tilde_bulk

        # Qs_tilde = -rho_fes*cp_fes * (4/3)*np.pi*(prm.r_cmb**3 - r_fes**3) * (core.T_cmb/core.Tcen)

        #Normlised heatflow at r_fes before we modify it.
        Q_tilde = (core.Q_rs-core.Qr)/core.dT_dt

        Q_rs = core.Q_cmb * (1 - Qs_tilde_fes/(Q_tilde + Qs_tilde_fes))

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

        #If no stable layer solution, Q_fes is set here. Also add on secular cooling of FeS layer.
        #These are handled in stable layer model instead if being used.
        if not prm.stable_layer:
            core.Q_fes = core.Q_cmb * (1 - Qs_tilde_fes/(Q_tilde + Qs_tilde_fes)) #heat flow at r_fes
            core.Qs += Qs_tilde_fes * core.dT_dt* (core.T_cmb/core.Tcen) #Secular cooling of FeS layer.

            #!# Add on Ek and Es contributions from FeS layer
            core.Ek += en.cond_entropy(Ta[idx:], Ta_grad[idx:], k[idx:], r[idx:], r.size-idx)
            core.Es += Es_tilde_fes*dT_dt


        

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