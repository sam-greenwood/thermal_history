import numpy as np
from ...utils.optimised_funcs import polyval

#Import redefined profiles functions
from . import profiles

#Define main functions based on leeds
from .. import leeds
update = leeds.update
progress = leeds.progress

def setup(model):

    prm = model.parameters
    core = model.core

    #Does not work with iron snow
    assert not prm.iron_snow, 'Method does not work with iron snow.'

    #Force initial stable layer to FeS layer size, then call leeds setup
    prm.layer_thickness = prm.FeS_size
    leeds.setup(model)

    #Set Tcen such that T_cmb is correct to specified value
    # taking into account the FeS layer
    r_fes = prm.r_cmb - prm.FeS_size
    core.Tcen =  core.T_cmb/polyval(prm.core_adiabat_params[::-1], r_fes)



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
    #so must be set here rather than in setup.
    if model.it == 1 and prm.stable_layer and prm.FeS_method == 'conducting':

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

        #Normalised secular cooling for FeS layer (normalised to Tcen).
        Qs_tilde = -rho_fes*cp_fes * (4/3)*np.pi*(prm.r_cmb**3 - r_fes**3) * (core.T_cmb/core.Tcen)

        #Normlised heatflow at r_fes before we modify it.
        Q_tilde = (core.Q_rs-core.Qr)/core.dT_dt

        Q_rs = core.Q_cmb * (1 - Qs_tilde/(Q_tilde + Qs_tilde))
        Q_rs = np.max([Q_rs, core.Q_rs]) #At least an adiabatic heat flow from the interior

        #Correct cooling rate for combined cooling of bulk and FeS layer
        dT_dt = Q_rs/Q_tilde

        #reset energies/entropies with new cooling rate
        for x in ['Qs', 'Ql', 'Qg', 'Es', 'El', 'Eg', 'Eg_mgo', 'dri_dt', 'dc_dt']:
            setattr(core, x, getattr(core, x)*(dT_dt/core.dT_dt))
        core.dT_dt = dT_dt

        #Reset Ej. Contribution from cooling FeS layer is added on when stable layer is evolved.
        core.Ej = core.Es + core.El + core.Eg + core.Eh + core.Er + core.Eg_mgo - core.Ek - core.Ea

        #Heat flow extracted from bulk
        core.Q_rs  = Q_rs #core.Q_cmb - Qs_tilde * core.dT_dt #* (core.T_cmb/core.Tcen)

        core.Q_fes = core.Q_cmb * (1 - Qs_tilde/(Q_tilde + Qs_tilde)) #heat flow at r_fes
        # core.Q_fes = Q_rs

        
        

    