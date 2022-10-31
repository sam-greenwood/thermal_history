'''New function for the melting curve to use in place of the one in the 'leeds' core model.
'''

#Core chemistry functions

from thermal_history.utils.optimised_funcs import polyval
import numpy as np

from ...leeds.routines.chemistry import *  #Import all functions from the leeds core model. Then redefine melting curve.


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

    'RI': Uses a polynomial fit to Rivoldini's EOS model. Polynomial coefficients are for both Pressure and Compositional dependence
    and are calculated in the included stand-alone script riv_pre_compute_poly.py. 

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


    elif melting_params[0] == 'RI':
        #Melting temperature based on fitting Rivoldini's EOS model to represent
        #Tm with polynomials in both pressure and composition. Polynomials are pre-computed
        #using the riv_pre_compute_poly.py script in the routines sub-package. Output coefficients
        #are saved into optimal_solution.txt and must be included into the parameters file.

        #Simple partitioning
        core.conc_s = core.conc_l * prm.partition_coeff

        params = melting_params[1:].astype('float64') #Don't use first string

        #melting points of iron and alloy
        Tm_fe = riv_Tm(P, 0, params)
        Tm = riv_Tm(P, core.conc_l[0], params)
        dTm = Tm-Tm_fe

        dTm_dP = riv_dTm_dP(P, core.conc_l[0], params)

    elif melting_params[0] == 'RU':
        #Melting temperature based on Rückriemen et al. (2018). See their apendix A2.

        #Simple partitioning
        core.conc_s = 0.3647 #FeS

        params = melting_params[1:].astype('float64') #Don't use first string

        T_eutectic = params[0] - params[1]*(P-params[2])
        dT_e_dP = -params[1]

        Tm_fes = params[3] + params[4]*P - params[5]*P**2
        dTm_fes_dP = params[4] - 2*params[5]*P

        x_eu = 0.11 + 0.187*np.exp(-0.065e-9*P)
        dx_dP = -0.065e-9*0.187*np.exp(-0.065e-9*P)

        Tm = Tm_fes + (Tm_fes-T_eutectic)/(0.3647-x_eu)  * (core.conc_l - 0.3647)

        dTm_dP = dTm_fes_dP + ((0.3647-x_eu)*dTm_fes_dP - Tm_fes*dx_dP)/((0.3647-x_eu)**2) - ((0.3647-x_eu)*dT_e_dP - T_eutectic*dx_dP)/((0.3647-x_eu)**2)

        for j in range(dTm_dP.size-1):
            if j > 1 and core.profiles['r'][j] == core.profiles['r'][j+1]:
                dTm_dP[j] = dTm_dP[j-1]
            else:
                dTm_dP[j] = (Tm[j+1]-Tm[j])/(P[j+1]-P[j])
        dTm_dP[j+1] = dTm_dP[j]

        Tm_fe = (Tm_fes + (Tm_fes-T_eutectic)/(0.3647-x_eu)  * (0-0.3647))
        dTm = Tm - Tm_fe


    else:

        raise ValueError(f'Incorrect string denoting which melting parameterisation to use. See this function (core_models.leeds.routines.chemisty.melting_curve) for valid options. Supplied core_melting_params={melting_params}')

    Tm = Tm_fe + dTm

    return Tm_fe, Tm, dTm_dP

