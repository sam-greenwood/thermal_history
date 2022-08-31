Description = 'Mantle model based on Knibbe and van Westrenen (2018) for Mercury.\
Stagnant lid plate tectonics is assumed with growth of crust by mantle melting.'

#List individually confirmed compatibility with other regions.
compatibility = {'core': [''],
                 'stable_layer': [''],
                 'bmo': ['']}

#Import statements
import logging
logger = logging.getLogger(__name__)

import numpy as np
import copy as cp
from scipy.integrate import trapezoid

from .routines import profiles as prof
from .routines import melting as melt

#import matlab.engine
#eng = matlab.engine.start_matlab()

#setup, evolve, and update functions
def setup(model):
    '''
    Include any code that should be run after parameters have been loaded in and before 1st iteration is run.
    E.g. set any initial conditions as attributes of the model.
    '''
    mantle = model.mantle
    prm = model.parameters


    #Assign initial conditions
    mantle.Tm = cp.deepcopy(prm.Tm)
    mantle.delta_upper, mantle.delta_lower = 0, 0 #Initialise TBL thickness
    mantle.Dl = cp.deepcopy(prm.stagnant_lid_thickness)
    mantle.D_crust = prm.r_surf - prm.r_crust


    #Calculate the total budget of radiogenic heating at present. Used for calculating Qr at each timestep.

    #Mass fractions of K, U, Th
    c_K, c_U, c_Th = prm.HPE_concentrations #Present day K, U, Th concentrations (wt %) #1288e-6, 90e-6, 155e-6  #Tosi et al. (2013)

    #Fraction of radioactive nuclei per atom of element (K40, U238, U235, Th232)
    # Order of isotopes must be consistent with other properties. Numbers from Sramek et al. 2013
    isotope_abundance = np.array([117e-6, 0.9927, 0.007204, 1])

    #Crust mass fraction of radioactive isotopes
    HPE_mass_fraction = np.array([c_K, c_U, c_U, c_Th]) * isotope_abundance
    HPE_atom_mass = np.array([39.9640, 238.051, 235.044, 232.038])*1.661e-27 #Atomic masses
    HPE_Q_per_atom = np.array([0.110, 7.648, 7.108, 6.475])*1e-12 #Energy from entire decay chain, J/atom
    HPE_half_lives = np.array([1.25e9, 7.04e8, 4.47e9, 1.40e10])*prm.ys  #Half lives in seconds

    H = HPE_mass_fraction/HPE_atom_mass * HPE_Q_per_atom * np.log(2)/HPE_half_lives
    r_crust = prm.r_surf-38e3
    
    #FIX TO KW18 NUMBERS FOR NOW. heating rate @ 4.3 Ga
    H = np.array([1e-11*(1288/265), 7e-12*(90/20), 3.5e-12*(90/20), 2e-12*(155/69)])
    H = H/np.exp(np.log(2)*(4.3e9*prm.ys)/HPE_half_lives)

    #Present day total radiogenic heating assuming a 38km crust.
    Q_total =   H*prm.crust_density                  *(4/3)*np.pi*(prm.r_surf**3 - r_crust**3) \
              + H*prm.HPE_fraction*prm.mantle_density*(4/3)*np.pi*(r_crust**3    - prm.r_cmb**3)

    mantle.Qr_present = Q_total
    mantle.H_present  = H


def evolve(model):
    '''
    Calculates one timestep of the model.
    '''

    #Read in variables from model and parameters file
    mantle = model.mantle
    prm    = model.parameters


    Tm   = mantle.Tm      #average mantle temperature
    time = model.time

    r_surf  = prm.r_surf   #upper radius
    rc      = prm.r_cmb    #lower radius
    rl      = r_surf-mantle.Dl
    r_crust = r_surf-mantle.D_crust
    r_reg   = r_surf-prm.regolith_size

    g = prm.surface_gravity  #Surface gravity, assumed constant in mantle

    k_reg    = prm.regolith_conductivity
    k_crust  = prm.crust_k #thermal conductivity of the crust
    k_upper  = prm.mantle_k_upper  #thermal conductivity in upper mantle
    k_lower  = prm.mantle_k_lower  #thermal conductivity in lower mantle

    T_surf   = prm.T_surf                   #Surface temperature
    rho_man  = prm.mantle_density          #Mantle density
    rho_crust= prm.crust_density           #Crust density
    alpha    = prm.mantle_alpha_T          #volumetric expansion
    latent_heat = prm.mantle_latent_heat
    cp_crust = prm.crust_cp

    kappa    = prm.mantle_diffusivity_T    #thermal diffusivity
    cp       = prm.mantle_cp               #specific heat capacity
    Rg       = prm.Rg                      #Gas constant
    Av       = prm.activation_energy       #Activation energy
    V        = prm.activation_volume       #Activation volume
    eta_ref  = prm.viscosity_ref           #Reference viscosity
    T_ref    = prm.viscosity_ref_temperature #Temperature for reference viscosity
    b_factor = prm.b_factor                #'b' factor for scalings

    delta_upper, delta_lower = mantle.delta_upper, mantle.delta_lower
    rb = rc + delta_lower
    rm = rl - delta_upper

    #Get Tcmb
    T_cmb = model.core.T_cmb

    #Temperatures
    T_upper = Tm
    T_lower = Tm + (rl-delta_upper-rc-delta_lower)*Tm*alpha*g/cp

    T_lid = T_upper - 2.21 * prm.Rg*T_upper**2 / Av

    dT_upper  = T_upper-T_lid
    dT_lower  = T_cmb-T_lower
    if model.it==1:
        dT_lower = T_cmb-Tm

    Tup_av = (T_lid+T_upper)/2
    Tlow_av = (T_cmb+T_lower)/2

    #Pressures at base of lid, CMB and mid-mantle (assuming constant rho and g)
    r, rho = prof.radial_density_grid(rc, rb, rm, rl, r_crust, r_reg, r_surf, rho_man, rho_crust)
    P = prof.pressure(r, rho, g)
    P_upper = np.interp(rl-delta_upper, r, P)
    P_lower = np.interp(rc+delta_lower, r, P)

    #Viscosities in the mantle. Scales eta_0 with activation energy (relative to T_ref) and volume (relative to 0 Pa).
    eta        = eta_ref * np.exp((Av+np.mean(P)*V)/Rg  * (1/Tm-1/T_ref))
    eta_upper  = eta_ref * np.exp((Av+P_upper*V)/Rg * (1/Tup_av -1/T_ref))
    eta_lower  = eta_ref * np.exp((Av+P_lower*V)/Rg * (1/Tlow_av-1/T_ref))

    # Rayleigh numbers.
    Rai = g*alpha*rho_man*(np.abs(Tm-T_surf)+np.abs(dT_lower))*(r_surf-rc)**3 /(kappa*eta)
    Rai_c = 0.28*Rai**0.21

    Ra = g*alpha*rho_man*(np.abs(dT_upper)+np.abs(dT_lower))*(rl-rc)**3 /(kappa*eta)
    if T_cmb-T_lower+T_upper-T_lid < 10:
        # breakpoint()
        Ra = g*alpha*rho_man*(10)*(rl-rc)**3 /(kappa*eta)
    Ra_c = prm.critical_rayleigh

    mantle.test = np.array([dT_lower, dT_upper])

    # breakpoint()

    #TBL sizes. TESTING
    delta_upper = (rl-rc)*(Ra_c/Ra)**b_factor

    # breakpoint()
    if np.abs(T_cmb-T_lower)>5:
        delta_lower = (kappa*eta_lower*Rai_c / (alpha*rho_man*g*(np.abs(dT_lower))))**b_factor #Use abs(dT_lower) otherwise it is not defined
    else:
        delta_lower = mantle.delta_lower #Use previous time steps value

    #Set boundaries the same for first timestep
    if model.it==1:
        delta_lower = delta_upper

    #My Way
    if rl-delta_upper < rc+delta_lower:
        delta_upper = rl - (rc+delta_lower)
    # #KW way
    # if rl-delta_upper < rc:
    #     delta_upper = rl - rc

    #Recalculate T_lower, rb/rm now deltas have been calculated
    T_lower = Tm + (rl-delta_upper-rc-delta_lower)*Tm*alpha*g/cp
    rb = rc + delta_lower
    rm = rl - delta_upper
    

    if np.isnan(delta_lower) or np.isnan(delta_upper):
        breakpoint()


    #Fluxes through TBL
    dT_upper  = T_upper-T_lid
    dT_lower  = T_cmb-T_lower

    F_m = k_upper * dT_upper / delta_upper
    F_cmb = k_lower * dT_lower / delta_lower
    if np.abs(dT_lower) < 5:
        F_cmb = 0

    # r_crust = r_surf-38e3 #Fixed 38km thick crust
    mass_crust  = 4/3 * np.pi * (r_surf**3-r_crust**3) * rho_crust
    mass_mantle = 4/3 * np.pi * (r_crust**3-rc**3) * rho_man

    #Radiogenic heating
    age  = (4.3e9*prm.ys - time)

    #Half lives in seconds. Must be the same as defined in setup() above
    HPE_half_lives = np.array([1.25e9, 7.04e8, 4.47e9, 1.40e10])*prm.ys

    #Total amount of radiogenic heating in the mantle+crust
    Q_K, Q_U235, Q_U238, Q_Th = mantle.Qr_present * np.exp(np.log(2)*age/HPE_half_lives)

    #Crustal heating rate, W/kg
    H = np.sum(mantle.H_present * np.exp(np.log(2)*age/HPE_half_lives))

    #Mantle heating rate, W/kg
    H_man = (np.sum([Q_K, Q_U235, Q_U238, Q_Th]) - H * rho_crust * (4/3)*np.pi*(r_surf**3-r_crust**3))/mass_mantle

    #Radiogenic heating totals
    mantle.Qr_crust = mass_crust  * H
    mantle.Qr       = rho_man * 4/3 * np.pi * (rl**3-rc**3) * H_man


    #Triple conductive layers
    solution = joint_conduction(T_surf, T_lid,
                                [r_surf, r_reg, r_crust,         rl],
                                [H*rho_crust, H*rho_crust, H_man*rho_man],
                                [k_reg,         k_crust,         k_upper])

    sol_reg   = solution[0]
    sol_crust = solution[1]
    sol_lith  = solution[2]

    #Heat fluxes through lid and surface
    F_lid = -k_upper*(-sol_lith[0]/rl**2 - 2*rl*(H_man*rho_man) / (6*k_upper))

    F_surf = -k_reg*(-sol_reg[0]/r_surf**2 - 2*r_surf*(H*rho_crust) / (6*k_reg))

    #Radial profiles
    r, rho = prof.radial_density_grid(rc, rb, rm, rl, r_crust, r_reg, r_surf, rho_man, rho_crust)
    T = np.zeros(r.size)
    
    #First 4 nodes are interfaces rc->rl (assume linear between)
    T[:4] = np.array([T_cmb, T_lower, Tm, T_lid])
    
    #Conductive layers
    T[4:] =  sol_lith[0]/r[4:] + sol_lith[1] - H_man*rho_man*r[4:]**2/(6*k_upper) #Lithosphere

    T[4+10:] =  sol_crust[0]/r[4+10:] + sol_crust[1] - H*rho_crust*r[4+10:]**2/(6*k_crust) #Crust

    T[4+20:] =  sol_reg[0]/r[4+20:] + sol_reg[1] - H*rho_crust*r[4+20:]**2/(6*k_reg) #Regolith

    # T[4:] = T_lid - (r[4:]-rl)*(T_lid-T_surf)/(r_surf-rl)

    #Mantle Melting
    D_ref = (0.2/3)*(prm.r_surf**3 - prm.r_cmb**3)/prm.r_surf**2

    #Calculate volume of mantle melting and volumetric average melt fraction
    V_melt, ma, melting_radii, dm_dT, profiles = melt.mantle_melt(D_ref, mantle.D_crust, r, rho, g, T, prm.mantle_solidus_params, prm.mantle_liquidus_params)

    mantle.profiles = profiles

    mantle.melting_radii = melting_radii

    if V_melt > 0:

        Vm     = (4/3)*np.pi*(rl**3 - rc**3)

        #Convection velocity
        U = prm.convection_velocity * (Ra/Ra_c)**(2*1/3)

        #Crustal growth
        d_crust_dt = U * ma * V_melt / (4*np.pi*r_surf**3)

        St = latent_heat*V_melt*dm_dT / (cp*Vm)

        #Thermal energy associated with crustal growth.
        crust_term = (rho_crust*latent_heat + rho_crust*cp_crust*(T_upper-T_lid))*d_crust_dt


    else:
        d_crust_dt = 0
        crust_term = 0
        St = 0

    mantle.dDl_dt = (F_lid-F_m + crust_term)/(rho_man*cp*(T_upper-T_lid)) #Lithosphere growth


    #Force crust to be thinner that lithosphere
    if mantle.D_crust + d_crust_dt*model.dt > mantle.Dl + mantle.dDl_dt*model.dt:
        d_crust_dt = ((mantle.Dl + mantle.dDl_dt*model.dt -100) - mantle.D_crust)/model.dt
        crust_term = (rho_crust*latent_heat + rho_crust*cp_crust*(T_upper-T_lid))*d_crust_dt


    #Mantle energy budget
    mantle.Q_surf = 4*np.pi*r_surf**2 * F_surf
    mantle.Q_cmb  = 4*np.pi*rc**2 * F_cmb
    mantle.Q_conv = 4*np.pi*rl**2 * F_m
    mantle.Q_lid  = 4*np.pi*rl**2 * F_lid
    mantle.Q_vol  = 4*np.pi*rl**2 * crust_term

    epsilon = Tm/T_upper
    epsilon = 1 #KW keep this as 1
    mass = (4*np.pi/3)*(rl**3-rc**3)*rho_man
    mantle.dT_dt = (mantle.Q_cmb - mantle.Q_conv - mantle.Q_vol + mantle.Qr)/(cp*mass*(1+St)*epsilon)

    #Save variables
    mantle.dD_crust_dt = d_crust_dt
    mantle.St = St

    mantle.mass = mass

    mantle.T_upper = T_upper
    mantle.T_lower = T_lower

    mantle.delta_upper = delta_upper
    mantle.delta_lower = delta_lower

    mantle.Ra = Ra
    mantle.Rai = Rai

    mantle.Ra_c = Ra_c
    mantle.Rai_c = Rai_c

    mantle.dT_upper = dT_upper
    mantle.dT_lower = dT_lower

    mantle.eta_upper = eta_upper
    mantle.eta_lower = eta_lower

    mantle.Qs = mantle.Q_cmb - mantle.Q_surf + mantle.Qr

    mantle.T_cmb  = T_cmb


def update(model):
    '''
    Update variables based on what was calculated in the evolve() function.
    '''
    mantle = model.mantle
    prm = model.parameters

    mantle.Tm += mantle.dT_dt*model.dt
    mantle.Dl += mantle.dDl_dt*model.dt
    mantle.D_crust += mantle.dD_crust_dt*model.dt

    if mantle.Dl == 0:
        raise ValueError('No Lithosphere')

    if mantle.Dl > (prm.r_surf-prm.r_cmb+mantle.delta_lower):
        mantle.Dl = prm.r_surf-prm.r_cmb+mantle.delta_lower

#Required parameters. Used to check given parameters before start. Give the name of the parameter as the key with a description as the associated value within the dictionary.

required_params = {'Tm': 'Initial temperature of the mid mantle [K]. Float',
                   'T_surf': 'Surface temperature [K]. Float',
                   'r_surf': 'Radius of the planet [m]. Float',
                   'r_crust': 'Radius of the base of the crust [m]. Float',
                   'stagnant_lid_thickness': 'Initial thickness of the stangant lid [m]. Float',
                   'mantle_density': 'Mantle density [kg/m^3]. Float',
                   'crust_density': 'Crust density [kg/m^3]. Float',
                   'mantle_k_upper': 'Thermal conductivity of the upper mantle [W/K/m]. Float',
                   'mantle_k_lower': 'Thermal conductivity of the lower mantle [W/K/m]. Float',
                   'crust_k': 'Thermal conductivity of the crust [W/K/m]. Float',
                   'mantle_alpha_T': 'Thermal expansivity of the mantle [/K]. Float',
                   'mantle_cp': 'Mantle specific heat capacity [J/kg/K]. Float',
                   'crust_cp': 'Erupted magma specific heat capacity [J/kg/K]. Float',
                   'mantle_diffusivity_T': 'Mantle thermal diffusivity [m^2/s]. Float',
                   'activation_energy': 'Activation energy for viscosity temperature scaling. Float',
                   'activation_volume': 'Activation volume for viscosity pressure scaling [m^3/mol]. Float',
                   'viscosity_ref': 'Reference viscosity at a specified temperature and atomospheric pressure [Pa s]. Float',
                   'viscosity_ref_temperature': 'Temperature of the reference viscosity [K]. Float',
                   'b_factor': 'Exponent used in scaling the thermal boundary layer thickness. Float',
                   'HPE_fraction': 'Fraction of heat producing elements in the bulk mantle relative to the crust. Float.',
                   'HPE_concentrations': 'Mass fractions of heat producing elements present in the crust in the order [K, U, Th]. List[Float]',
                   'regolith_conductivity': 'Thermal conductivity of the surface regolith. Float',
                   'regolith_size': 'Size of the regolith layer [m]. Float',
                   'surface_gravity': 'Gravity at the surface, only used if not using a core model.',
                   'critical_rayleigh': 'Critical Rayleigh number',
                   'mantle_solidus_params': 'Mantle solidus pressure polynomials. List',
                   'mantle_liquidus_params': 'Mantle liquidus pressure polynomials. List',
                   'convection_velocity': 'Typical convection velocity for scaling volume of erupted magma [m/s]. Float',
                   'mantle_latent_heat': 'Latent heat of melting for mantle [J/kg]'
}


def joint_conduction(T_upper, T_lower, r, H, k):
    '''
    Steady state spherical conduction solution for n shells.

    ---------- r0, T_upper
    k1, H1
    ---------- r1
    k2, H2
    ---------- r2
    ...

    kn, Hn
    ---------- rn, T_lower

    Shells are defined by an array of size n for their thermal conducitivities (W/m/K) and
    internal volumetric heating rates (W/m^3). Radii are specified by an array of size n+1
    of the boundaries and interfaces. Index 0 of arrays refer to the values of the outer-most shell, proceeding
    to the inner-most shell. Temperature boundary conditions for the radii r0 and rn are required, with continuity
    of temperature and heat flux assumed at each interface.

    Returns an array of coefficients (A1,B1...An,Bn) to the solution

    T(r) = Ai/r + Bi - Hi/(6 * ki) * r**2

    for each shell where i is in the range [1,n].

    ______________________________________
    Example for 2 shells

    #Temperature Boundary conditions
    T_upper=1
    T_lower=1

    #Radii of 2 shells, outer shell between r=0.5:1, inner shell between r=0:0.5
    r = [1,0.5,0]

    #Internat heating rates and conductivities for outer and inner shells.
    H = [0,0]
    k = [1,1]

    solution = joint_conduction(T_upper, T_lower, r, H, k)
    A1, B1 = solution[:2]
    A2, B2 = solution[2:]
    '''

    r, H, k = np.array(r), np.array(H), np.array(k)

    #Define variables l1/2 for simplicity
    L = 1/(6*k)

    #Initialise sizes of A matrix and b vector
    n = k.size*2
    A = np.zeros((n,n))
    b = np.zeros(n)

    #Boundary conditions at r[0] and r[-1]
    A[0,:2]    = [1/r[0] , 1]
    A[-1, -2:] = [1/r[-1], 1]

    b[0]  = T_upper + H[0]  * L[0]  * r[0]**2
    b[-1] = T_lower + H[-1] * L[-1] * r[-1]**2

    #Fixed flux and temperature at interfaces r[1:-1]
    i = 1
    for c in range(k.size-1):
        j = c*2
        w = j+4

        A[i,   j:w] = [-1/r[c+1], -1, 1/r[c+1], 1] #Temperature continuity
        A[i+1, j:w] = [k[c],       0, -k[c+1],  0] #Flux continuity

        b[i]   = (H[c+1]*L[c+1] - H[c]*L[c]) * r[c+1]**2 #Temperature continuity
        b[i+1] = (1/3)*r[c+1]**3 * (H[c+1]-H[c])         #Flux continuity

        i+=2

    #Solve the system for [A1,B1,...,An,Bn]
    solution = np.linalg.solve(A,b)

    #Return list where 0th index gives array(A1,B1) and (n-1) index gives array(An,Bn)
    return [solution[i:i+2] for i in range(0,len(solution)-1,2)]



# def joint_conduction_old(T0, T2, r0, r1, r2, H1, H2, k1, k2):
#     '''
#     Steady state spherical conduction solution for n shells.
#
#     ---------- r0, T_upper
#     k1, H1
#     ---------- r1
#     k2, H2
#     ---------- r2
#     ...
#
#     kn, Hn
#     ---------- rn, T_lower
#
#     Shells are defined by an array of size n for their thermal conducitivities (W/m/K) and
#     internal volumetric heating rates (W/m^3). Radii are specified by an array of size n+1
#     of the boundaries and interfaces. Index 0 of arrays refer to the values of the outer-most shell, proceeding
#     to the inner-most shell. Temperature boundary conditions for the radii r0 and rn are required, with continuity
#     of temperature and heat flux assumed at each interface.
#
#     Returns an array of coefficients (A1,B1...An,Bn) to the solution
#
#     T(r) = Ai/r + Bi - Hi/(6 * ki) * r**2
#
#     for each shell where i is in the range [1,n].
#
#     ______________________________________
#     Example for 2 shells
#
#     #Temperature Boundary conditions
#     T_upper=1
#     T_lower=1
#
#     #Radii of 2 shells, outer shell between r=0.5:1, inner shell between r=0:0.5
#     r = [1,0.5,0]
#
#     #Internat heating rates and conductivities for outer and inner shells.
#     H = [0,0]
#     k = [1,1]
#
#     solution = joint_conduction(T_upper, T_lower, r, H, k)
#     A1, B1 = solution[:2]
#     A2, B2 = solution[2:]
#
#
#
#     '''
#
#     #Define variables l1/2 for simplicity
#     l1, l2 = 1/(6*k1), 1/(6*k2)
#
#     #Set up linear equations in form Ax=b
#     b = np.array([T0 + H1 * l1 * r0**2,
#                   T2 + H2 * l2 * r2**2,
#                   (H2*l2 - H1*l1) * r1**2,
#                   2*r1**3 * (l1*H1*(k1/k2) - l2*H2)])
#
#     A = np.array([[1/r0, 1, 0, 0],
#                   [0, 0, 1/r2, 1],
#                   [-1/r1, -1, 1/r1, 1],
#                   [-k1/k2, 0, 1, 0]])
#
#
#     # Ac, Bc, Al, Bl = np.linalg.solve(A,b)
#
#     # rc = np.linspace(r1,r0)
#     # rl = np.linspace(r2,r1)
#     # Tc = Ac/rc + Bc - H1*l1*rc**2
#     # Tl = Al/rl + Bl - H2*l2*rl**2
#     # plt.plot(rl,Tl)
#     # plt.plot(rc,Tc)
#     # plt.show()
#     # breakpoint()
#
#     #Solve the system for [A1,B1,A2,B2]
#     return np.linalg.solve(A,b)
