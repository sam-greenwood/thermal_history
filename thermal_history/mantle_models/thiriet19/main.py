Description = 'Mantle model for stagnant lid convection adapted from the Mars model in Thiriet et al. (2019). \
Here stagnant lid size is assumed consant in time'

from copy import deepcopy
import numpy as np

def setup(self):

    mantle = self.mantle
    prm = self.parameters

    mantle.Tm = deepcopy(prm.Tm)  #Average mantle temperature

    #Set stagnant lid size
    mantle.lid_thickness = prm.stagnant_lid_thickness


def evolve(self):

    #Read in variables from model and parameters file
    mantle = self.mantle
    prm = self.parameters


    Tm   = mantle.Tm      #average mantle temperature
    time = self.time
    it   = self.it
    dt   = self.dt

    r_surf  = prm.r_surf   #upper radius
    r_cmb  = prm.r_cmb    #lower radius


    k_upper  = prm.mantle_k_upper  #thermal conductivity in upper mantle
    k_lower  = prm.mantle_k_lower  #thermal conductivity in lower mantle

    T_surf   = prm.T_surf                   #Surface temperature
    rho      = prm.mantle_density          #Mantle density
    alpha    = prm.mantle_alpha_T          #volumetric expansion
    g        = prm.g                       #gravity

    #Get gravity from core model if it exists.
    if prm.core:
        g_cmb    = self.core.profiles['g'][-1]
    else:
        g_cmb = g

    kappa    = prm.mantle_diffusivity_T    #thermal diffusivity
    cp       = prm.mantle_cp               #specific heat capacity
    Rac      = prm.Rac                     #critical Rayleigh Number (upper)
    Rg       = prm.Rg                      #Gas constant
    Av       = prm.Av                      #Activation energy
    V        = prm.activation_volume       #Activation volume
    eta_ref  = prm.viscosity_ref           #Reference viscosity
    T_ref    = prm.viscosity_ref_temperature #Temperature for reference viscosity
    a_factor = prm.a_factor                #'a' factor for scalings
    b_factor = prm.b_factor                #'b' factor for scalings

    #Get Tcmb and lid size
    T_cmb         = self.core.T_cmb
    lid_thickness = mantle.lid_thickness

    r_lid = r_surf - lid_thickness   #radius of stagnant lid

    mass = (4*np.pi/3)*(r_lid**3-r_cmb**3)*rho

    T_upper = Tm - a_factor * Rg*Tm**2/Av     #Temperature of upper convecting mantle

    T_lower     = Tm #+ alpha_m*Tm*g*dR/cp_m      #T drop adiabatic increase, see below14 #Temperature of lower convecting mantle
    T_lower_av  = (T_cmb+Tm)/2    #Average temperature in the lower boundary layer for setting the viscosity

    dT_upper  = Tm-T_upper #Chk - their stuff below 14 cant be right!  #I think they got the two equations for dTu the wrong way round. Meant to say it is =Tm-Tl which can be approximated as Tc-Tl.
    dT_lower  = T_cmb-T_lower

    #Pressures
    P_upper = rho*g*lid_thickness
    P_lower = rho*np.mean([g,g_cmb])*(r_surf-r_cmb)

    #Mantle Viscosity
    eta_upper  = eta_ref * np.exp((Av/Rg)*(1/Tm -1/T_ref)        + (P_upper*V)/(Tm*Rg))
    eta_lower  = eta_ref * np.exp((Av/Rg)*(1/T_lower_av-1/T_ref) + (P_lower*V)/(T_lower_av*Rg))

    # Assumes g and alpha and kappa are constant
    Ra_upper   = alpha * rho * g     * dT_upper * (r_lid-r_cmb)**3 / (kappa * eta_upper)  #Upper and lower Rayleigh numbers
    Ra_lower   = alpha * rho * g_cmb * dT_lower * (r_lid-r_cmb)**3 / (kappa * eta_lower)

    #Must have positive Rayleigh number
    if Ra_lower < 0:
        Ra_lower = -Ra_lower

    Rac_lower = 0.28*(alpha*rho*g*(T_cmb-T_surf)*(r_surf-r_cmb)**3 / (kappa * eta_lower))**0.21
    # Rac_lower = Rac*eta_upper/eta_lower
    Rac_upper = Rac
    # Thiriet never say where eta is evaluated in Ra_int
    # I also assume dT and D here refer to whole shell.

    delta_upper = (r_lid-r_cmb)*(Rac_upper/Ra_upper)**b_factor
    delta_lower = (r_lid-r_cmb)*(Rac_lower/Ra_lower)**b_factor


    # k = k_upper #Same conductivity for lower/upper
    F_lid = k_upper * dT_upper / delta_upper
    F_cmb = k_lower * dT_lower / delta_lower

    #Radiogenic heating
    age  = (4.5e9*prm.ys - time)/prm.ys

    hlK40, hlU238, hlU235, hlT232 = 1.26e9, 4.47e9, 7.04e8, 1.40e10  #Half lives, in years

    #Mars heating rates as used in Thiriet et al. (2019)
    mantle.Qr = (1.9*5.56e-13  *np.exp(np.log(2)*age/hlK40)+
                 1.50e-12      *np.exp(np.log(2)*age/hlU238)+
                 6.46e-14      *np.exp(np.log(2)*age/hlU235)+
                 0.875*1.69e-12*np.exp(np.log(2)*age/hlT232))*mass

    mantle.Q_surf = 4*np.pi*r_lid**2 * F_lid
    mantle.Q_cmb  = 4*np.pi*r_cmb**2 * F_cmb

    mantle.dT_dt = (mantle.Q_cmb - mantle.Q_surf + mantle.Qr)/(cp*mass) #GB08 eq 1


    mantle.mass = mass

    mantle.T_upper = T_upper
    mantle.T_lower = T_lower

    mantle.delta_upper = delta_upper
    mantle.delta_lower = delta_lower

    mantle.Ra_upper = Ra_upper
    mantle.Ra_lower = Ra_lower

    mantle.Rac_upper = Rac_upper
    mantle.Rac_lower = Rac_lower

    mantle.dT_upper = dT_upper
    mantle.dT_lower = dT_lower

    mantle.eta_upper = eta_upper
    mantle.eta_lower = eta_lower

    mantle.Qs = mantle.Q_cmb - mantle.Q_surf + mantle.Qr

    mantle.T_cmb  = T_cmb


def update(self):
    dt = self.dt
    mantle = self.mantle
    mantle.Tm += mantle.dT_dt * dt




#Dictionary of required parameters by this model. 'Name': 'Description'
required_params = {'Tm': 'Initial mantle temperature, float',
                   'r_surf': 'Surface radius, float',
                   'r_cmb': 'CMB radius, float',
                   'stagnant_lid_thickness': 'Thickness of stagnant lid, float',
                   'mantle_k_upper': 'Upper mantle thermal conducitivity, float',
                   'mantle_k_lower': 'Lower mantle thermal conductivity',
                   'T_surf': 'Surface temperature, float',
                   'mantle_density': 'Mantle density, float',
                   'mantle_alpha_T': 'Mantle thermal expansivity, float',
                   'g': 'Surface gravity, float',
                   'mantle_diffusivity_T': 'Mantle thermal diffusivity, float',
                   'mantle_cp': 'mantle specific heat capacity, float',
                   'Rac': 'Critical Rayleigh number for mantle convction, float',
                   'Rg': 'Gas constant',
                   'Av': 'Activation energy for mantle viscosity, float',
                   'viscosity_ref': 'Reference mantle viscosity, float',
                   'viscosity_ref_temperature': 'Temperature of reference viscosity, float',
                   'a_factor': 'Factor controlling temperature drop across upper thermal boundary layer (Eq 15 of Thiriet et al. 2019)',
                   'b_factor': 'Factor scaling boundary layer thickness with Rayleigh number.',
                   'activation_volume': 'Activate volume for mantle viscosity, float'
                   }


#Optional parameters that the user may omit, in which event are set to the specified default values in this dictionary.
#Use the format {'Name': ('Description', default_value)}
optional_params = {}


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

    text = f'    Tm: {mantle.Tm:.2f} ˚K'

    return text

