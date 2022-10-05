#Automatically generated parameters file for the following methods:
#core: leeds
#mantle: driscol_bercovici14

# All in SI units
#Physical Constants.
#These are also hard coded so you cannot change their values but may be useful for defining parameters below. 
ys = 60*60*24*365     #Seconds in a year
ev = 1.602e-19        #Electron volt
kb = 1.3806485e-23    #Boltzmanns constant
G  = 6.67e-11         #Gravitational Constant
Na = 6.022140857e23   #Avogadros Constant
Rg = 8.31446261815324 #Gas constant

core         = True  #True/False. Include core in solution
stable_layer = False #True/False. Include stable layer in solution
mantle       = True  #True/False. Include mantle in solution

#core: leeds parameters

T_cmb                      = 4000  # Initial temperature of the CMB. Can specify Tcen (temperature at the center of the core) instead. 
                                   # Float.

                                    #LE concentration and properties are irrelevant for DB14 model.
conc_l                     = [0.1]  # Initial light element mass fraction of the core. Preserve order consistency with all other light 
                                    # element property lists. List(float).

core_solid_density_params  = [13088.5, 0, -8838.1/(6371e3)**2, 0]
                               # Inner core density polynomials (radial). List(float)

core_liquid_density_params = [12581.5, -1263.8/6371e3, -3642.6/(6371e3)**2, -5528.1/(6371e3)**3]
                               # Outer core density polynomials (radial). List(float).

ri                         = 1221e3  # Initial inner core radius. Float.

r_cmb                      = 3480e3  # CMB radius. Float.

core_alpha_T_params        = [1e-5]  # Core thermal expansivity pressure polynomial coefficients. List(Float).

core_cp_params             = [840]  # Core specific heat capacity pressure polynomial coefficients. List(Float).

core_conductivity_params   = ['r', 156.5027106118087, -1.1139691727354645e-06, -4.035495717259042e-12, -7.580987797074751e-20]
                               #Core thermal conductivity, radial 
                               # List, first element is a string followed by the core thermal conductivity polynomial coefficients. 
                               # The string (either 'r'/'T'/'P') indicates if the polynomials are radial, temperature, or pressure 
                               # coefficients. Can instead provide a single number which will be taken as the constant conductivity 
                               # value. List(string, Float, Float...) or List(Float).

core_melting_params        =  ['AL', 2.92154051e+03,  6.12222119e-09,  7.09821494e-21, -1.06127403e-32] #FIT TO DB14 EQ 28.
                               # List, first element is a string followed by the constants used for the melting temperature 
                               # parameterisation. The string element indicates the parameterisation to be used, if no string is 
                               # given, 'AL' (method of Alfe et al. 2002) is assumed. See melting_curve() in 
                               # core_models/leeds/routines/chemistry.py for possible options.

entropy_melting_params     =  [0] # Change of entropy on freezing pressure polynomial coefficients. List(Float).

mm                         =  [56, 16] # Molar masses (g/mol) of Iron followed by alloying light elements. List(float)

alpha_c                    =  [0] # Chemical expansivity of alloying light elements. List(float)

diffusivity_c              =  [0] # Chemical diffusivities of alloying light elements.  List(float)

use_partition_coeff        =  True # Boolean dictating if constant partition coefficients shoulf be used (True) or if partitioning based 
                                   # on chemical equilibrium should be used (False)

core_h0                    = 2e12/1.9e24  # Present day radiogenic heating per unit mass in the core

half_life                  = 1.2e9*ys  # Half life of radioactive isotope in the core.

partition_coeff            = [1]  # Partition coefficients (x_i/x_o) for mass fraction for each light element. Defined as the ratio of 
                               # inner core concentration over outer core. List(float)

lambda_sol                 = [0]  # Linear corrections to chemical potentials (solid phase). List(float)

lambda_liq                 = [0]  # Linear corrections to chemical potentials (liquid phase). List(float)

dmu                        = [0]  # Change in chemical potential between solid and liquid. List(float)

n_profiles                 = 500  # Number of nodes used in radial profiles for core properties. Float

P_cmb                      = 139e9  # CMB pressure. Float

precip_temp                = 1e10  # Temperature below which precipitation of MgO begins. Float

Cm_mgo                     = 0  # mass fraction of MgO in the core

alpha_c_mgo                = 0  # MgO chemical expansivity



#Optional parameters that if not specified take their default values

#core: leeds optional parameters

#include_baro_diffusion =   # (Default: True) Boolean, if True, barodiffusion will be included for both the entropy budget and 
                            # also chemical gradients in any chemcially stratified layer.

core_adiabat_params    =   [9.99551879e-01,  2.70599959e-09, -2.85961010e-14,  1.82335351e-21] #FIT TO DB14 EQ 27
                            # (Default: None) None or List(float). Radial polynomial coefficients for the core adiabat normalised 
                            # to T(r=0). If left as None, these will be calculated by fitting to a numerically integrated 
                            # adiabatic gradient using given alpha/cp/g profiles.

#set_cooling_rate       =   # (Default: False) Boolean, if True, the core cooling rate is set to the parameter core_dTa_dt and the 
                            # CMB heat flow is overwritten with the required value. Only used if no stable layer is included
#core_dTa_dt            =   # (Default: 0) Float, the cooling rate of the core if set_cooling_rate is True
#iron_snow              =   # (Default False) Boolean, if True, an top down core crystallisation rather than bottom up is assumed
#Ej_fixed               =   # (Default None) None or Float, if set to a value, Ej will be fixed to this value at all times. Q_cmb 
                            # will be altered to ensure this condition
#Ej_lower_bound         =   # (Default None) None or Float, if set to a value, Ej will be limited to this value if it were to fall 
                            # below it. Q_cmb will be altered to ensure this condition
#Ej_fixed_pre_ic        =   # (Default False). Boolean, when dt<0 and if true, Ej is fixed prior to ICN to the value immedietly 
                            # post ICN. Q_cmb will be altered to ensure this condition.

core_latent_heat       = 750e3 + 3e5 #Combines both latent heat and Gravitational energy in one term. Entropy will subsequently not be correct
                                     #but DB14 don't calculate entropy, just add 3e5 J/kg of energy on for gravitational release.
                            # (Default None). Latent heat is calculated using the change of entropy on melting unless this is set 
                            # to a fixed number (J/kg). None or Float

#mf_l                   =   # (Default None) Starting mole fraction of light element in the core. If specified it will overrule 
                            # conc_l.
#Tcen                   =   # (Default None) Initial temperature at the center of the core. If specified it will overrule T_cmb 
                            # (but may still be overrulled by contraint of T=Tm at ri.)
#use_new_Cr             =   # (Default False) Use new Cr factor that takes into account depression of Tm with changing LE 
                            # concentration

#mantle: driscol_bercovici14 parameters

T_surf               = 290        # Surface temperature. Float
mantle_alpha_T       = 3e-5       # Thermal expansivity. Float
g                    = 10         # Gravity (assumed constant throughout mantle. Float
mantle_diffusivity_T = 1e-6       # Thermal diffusivity. Float
mantle_cp            = 1265       # Mantle specific heat capacity. Float
Rac                  = 660        # Critical Rayleigh number. Float
activation_energy    = 3e5        # Activation energy for Arrheniun viscosity. Float
r_surf               = 6371e3     # Planet radius. Float
viscosity_ref        = None       # Reference viscosity of the mantle. If Q_surface is set, viscosity will be adjusted to give the value 
                                  # of Q_surface on the first iteration.
Tm                   = 1630/0.7   # Initial bulk mantle temperature. Float
mantle_k_upper       = 4.2        # Upper mantle thermal conductivity. Float
mantle_k_lower       = 10         # Lower mantle thermal conductivity. Float
f_visco              = 3          # Viscosity ratio between upper and lower mantle.
Qr0                  = 13e12      # Present day radiogenic heating rate. Float.
Trad                 = 2.94e9*ys  # Radiogenic half life. Float.
mantle_mass          = 4.06e24    # Mass of mantle. Float


#Following are required but are already defined above:
# r_cmb


#Optional parameters that if not specified take their default values

#mantle: driscol_bercovici14 optional parameters

Q_surface         = 39e12  # (Default: None). Initial conductive heat flow through upper mantle thermal boundary layer (inc. 
                           # crust), used to define viscosity_ref of the mantle.

#basal_magma_ocean =   # (Default: False) Toggle to include a basal magma ocean
#r_bmo             =   # (Default: None). Initial radius of the top of the basal magma ocean. Float
#bmo_liquidus      =   # (Default: None). Initial bmo liquidus temperature . Float
#FeO_mf_bmo        =   # (Default: 0) Initial mole fraction of FeO in bmo. Float
#partition_FeO     =   # (Default: 0) Partition coefficient between mantle bmo and core
#tres              =   # (Default 0) Residence time for chemical boundary layer. Float
#ds_bmo            =   # (Default 0) Change in entropy on freezing for the BMO. Float
#dc_bmo            =   # (Default 0) Difference in FeO between mantle and bmo
#rho_bmo           =   # (Default 5000) BMO density. Float

