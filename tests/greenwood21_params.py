#Automatically generated parameters file for the following methods:
#core: leeds
#stable_layer: leeds_thermal
#mantle: thiriet19

# All in SI units
#Physical Constants.
#These are also hard coded into thermal_history.model.Parameters so you will generate an error
#message if you change their values but they may be useful for defining parameters below.
ys = 60*60*24*365     #Seconds in a year
ev = 1.602e-19        #Electron volt
kb = 1.3806485e-23    #Boltzmanns constant
G  = 6.67e-11         #Gravitational Constant
Na = 6.022140857e23   #Avogadros Constant
Rg = 8.31446261815324 #Gas constant

core         = True  #True/False. Include core in solution
stable_layer = True  #True/False. Include stable layer in solution
mantle       = True  #True/False. Include mantle in solution

#core: leeds parameters

T_cmb                      =  2327+182 # Initial temperature of the CMB. Can specify Tcen (temperature at the center of the core) instead. 
                               # Float.

conc_l                     = [0.16]  # Initial light element mass fraction of the core. Preserve order consistency with all other light 
                               # element property lists. List(float).

core_solid_density_params  = [6200]  # Inner core density polynomials (radial). List(float)

core_liquid_density_params = [6200]  # Outer core density polynomials (radial). List(float)

ri                         = 0       # Initial inner core radius. Float

r_cmb                      = 1830e3  # CMB radius. Float

core_alpha_T_params        = [1e-5]  # Core thermal expansivity pressure polynomial coefficients. List(Float)

core_cp_params             = [1089]  # Core specific heat capacity pressure polynomial coefficients. List(Float)

core_conductivity_params   = [24]  # List, first element is a string followed by the core thermal conductivity polynomial coefficients. 
                               # The string (either 'r'/'T'/'P') indicates if the polynomials are radial, temperature, or pressure 
                               # coefficients. Can instead provide a single number which will be taken as the constant conductivity 
                               # value. List(string, Float, Float...) or List(Float).

core_melting_params        = [0]  # List, first element is a string followed by the constants used for the melting temperature 
                               # parameterisation. The string element indicates the parameterisation to be used, if no string is 
                               # given, 'AL' (method of Alfe et al. 2002) is assumed. See melting_curve() in 
                               # core_models/leeds/routines/chemistry.py for possible options.

entropy_melting_params     = [kb* 1.99731, kb* -0.0082/1e9]  # Change of entropy on freezing pressure polynomial coefficients. List(Float).

mm                         = [56,32]  # Molar masses (g/mol) of Iron followed by alloying light elements. List(float)

alpha_c                    = [0.64]  # Chemical expansivity of alloying light elements. List(float)

diffusivity_c              = [0.5e-8]  # Chemical diffusivities of alloying light elements.  List(float)

use_partition_coeff        = True  # Boolean dictating if constant partition coefficients shoulf be used (True) or if partitioning based 
                               # on chemical equilibrium should be used (False)

core_h0                    = 0  # Present day radiogenic heating per unit mass in the core

half_life                  = 0  # Half life of radioactive isotope in the core.

partition_coeff            = [0]  # Partition coefficients (x_i/x_o) for mass fraction for each light element. Defined as the ratio of 
                               # inner core concentration over outer core. List(float)

lambda_sol                 = [0]  # Linear corrections to chemical potentials (solid phase). List(float)

lambda_liq                 = [0]  # Linear corrections to chemical potentials (liquid phase). List(float)

dmu                        = [0]  # Change in chemical potential between solid and liquid. List(float)

n_profiles                 = 500  # Number of nodes used in radial profiles for core properties. Float

P_cmb                      = 18.7e9  # CMB pressure. Float

precip_temp                = 0  # Temperature below which precipitation of MgO begins. Float

Cm_mgo                     = 0  # mass fraction of MgO in the core

alpha_c_mgo                = 0  # MgO chemical expansivity



#Optional parameters that if not specified take their default values

#core: leeds optional parameters

#include_baro_diffusion =   # (Default: True) Boolean, if True, barodiffusion will be included for both the entropy budget and 
                           # also chemical gradients in any chemcially stratified layer.

core_adiabat_params    = [1, -2.88429493500129e-12/1e3, -7.60393797523802e-08/1e6]  # (Default: None) None or List(float). Radial polynomial coefficients for the core adiabat normalised 
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

#core_latent_heat       =   # (Default None). Latent heat is calculated using the change of entropy on melting unless this is set 
                           # to a fixed number (J/kg). None or Float

#mf_l                   =   # (Default None) Starting mole fraction of light element in the core. If specified it will overrule 
                           # conc_l.

#Tcen                   =   # (Default None) Initial temperature at the center of the core. If specified it will overrule T_cmb 
                           # (but may still be overrulled by contraint of T=Tm at ri.)

#use_new_Cr             =   # (Default False) Use new Cr factor that takes into account depression of Tm with changing LE 
                           # concentration

#stable_layer: leeds_thermal parameters


#Following are required but are already defined above:
# core_liquid_density_params
# r_cmb
# core_alpha_T_params
# core_cp_params
# core_conductivity_params
# P_cmb


#Optional parameters that if not specified take their default values

#stable_layer: leeds_thermal optional parameters

#entrainment_T           =   # (Default: 0) Float in the range 0 <= x < 1. Entrainment coefficient, modifies the adiabatic heat 
                            # flow out the convecting region by the factor (1-x) to represent entrainment of stable fluid at the 
                            # base of the stratified layer.

#max_resolution          =   # (Default: 100) Maximum number of nodes in grid. Int

#min_resolution          =   # (Default: 10) Minumum number of nodes in grid. Int

#resolution              =   # (Default: 1/1000) Target number of nodes per meter when constructing the grid. Float.

#depth_tolerance         =   # (Default: 1000) Layer size below which the layer is assumed full eroded. Note this default value 
                            # will be halved when considering chemically stable layers fot stability. Float

#init_size               =   # (Default: 5000) Size the stable layer should be initialised if conditions promote layer growth and a 
                            # layer does not yet exist. Note this default value will be halved when considering chemically stable 
                            # layers fot stability. Float

#mix_layer               =   # (Default: False) Boolean, if True, uses an experimental function to calculate the T/c profile of a 
                            # layer that has the top portion mixed.

#thermal_stratification  =   # (Default: True). Boolean, internal flag to tell the code that thermal stratification is being used

#chemical_stratification =   # (Default: False). Boolean, internal flag to tell the code that chemical stratification is not being 
                            # used

#layer_thickness         =   # (Default: 0). Float, Initial stable layer thickness.

#mantle: thiriet19 parameters

Tm                        = 2327  # Initial mantle temperature, float

r_surf                    = 3400e3  # Surface radius, float

stagnant_lid_thickness    = 300e3  # Thickness of stagnant lid, float

mantle_k_upper            = 2  # Upper mantle thermal conducitivity, float

mantle_k_lower            = 5  # Lower mantle thermal conductivity

T_surf                    = 273  # Surface temperature, float

mantle_density            = 3550  # Mantle density, float

mantle_alpha_T            = 2.65e-5  # Mantle thermal expansivity, float

g                         = 3.7  # Surface gravity, float

mantle_cp                 = 798  # mantle specific heat capacity, float

mantle_diffusivity_T      = mantle_k_lower/(mantle_density*mantle_cp)  # Mantle thermal diffusivity, float

Rac                       = 450  # Critical Rayleigh number for mantle convction, float

Av                        = 3e5  # Activation energy for mantle viscosity, float

viscosity_ref             = 2.5e20  # Reference mantle viscosity, float

viscosity_ref_temperature = 1600  # Temperature of reference viscosity, float

a_factor                  = 2.54  # Factor controlling temperature drop across upper thermal boundary layer (Eq 15 of Thiriet et al. 
                              # 2019)

b_factor                  = 0.335 # Factor scaling boundary layer thickness with Rayleigh number.

activation_volume         = 0  #Activation volume for mantle viscosity, Float.


#Following are required but are already defined above:
# r_cmb


