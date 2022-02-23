# All in SI units
#Constants
ys = 60*60*24*365     #Seconds in a year
ev = 1.602e-19        #Electron volt
kb = 1.3806485e-23    #Boltzmanns constant
G  = 6.67e-11         #Gravitational Constant
Na = 6.022140857e23   #Avogadros Constant
Rg = 8.31446261815324 #Gas constant

core = True  #True/False. Include core in solution
stable_layer = False #True/False. Include stable layer in solution
mantle = False #True/False. Include mantle in solution

#core: simple_test parameters

rc      =  3480e3 # Radius of the core
density =  6000   # Core density
cp      =  800    # Specific heat capacity
T       =  4000   # Initial temperature of the core

#Optional parameters that if not specified take their default values

#core: simple_test optional parameters

#example_variable =   # An example optional variable for demonstration


