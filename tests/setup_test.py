import thermal_history as th

#Create parameters class with just those necessary for 'simple_test' core method.
#Inherit from th.model.Parameters to get built in constants
class test_parameters(th.model.Parameters):

    #Redefine init function as we are not specifying any input file, just including
    #the params needed for this test.
    def __init__(self):
        self.core         = True
        self.mantle       = False
        self.stable_layer = False

        self.rc      = 1  #Radius of the core
        self.density = 1  #Core density
        self.cp      = 1  #Specific heat capacity
        self.T       = 1  #Initial temperature of the core

prm = test_parameters()

model = th.model.setup_model(prm, core_method='simple_test', verbose=False)