"""
Basic test, check whether a simple example model can run.
"""

def test_basic():

    import thermal_history as th

    #Create parameters class with just those necessary for 'simple_test' core method.
    #This shortcut saves needing a parameters file for this test case.
    #Inherit from th.model.Parameters to get built in constants
    class example_parameters(th.model.Parameters):

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

    prm = example_parameters()

    #Setup the model
    model = th.model.setup_model(prm, core_method='simple_test', verbose=False)

    #Iterate the model 1s
    model.mantle.Q_cmb = 1
    model.evolve(1)

    #Save model output
    model.write_data('example')

if __name__ == '__main__':
    test_basic()