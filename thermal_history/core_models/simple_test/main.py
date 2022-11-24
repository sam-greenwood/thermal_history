#Description of this model3
Description = 'Simple example core model.'

#import our logger
import logging
logger = logging.getLogger(__name__)

import numpy as np

#3 required functions for any method: setup, evolve, and update

def setup(model):
    #Sets up any inital conditions

    core = model.core
    prm = model.parameters

    core.T = prm.T #Set initial core temperature

def evolve(model):
    #Main model

    core = model.core
    prm = model.parameters

    Q = model.mantle.Q_cmb

    mass = (4/3)*np.pi*prm.rc**3 * prm.density

    dT_dt = Q/(-mass*prm.cp)

    #Variables to be saved need to attributes of model.core
    core.dT_dt = dT_dt


def update(model):
    #Apply any changes between timesteps. Attributes of model.core will be saved to 

    core = model.core

    #Update core temperature
    core.T += core.dT_dt*model.dt

    #Example use of the logger
    if core.T < 0:
        logger.critical('Core temperature is negative!')

def progress(model):
    #Text to be printed to STDOUT on each iteration

    core = model.core

    text = f'    temperature = {core.T:.2f} K'

    return text

#list required parameters {'name': 'description'}
required_params = {'rc': 'Radius of the core',
                    'density': 'Core density.',
                    'cp': 'Specific heat capacity',
                    'T': 'Initial temperature of the core'}

#list optional parameters that if not specified in model.parameters, will get set to a default value.
#{'name': ('description', default_value)}
optional_params = {'example_variable': ('An example optional variable for demonstration purposes', True)}