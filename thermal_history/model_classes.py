from stat import filemode
import numpy as np
import os
import datetime
import pickle
from copy import deepcopy
import logging
logger = logging.getLogger(__name__)
import thermal_history.utils as utils

#Defined here are the classes for the 3 different regions that can be modelled, plus the main class that the user will interact with (which will contain the other 3)


class BaseModel:
    '''Base class from which each class associated with regions of the planet (core, stable layer, mantle etc)
    inherits from.

    Attributes
    ----------
    model_type : string
        String denoting the model type for internal checking
    ys : float
        Number of seconds in a year, commonly used constant.
    parameters: None
        Placeholder for parameters class, unless the region is not modelled in which
        case it is left as None.
    '''

    ys = 60*60*24*365     #seconds in a year

    model_type = 'base'  #Type of model class
    parameters = None

    def __init__(self, method_module=None):
        '''Initialises class with given method

        Parameters
        ----------
        method_module : module, optional
            Imported module for given method, by default None
        '''        

        #Set setup, evolution, and update functions
        self._method_module = method_module


    def setup(self, *args, **kwargs):
        '''Defines setup() based on the specific method in __init__().
        '''        
        return self._method_module.setup(*args,*kwargs)

    def evolve(self, *args, **kwargs):
        '''Defines evolve() based on the specific method in __init__().
        '''   
        return self._method_module.evolve(*args,*kwargs)

    def update(self, *args, **kwargs):
        '''Defines update() based on the specific method in __init__().
        '''   
        return self._method_module.update(*args,*kwargs)


    def state(self):
        ''' Creates dictionary with all attributes.
        The current 'state' of the model. This retrieves all attributes of the model class,
        with exceptions listed in the 'ignore_list' variable, and returns them in a
        dictionary. Used primarily for producing a dictionary of values to save.

        Returns
        -------
        state
            Dictionary containing model attributes (minus those ignored)
        '''

        ignore_list = ['model_type','ys','filename','parameters','profiles','_next_profiles']
        keys = self.__dict__.keys()
        keys = [key for key in keys if (not key in ignore_list and not key[0]=='_')]

        return {key: self.__dict__[key] for key in keys}

class StableLayer(BaseModel):
    '''Stable layer model class which inherits the BaseModel class.

    Attributes
    ----------
    model_type : string
        String denoting the model type for internal checking
    '''

    model_type = 'stable_layer'

    ##################

class Mantle(BaseModel):
    '''
    Mantle model class which inherits the BaseModel class.

    Attributes
    ----------
    model_type : string
        String denoting the model type for internal checking
    '''

    model_type = 'mantle'

    #DEFAULTS
    cmb_mass_flux = 0 #Mass flow

    ##################

class Core(BaseModel):

    '''
    Core model class which inherits the BaseModel class.

    Attributes
    ----------
    model_type : string
        String denoting the model type for internal checking
    '''

    model_type='core'

    ##################

class ThermalModel(BaseModel):
    '''
    Main thermal history model class. Contains the primary methods for calculating the evolution of the planet. This is the main class the user interacts with.

    Attributes
    ----------
    model_type : string
        String denoting the model type for internal checking

    time: float
        Current time of the model. Initialised at zero by default and incrimented by the time step whenever evolve() is called.
    it : int
        Iteration number (starts at 1)
    '''

    model_type = 'thermal'
    time = 0
    it = 1

    def __init__(self, parameters, methods, verbose=True, log_file='out.log', log_level=30):
        '''Setup the model based on given parameters and methods

        Parameters
        ----------
        parameters : Parameters class  
            Class containing all the parameters
        methods : Dict
            Dictionary containting the imported method modules for each region (core/stable_layer/mantle)
        verbose : bool, optional
            Print progress to STDOUT, by default True
        log_file : String, optional
            Name of the log file
        log_level : int, optional
            Numerical code for setting the logging level (10: DEBUG, 20: INFO, 30: WARNING, 40: ERROR, 50: CRITICAL), by default 30.
        '''        


        self.parameters = parameters
        if hasattr(parameters, 'time'):
            self.time = self.parameters.time

        #---------------------------------------------
        #Initialise models for regions of the planet.
        self.mantle       = Mantle(method_module=methods['mantle'])
        self.core         = Core(method_module=methods['core'])
        self.stable_layer = StableLayer(method_module=methods['stable_layer'])

        regions = ['mantle', 'core', 'stable_layer'] #Order they are calculated
        setup_order = ['core', 'stable_layer', 'mantle'] #setup must be run in radius order to calculate gravity field.
        #Add any additional regions to the planet here^
        #----------------------------------------------

        self.regions = [r for r in regions if not methods[r]==None] #list of just regions being solved.

        #Set verbose
        self.verbose=verbose
        if verbose:
            print('\nRegions to be modelled:')
            for r in self.regions:
                print('-',r)

        #Set intitial conditions
        for r in [x for x in setup_order if not methods[x]==None]:
            getattr(self,r).setup(self)


        # #Set log level and log file
        # import logging
        # if debug:
        #     logging.basicConfig(filename=self.parameters.filename+'.log', level=logging.DEBUG)
        # else:
        #     logging.basicConfig(filename=self.parameters.filename+'.log', level=logging.INFO)

        # self.logger = logging.getLogger(__name__)

        #Print out initial conditions
        if verbose:
            print('\n~~~~~ Initial Conditions ~~~~~')
            print(f'\ntime (Myr) = {self.time/(1e6*self.ys):.4e}')
            for r in self.regions:
                region = getattr(self,r)
                print('\n'+r+'\n-------')

                for key in region.state().keys():
                    value = getattr(region,key)

                    if not key in ['time','profiles','profile_names'] and not key[0]=='_':
                        if type(value) in [list,tuple,np.ndarray]:
                            print_string = ', '.join([f'{x:.4e}' for x in value])
                            print(f'{key:<10}= {print_string}')
                        elif not value == None:
                            print_string = f'{value:.4e}'
                            print(f'{key:<10}= {print_string}')


            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

        self.log_file = log_file
        handler = logging.FileHandler(filename=log_file, mode='w')
        handler.addFilter(utils.PackagePathFilter)

        logging.basicConfig(level=log_level,
                        format='%(asctime)s [%(packagepath)s:%(lineno)d] %(levelname)s: %(message)s', datefmt='%d-%b-%y %H:%M:%S',
                        handlers=[handler])

    def state(self):
        '''Create dictionary of current variables in the model.

        This retrieves all attributes of all model classes
        contained with the ThermalModel by calling their state() methods and returns them
        in a dictionary.

        Returns
        -------
        Dict
            Dictionary containing all attributes for all regions being modelled
        '''        

        full_state = {}

        for r in self.regions:
            model = getattr(self,r)
            full_state[r] = model.state()

        return full_state

    def evolve(self, dt, print_freq=100):
        '''Evolve the model one timestep.

        When this is called, the evolve method for each region is called, results are appended to the save_dict and the update method for each region is called.

        Parameters
        ----------
        dt : Float
            Time step for the model iteration (seconds)
        print_freq : int, optional
            Number of iterations between progress messages printed to screen, by default 100

        Raises
        ------
        ValueError
            If a negative step time is used when a stable layer is being solved for. Diffusion requires dt>0.
        '''        


        if self.parameters.stable_layer and dt<0:
            raise ValueError('Cannot have a negative timestep with stable layer calculation!')

        self.dt = dt


        for r in self.regions:
            if r in self.regions:
                region = getattr(self, r)
                region.time = self.time
                region.it   = self.it
                region.evolve(self)

        #Set up save_dict if not already setup
        if not hasattr(self,'save_dict'):
            self.setup_save_dict()
        else:
            self.append_to_save_dict()


        #Update values and print progress
        self.update()
        if self.verbose and self.it%print_freq == 0:
            self.print_progress()

    def update(self):      
        '''
        Calls update methods for each region and save to model classes
        '''

        #Update iteration and time
        self.it += 1
        self.time += self.dt

        #Update each calculated region
        for r in self.regions[::-1]:
            getattr(self, r).update(self)

    def setup_save_dict(self):
        '''Setup the dictionary that holds all data to be saved, stored as the save_dict attribute.

        Notes
        -----
            Only integers, floating point numbers and numpy arrays will be saved. Furthermore, the profiles
            are not saved since they take a lot of memory. If you would like to save the profiles at each time-step,
            this can be done manually after the users calls ThermalModel.evolve() at each timestep.
        '''

        save_dict = {}

        #Accepted types in save_dict
        types = (int,float,np.integer,np.floating,np.ndarray)
        #Ignore list. Only final profiles are added when write_data() method is called.
        ignore = ['profiles','next_profiles']

        #Iterate through each calculated region
        for r in self.regions:

            #Setup dictionary for specific region
            save_dict[r]={}

            model = getattr(self,r)
            state = model.state()

            keys = [key for key in state.keys() if isinstance(state[key],types) and not key in ignore]

            #Add all attributes of region that have a type given in types and are not in the ignore list.
            # for key in keys:
            #     value = state[key]
            #     if type(value) == np.ndarray:
            #         s = (len(value),)
            #     else:
            #         s = (1,)

            #     save_dict[r][key] = np.zeros(s)
            #     save_dict[r][key][:] = value

            for key in keys:
                value = state[key]
                # if type(value) == np.ndarray:
                #     s = list(value)
                # else:
                #     s = value

                save_dict[r][key] = [deepcopy(value)]


        self.save_dict = save_dict

    def append_to_save_dict(self):
        '''Append latest model data to the save_dict
        '''
        #Append latest results to the save_dict
        for r in self.regions:
            model = getattr(self,r)
            state = model.state()
            for key in state.keys():
                try:
                    # self.save_dict[r][key] = np.row_stack((self.save_dict[r][key],state[key]))

                    value = state[key]
                    # if type(value) == np.ndarray:
                    #     s = list(value)
                    # else:
                    #     s = value
                    self.save_dict[r][key].append(deepcopy(value))

                except:
                    logger.error(f'Could not append {r}.{key} to save_dict')
                    if not key in self.save_dict[r].keys():
                        logger.error(f"{key} not pre-exisiting in save_dict['{r}'].\
                            Check this variable is created on the first iteration of the model and is included when ThermalModel.setup_save_dict() is called.")
                    pass

    def print_progress(self):
        '''Print progress to STDOUT
        '''
        #Print Progress
        text = f'iteration: {self.it}'

        for r in self.regions:
            if hasattr(getattr(self, r)._method_module, 'progress'):
                text += getattr(self, r)._method_module.progress(self)

        print(text)

    def save_dict_to_numpy_array(self):
        '''Convert the lists in save_dict to numpy arrays

        Returns
        -------
        Dict
            Copy of save_dict except that all values are converted to numpy arrays.
        '''        

        x = {}
        for key1 in self.save_dict.keys():
            x[key1] = {}
            for key2 in self.save_dict[key1].keys():
                try:
                    x[key1][key2] = np.squeeze(self.save_dict[key1][key2])
                except:
                    print(f'Unable to convert save_dict[\'{key1}\'][\'{key2}\'] to a numpy array, leaving it as a list')
                    x[key1][key2] = deepcopy(self.save_dict[key1][key2])

        return x

    def write_data(self, fname):
        '''Writes the model data to a binary file.
        
        The parameters and final profiles will be added to the save_dict which is then written to
        a binary file using pickle. The time at which this function is called is also added to the
        parameters. The save_dict is also passed through save_dict_to_numpy_array() to convert all
        values to numpy arrays first.

        Parameters
        ----------
        fname : Str
            Name of the output file (.pik will be appended if not already included)

        Raises
        ------
        ValueError
            If a file with the same name already exists.
        '''

        if not '.pik' in fname and fname[-4:] == '.pik':
            fname += '.pik'

        if os.path.isfile(fname) :
            raise ValueError(f'{fname} already exists!')

        #Convert to numpy array
        numpy_dict = self.save_dict_to_numpy_array()

        #Get current time for metadata
        save_time = str(datetime.datetime.now())

        #Get all relevant attributes from parameters to add to save_dict
        keys = [key for key in dir(self.parameters) if '__' not in key]
        keys = [key for key in keys if not hasattr(getattr(self.parameters,key), '__call__')]

        #Add dictionary containing parameters
        numpy_dict['parameters'] = {'save_time': save_time}
        for key in keys:
            numpy_dict['parameters'][key] = getattr(self.parameters,key)

        #Add in final timestep profiles if they exist
        for r in ['core', 'stable_layer', 'mantle']:
            if hasattr(getattr(self, r), 'profiles'):
                numpy_dict[r]['profiles'] = getattr(self, r).profiles
     
        with open(fname, 'wb') as f:
            pickle.dump(numpy_dict, f)
            f.close()
            print('\nWrote {} to disk'.format(fname))


