import importlib
import numpy as np

from copy import deepcopy
from types import ModuleType
# from .utils import get_logger
# logger = get_logger(__name__)
import logging
logger = logging.getLogger(__name__)

import thermal_history as th
from .model_classes import ThermalModel


def setup_model(parameters, core_method =         None,
                            stable_layer_method = None,
                            mantle_method =       None,
                            verbose=True,
                            log_file='out.log'):
    '''Setup main model

    Parameters
    ----------
    parameters : Parameters class
        Instance of Parameters class
    core_method : String, optional
        Name of method to solve for the core, by default None
    stable_layer_method : String, optional
        Name of method to solve for the stable layer, by default None
    mantle_method : String, optional
        Name of method to solve for the mantle, by default None
    verbose : bool, optional
        If True, progress information will be printed to STDOUT, by default True
    log_file : String, optional
        name of the log file, by default out.log

    Returns
    -------
    ThermalModel
        Instance of the main model class using the given methods and parameters

    Raises
    ------
    ValueError
        If a region is specified in the parameters but no method is given in the corresponding keyword argument
    ValueError
        If a specified method cannot be imported
    '''

    methods = {'core':         core_method,
               'stable_layer': stable_layer_method,
               'mantle':       mantle_method}

    #Check that supplied methods are consistent with regions specified in parameters
    for key, r in methods.items():
        if r == None and getattr(parameters, key):
            raise ValueError(f'{key} is set to True in parameters but a method has not been specified for it')
        elif type(r)==str and not getattr(parameters, key):
            if verbose:
                print(f'{key} is set to False in parameters but a method has been specified for that region. Ignoring method specified.')
            methods[key] = None

    required_params = {}
    optional_params = {}


    regions = [x for x in methods.keys() if getattr(parameters, x)]

    #Check each region and import the relevant python module and add to methods dict.
    for r in regions:

        required_params[r] = {}
        optional_params[r] = {}

        if hasattr(parameters, r):
            if getattr(parameters, r):
                try:
                    methods[r] = importlib.import_module(f'thermal_history.{r}_models.{methods[r]}.main')
                except Exception as e:
                    raise ValueError(f'{e}\nCannot import thermal_history.{r}_models.{methods[r]}.main')

                #Get required parameters from method 'r' and append to full required params dict
                try:
                    required_params[r].update(methods[r].required_params)
                except:
                    raise ValueError(f'No required_params defined in {methods[r]}')

                #Add optional parameters as well if they exist
                if hasattr(methods[r], 'optional_params'):
                    optional_params[r].update(methods[r].optional_params)

            else:
                raise ValueError(f'{r} is set to False in parameters')
        else:
            raise ValueError(f'{r} not set to True or False in parameters')

    assert len(regions) > 0, 'No models specified'

    #Check that all necessary inputs have been specified in parameters.
    th.utils.check_parameters(parameters, regions, required_params, optional_params, verbose=verbose)

    #Set optional parameters to their default values if not set in parameter file
    for r in regions:
        for key, value in optional_params[r].items():
            if not hasattr(parameters, key):
                setattr(parameters, key, value[1])

    return ThermalModel(parameters, methods, verbose=verbose, log_file=log_file)


class Parameters:
    '''Parameters class

    An instance of this class contains all of the parameters and is accesible to the main model
    '''    

    ys = 60*60*24*365     #Seconds in a year
    ev = 1.602e-19        #Electron volt
    kb = 1.3806485e-23    #Boltzmanns constant
    G  = 6.67e-11         #Gravitational Constant
    Na = 6.022140857e23   #Avogadros Constant
    Rg = 8.31446261815324 #Gas constant

    def __init__(self, parameters, folder='', copy=False):

        '''Initialises the class with input files

        Parameters
        ----------
        parameters : str/list/tuple
            A string or list/tuple of strings with the filenames of the parameter files. If copy is True then this can be another Parameters instance.
        folder : str, optional
            folder the parameters files exist in, to save typing out the full relative file path for each.
        copy : bool, optional
            If True, the parameters argument is instead taken to be another Parameters instance and a copy of it is returned, by default False

        Raises
        ------
        ValueError
            If a parameter is specified twice with different values
        '''        

        # if copy:
        #     for key in parameters.__dict__.keys():
        #         setattr(self, key, deepcopy(getattr(parameters, key)))

        # else:
        if not copy:
            assert type(parameters) in [str,list,tuple], 'input parameters must be a single string or list/tuple'

            if type(parameters) == str:
                parameters = [parameters]

            params = list(parameters)

            #Make sure folder has a trailing slash
            if len(folder)>1 and not folder[-1] == '/':
                folder+='/'

            for i,p in enumerate(params):
                if type(p) == str:
                    spec = importlib.util.spec_from_file_location('params'+str(i+1), folder+p)
                    params[i] = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(params[i])

        else:
            #An exisiting parameters class has been provided
            params = [parameters]


        for prm in params:
            keys = [key for key in prm.__dict__.keys() if ('__' not in key and not type(getattr(prm,key))==ModuleType)]
            for key in keys:
                value = deepcopy(getattr(prm,key))
                if type(value)==list: #Lists should be converted to arrays
                    value = np.array(value)

                if not hasattr(self,key):
                    setattr(self, key, value)

                elif not value == getattr(self,key):
                    v1, v2 = getattr(self,key), value
                    raise ValueError('multiple instances of {} in provided parameters files: \n{}\n{}'.format(key,v1,v2))
