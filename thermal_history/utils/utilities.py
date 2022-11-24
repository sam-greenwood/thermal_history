import importlib
import os
from types import ModuleType

# #Set numba logger to warning only
# numba_logger = logging.getLogger('numba')
# numba_logger.setLevel(logging.WARNING)

#Custom filter to get logger only paths within the package thermal_history rather than absolute paths
import logging
logger = logging.getLogger(__name__)

class PackagePathFilter(logging.Filter):
    '''Custom filter to extract just the pathname realtive to the thermal_history package top level.
    '''
    def filter(record):
        pathname = record.pathname

        split = pathname.split('/')

        for i in range(len(split)-1, -1, -1):
            if split[i] == 'thermal_history':
                break
        record.packagepath = '/'.join(split[i:])

        return True



def check_parameters(parameters, method_names, required_params, optional_params, verbose=True):

    error = False

    message = '------- Parameter Check -------'
    n1 = len(message)
    if verbose:
        print('\n'+message)

    result = ' Passed '

    #Check required parameters have been given
    for f in method_names:
        for key, description in required_params[f].items():

            if not hasattr(parameters, key):
                print('{} method requires {} parameter. \n{} : {}\n'.format(f, key, key, description))
                error = True
                result = ' Failed '

    #Parameters given but not explicitly required

    All_params=[] #All required/optional parameters
    for f in method_names:
        All_params+=[x for x in required_params[f].keys() if not x in All_params]
        All_params+=[x for x in optional_params[f].keys() if not x in All_params]

    given_params = [key for key in parameters.__dict__.keys() if ('__' not in key and not type(getattr(parameters,key))==ModuleType)] #Params specified
    not_required = [key for key in given_params if not key in All_params] #Exclude required/optional parameters
    not_required = [key for key in not_required if not key in ['core', 'mantle', 'stable_layer']] #Exclude 'core'/'mantle'/'stable_layer'. These are assumed required

    if len(not_required) > 0:
        if verbose:
            print('Following parameters have been specified but are not specifically listed as required by the given methods:')
            for key in not_required:
                print(f'{key} = {getattr(parameters, key)}')


    #Check regions are defined
    regions = {'core': 'True/False. Toggle core in solution',
               'stable_layer': 'True/False. Toggle stable layer in solution',
               'mantle': 'True/False. Toggle mantle in solution'}

    for key, description in regions.items():
        if hasattr(parameters, key):
            if getattr(parameters, key) and not key in method_names:
                print(f'{key} is True in parameters but has no method specified')
                error = True
            elif not getattr(parameters, key) and key in method_names:
                print(f'{key} is False in parameters but a method for {key} has been specified')

        else:
            print(f'{key} must be defined\n{key}:{description}')
            error = True

    #Check physical constants are defined:
    constants = {'ys': 'Seconds in a year',
                 'ev': 'Electron volt (J)',
                 'kb': 'Boltzmann constant (J/K)',
                 'G': 'Gravitational Constant',
                 'Na': 'Avogadros Constant',
                 'Rg': 'Gas constant (J/K/mol)'}

    for key, description in constants.items():

        if not hasattr(parameters, key):
            print(f'Physical constant {key}:{description} must be defined')
            error = True
            result = ' Failed '


    n2 = int((n1-len(result))/2)
    if verbose:
        print('-'*n2 + result + '-'*n2 + '\n')

    if error == True:
        raise ValueError('\n\nParameters missing/inconsistent with declared methods (see above parameter check message)')


def create_parameters_file(fname,
                           core_method =         None,
                           stable_layer_method = None,
                           mantle_method =       None):

    '''
    Generate an input parameters file for the specified methods.

    args
    ------
    fname: file name for generated parameters file

    kwargs
    ------
    core_method: (None) String of core method name
    stable_layer_method: (None) String of stable layer method name
    mantle_method: (None) String of mantle method name
    '''


    if not '.py' in fname:
        fname += '.py'

    assert not os.path.isfile(fname), f'{fname} already exists!'

    methods = {'core': core_method,
               'stable_layer': stable_layer_method,
               'mantle': mantle_method}

    required_params = {'core': 'True/False. Include core in solution',
                       'stable_layer': 'True/False. Include stable layer in solution',
                       'mantle': 'True/False. Include mantle in solution'}

    Physical_constants = ['#Physical Constants.',
                          '#These are also hard coded so you cannot change their values but may be useful for defining parameters below. ',
                          'ys = 60*60*24*365     #Seconds in a year',
                          'ev = 1.602e-19        #Electron volt',
                          'kb = 1.3806485e-23    #Boltzmanns constant',
                          'G  = 6.67e-11         #Gravitational Constant',
                          'Na = 6.022140857e23   #Avogadros Constant',
                          'Rg = 8.31446261815324 #Gas constant']


    regions = [x for x in methods.keys() if not methods[x] == None]
    assert len(regions) > 0, 'No models specified'

    #Open file and write header text and physical constants
    f = open(fname, 'w')
    f.write('#Automatically generated parameters file for the following methods:\n')
    for r in regions:
        f.write(f'#{r}: {methods[r]}\n')
    f.write('\n')
    f.write('# All in SI units\n')

    #Write out physical constants
    for x in Physical_constants:
        f.write(x+'\n')
    f.write('\n')

    #Write out toggles to control which regions are included.
    max_l = max([len(x) for x in methods.keys()])
    for key, value in required_params.items():
        if key in regions:
            f.write(f'{key:<{max_l}}'+ f' = True  #{value}\n')
        else:
            f.write(f'{key:<{max_l}}'+ f' = False #{value}\n')
    f.write('\n')

    method_modules = {}
    required_params = {}
    optional_params = {}

    #Check each region and import the relevant python module and add to methods_modules dictionary.
    for i in range(len(regions)):
        r = regions[i]
        required_params[r] = {}
        optional_params[r] = {}
        method_modules[r] = {}

        #Header for method
        f.write(f'#{r}: {methods[r]} parameters\n\n')

        #Try to import main.py from given methods and get the required/optional parameters
        try:
            method_modules[r] = importlib.import_module(f'thermal_history.{r}_models.{methods[r]}.main')
        except Exception as e:
            raise ValueError(f'\n\nCannot import thermal_history.{r}_models.{methods[r]}.main\n\n{e}')

        try:
            required_params[r].update(method_modules[r].required_params)
        except:
            logger.warning(f'No \'required_params\' dictionary from main.py of method {methods[r]} exists.\nSkipping...')
            pass
        try:
            optional_params[r].update(method_modules[r].optional_params)
        except:
            logger.warning(f'No \'optional_params\' dictionary from main.py of method {methods[r]} exists.\nSkipping...')
            pass

    
        #Write required parameters

        if len(required_params[r].keys())>0:

            #Track maximum length of parameter names to neatly format
            max_l = max([len(x) for x in required_params[r].keys()])

            #Write out parameters
            repeated=[]
            for key, value in required_params[r].items():

                #Track any repeated parameters, sometimes multiple methods require the same one.
                repeat=False
                if i>0:
                    for j in range(i):
                        if key in required_params[regions[j]].keys():
                            repeat=True

                if repeat:
                    repeated.append(key)
                else:
                    #write parameter
                    write_line(f, key, value, max_l)
                    f.write('\n')

            #Write out any repeated parameters commented out.
            if len(repeated)>0:
                f.write('\n#Following are required but are already defined above:\n')
                for x in repeated:
                    f.write(f'# {x}\n')

        else:

            f.write(f'No required parameters are specified in main.py for {methods[r]}')

        f.write('\n\n')

        if len(optional_params[r].keys()) > 0:

            #Write any optional parameters to file but all commented out.
            f.write('#Optional parameters that if not specified take their default values\n\n')

            #Header for method
            f.write(f'#{r}: {methods[r]} optional parameters\n\n')

            #Max length of strings for neat formatting
            max_l = max([len(x) for x in optional_params[r].keys()])
 
            for key, value in optional_params[r].items():

                #Write parameter
                write_line(f, key, value[0], max_l, comment_out=True)
                f.write('\n')

    f.close()


def write_line(f, key, description, max_l, comment_out=False):
    '''Write parameter to file in a formatted manner.
    
    Used by create_parameters_file().

    Parameters
    ----------
    f : _io.TextIOWrapper
        File handle for parameters 
    key : string
        name of parameter
    description : string
        description of parameter
    max_l : int
        Number of characters before '=' sign for neat formatting.
    comment_out : bool, optional
        prefix line with a comment ('#'), by default False
    '''

    n_char = 100

    variable_text = f'{key:<{max_l}} =   '
    if comment_out:
        variable_text='#'+variable_text

    f.write(variable_text)

    remaining_text = description

    start = True
    while len(remaining_text) > 0:

        if len(remaining_text) > n_char:
            end = min([n_char, len(remaining_text)])
            for c in range(100,0,-1):
                if remaining_text[c]==' ':
                    end = c+1
                    break
        else:
            end = len(remaining_text)

        text = '# '+remaining_text[:end]
        if start:
            start = False
        else:
            space = ' '
            text = f'{space:<{max_l+5}}'+text


        f.write(text+'\n')
        remaining_text = remaining_text[end:]


def get_available_models():


    #Top level filepath
    path = __file__.split('utils/utilities.py')[0]

    print('\nAvailable models:')

    errors = []
    for region in ['core', 'stable_layer', 'mantle']:

        print('\n*******************')
        print(f'{region}_models')
        print('*******************')
        # Get list of available directories present
        methods = [x for x in os.listdir(path+region+'_models') if os.path.isdir(path+region+'_models/'+x) and not '__pycache__' in x]

        for method in methods:

            #Try and import main.py for each method
            try:
                _ = importlib.import_module(f'thermal_history.{region}_models.{method}.main')
                if hasattr(_, 'Description'):
                    print(f'{method: <20} - {_.Description}')
                else:
                    print(f'{method: <20} - (No description available)')

            except Exception as e:
                errors.append([region, method, e])
    
    #Print any error messages from failed imports
    if len(errors)>0:
        print('\n\nErrors:')
        for error in errors:
            print('******************')
            print(f'thermal_history.{error[0]}_models.{error[1]}.main cannot be imported. Error message:\n{error[2]}')
            print('******************\n')