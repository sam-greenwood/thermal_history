import importlib
import os

# #Set numba logger to warning only
# numba_logger = logging.getLogger('numba')
# numba_logger.setLevel(logging.WARNING)

#Custom filter to get logger only paths within the package thermal_history rather than absolute paths
import logging
class PackagePathFilter(logging.Filter):
    '''Custom filter to extract just the pathname within the thermal_history package
    '''
    def filter(record):
        pathname = record.pathname

        split = pathname.split('/')

        for i in range(len(split)-1, -1, -1):
            if split[i] == 'thermal_history':
                break
        record.packagepath = '/'.join(split[i:])

        return True

        

def check_parameters(parameters, method_names, required_params):

    error = False

    message = '------- Parameter Check -------'
    n1 = len(message)
    print('\n'+message)

    result = ' Passed '
    All_params=[] #Keep track of params given but not required
    for f in method_names:
        All_params+=[x for x in required_params[f].keys() if not x in All_params]

    #Check required parameters have been given
    for f in method_names:
        for key, description in required_params[f].items():

            if not hasattr(parameters, key):
                print('{} method requires {} parameter. \n{} : {}\n'.format(f, key, key, description))
                error = True
                result = ' Failed '

            All_params=[x for x in All_params if not x==key]

    if len(All_params)>0:
        print('Following parameters have been specified but are not specifically listed as required by the given methods')
        for key in All_params:
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

    Physical_constants = ['#Constants',
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
    max_l = max([len(x) for x in regions])
    for key, value in required_params.items():
        if key in regions:
            f.write(f'{key:<{max_l}}'+ f' = True  #{value}\n')
        else:
            f.write(f'{key:<{max_l}}'+ f' = False #{value}\n')
    f.write('\n')

    required_params = {}
    method_modules = {}

    #Check each region and import the relevant python module and add to methods dict.
    for i in range(len(regions)):
        r = regions[i]
        required_params[r] = {}
        method_modules[r] = {}

        #Header for method
        f.write(f'#{r}: {methods[r]} parameters\n\n')

        #Try to import required_params from methods.
        try:
            method_modules[r] = importlib.import_module(f'thermal_history.{r}_models.{methods[r]}.main')
        except Exception as e:
            raise ValueError(f'\n\nCannot import thermal_history.{r}_models.{methods[r]}.main\n\n{e}')

        required_params[r].update(method_modules[r].required_params)

        #Track maximum lenght of parameter name to neatly format
        max_l = max([len(x) for x in required_params[r].keys()])

        #Write out parameter and track any repeated parameters, sometimes multiple methods require the same one.
        repeated=[]
        for key, value in required_params[r].items():
            repeat=False
            if i>0:
                for j in range(i):
                    if key in required_params[regions[j]].keys():
                        repeat=True
            if repeat:
                repeated.append(key)
            else:
                line = f'{key:<{max_l}}'+f' =   # {value}\n'
                if len(line) > 72: #Try to cut down long lines across 2 lines.
                    end = 72
                    for c in range(72,50,-1):
                        if line[c]==' ':
                            end = c
                            break

                    space = ' '
                    f.write(line[:end]+'\n')
                    f.write(f'{space:<{max_l}}'+f'     # {line[end:]}')
                else:
                    f.write(line)

        if len(repeated)>0:
            f.write('\n#Following are required but are already defined above:\n')
            for x in repeated:
                f.write(f'# {x}\n')

        f.write('\n\n')

    #Write any optional parameters to file but all commented out.
    optional_params = {}
    f.write('#Optional parameters that if not specified take their default values\n\n')

    for r in regions:
        try:
            optional_params[r] = {}
            optional_params[r].update(method_modules[r].optional_params) #Get optional_params if there are any
        except:
            pass


        #Header for method
        f.write(f'#{r}: {methods[r]} optional parameters\n\n')

        max_l = max([len(x) for x in optional_params[r].keys()])
        for key, value in optional_params[r].items():
            line = f'#{key:<{max_l}} =   # {value[0]}\n'
            if len(line) > 72: #Try to cut down long lines across 2 lines.
                end = 72
                for c in range(72,40,-1):
                    if line[c]==' ':
                        end = c
                        break
                space = ' '
                f.write(line[:end]+'\n')
                f.write(f'{space:<{max_l}}'+f'      #{line[end:]}')
            else:
                f.write(line)

        f.write('\n\n')

    f.close()
