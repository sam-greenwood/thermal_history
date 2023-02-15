import numpy as np

def tanh_conductivity(r, k_core, k_layer, r_fes, transition_width=1000):

    #If only FeS layer exists
    if r_fes <= r[0]:
        k = np.full(r.size, k_layer)
        dk_dr = np.zeros(r.size)

    #No Fes layer exists
    elif r_fes >= r[-1]:
        k = np.full(r.size, k_core)
        dk_dr = np.zeros(r.size)

    else:

        theta = np.arctanh(0.99) #value at which tanh function reaches 0.99
        x = (r - r_fes)*theta/transition_width #Normalise to distance from interface and width of transition zone

        k = k_core - ((k_core-k_layer)/2)*(1+np.tanh(x))
        dk_dx = -((k_core-k_layer)/2)*4/(np.exp(-x) + np.exp(x))**2
        dk_dr = dk_dx*theta/transition_width

    return k, dk_dr



def spacing(r, r_fes, transition_width):
    theta = np.arctanh(0.99) #value at which tanh function reaches 0.99
    return 1-np.tanh(np.abs(r - r_fes)*theta/transition_width)

def adaptive_grid(r_lower, r_fes, r_upper, coarse, fine, transition_width=50e3):

    r2 = [r_fes]

    if r_upper > r_fes:
        dr = fine
        while r2[-1] + dr < r_upper:
            r2.append(r2[-1]+dr)
            dr = coarse - (coarse-fine)*(np.max([spacing(r2[-1], r_fes, transition_width), spacing(r2[-1], r_upper, transition_width)]))

        r2[-1] = r_upper
    
    r1 = []
    if r_lower < r_fes:
        r1 = [r_fes]
        dr = fine
        while r1[-1] - dr > r_lower:
            r1.append(r1[-1]-dr)
            dr = coarse - (coarse-fine)*(np.max([spacing(r1[-1], r_fes, transition_width), spacing(r1[-1], r_lower, transition_width)]))

        r1[-1] = r_lower

    r1 = r1[::-1][:-1]

    return np.array(r1+r2)



