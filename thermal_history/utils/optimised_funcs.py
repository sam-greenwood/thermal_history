import numpy as np
from numba import njit
 
@njit
def linspace(x1,x2,n):
    '''
    Jit compiled version of numpy.linspace().
    '''
    x = np.zeros(n)
    x[0] = x1
    x[-1] = x2
    for i in range(1,n-1):
        x[i] = x[0] + i*(x2-x1)/(n-1)
    return x

@njit
def trapezoid(x,y):
    '''
    Jit compiled version of scipy.integrate.cumulative_trapezoid().
    Assumes initial=0.
    '''

    integrand = np.zeros(len(x))

    for i in range(1, len(x)):
        integrand[i] = integrand[i-1] + (x[i]-x[i-1])*0.5*(y[i]+y[i-1])

    return integrand

@njit
def jit_polyval(poly, x):
    '''
    Jit compiled version of np.polyval().
    '''

    n = len(poly)
    y = np.ones(len(x))*poly[-1]

    exponent = n-1
    for i in range(n-1):
        y = poly[i] *x**exponent + y
        exponent -= 1

    return y

def polyval(poly,x):
    if not hasattr(x, '__iter__'):
        y = float(jit_polyval(poly,np.array([x])))
    else:
        y = jit_polyval(poly,x)

    return y