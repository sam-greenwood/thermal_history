import numpy as np
from thermal_history.utils.optimised_funcs import linspace, trapezoid, polyval


def radial_density_grid(rc, rb, rm, rl, r_crust, r_reg, r_surf, rho_man, rho_crust):

    n = 10 #Number of grid points in triple conductive layers

    r = np.zeros(4 + 3*n)

    r[:4] = [rc, rb, rm, rl]

    r[4:4+n]     = linspace(rl, r_crust, n)
    r[4+n:4+2*n] = linspace(r_crust, r_reg, n)
    r[4+2*n:]    = linspace(r_reg, r_surf, n)

    rho = np.ones(r.size)*rho_man
    rho[-2*n:] = rho_crust

    return r, rho

def pressure(r, rho, g):
    P = trapezoid(r[::-1], -rho[::-1]*g)[::-1]
    return P

def liquidus(pressure, liquidus_params):

    T = polyval(np.array(liquidus_params[::-1]), pressure)

    return T

def solidus(pressure, solidus_params):

    T = polyval(np.array(solidus_params[::-1]), pressure)

    return T

