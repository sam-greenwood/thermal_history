import numpy as np
from scipy.integrate import cumulative_trapezoid

def radial_density_grid(radii, density_params):
    '''
    Calculates density profile across n shells, each with their own specified
    density polynomials.

    radii: list of interface radii of shells (size n+1), starting with CMB radius and ending
    with radius at the surface.

    density params: list of arrays, where each array gives the polynomial coefficients for the
    density.
    '''

    #100 radial points for each layer

    r   = np.array([])
    rho = np.array([])
    for i in range(len(radii) -1):
        r_temp = np.linspace(radii[i], radii[i+1], 100)
        
        rho = np.append((rho,np.polyval(density_params[::-1], r_temp)))
        r   = np.append((r, r_temp))

    return r, rho

def gravity(r, rho, g_cmb):
    '''
    r = radius (singluar value or array)
    rho = best-fit density polynomials (radial)
    g_cmb = CMB gravity

    returns: gravity at r
    '''
    G = 6.67e-11
    g = (4*np.pi*G)*cumulative_trapezoid(rho*r**2, x=r)/r**2 + g_cmb/(r**2-r[0]**2)

    return g

def pressure(r_surf, r_crust, r_cmb, rho_crust, rho_mantle, g):
    '''
    Pressure in the mantle. Assumed constant gravity throughout and constant
    density within the crust and mantle. Pressure is therefore linear between
    the surface and crustal base, then the crustal base and CMB. Just the pressures
    at r_crust, and r_cmb are returned.
    '''

    P0 = 0                                      #Surface pressure
    P1 = P0 + rho_crust * g * (r_surf-r_crust)  #Pressure at base of crust
    P2 = P1 + rho_mantle* g * (r_crust-r_cmb)   #Pressure at CMB

    return P1, P2


def liquidus(pressure, liquidus_params):

    T = np.polyval(liquidus_params[::-1], pressure)

    return T

def solidus(pressure, solidus_params):

    T = np.polyval(solidus_params[::-1], pressure)

    return T

def setup_profiles(model):

    mantle = model.mantle
    prm    = model.parameters

    mantle.profiles = {}
    mantle.profiles['r'], mantle.profiles['rho'] = radial_density_grid([prm.r_cmb, prm.r_crust, prm.r_surf],
                                                                       [])
