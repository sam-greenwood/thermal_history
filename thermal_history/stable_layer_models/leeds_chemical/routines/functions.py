import numpy as np


from scipy.interpolate import interp1d
from scipy.special import erfc
from scipy.optimize import bisect
from scipy.special import erfcinv

import numba
from numba import jit, njit

from ....core_models.leeds.routines import profiles as prof

import pdb


def primordial_layer_init(conc_l, dc, dc_dr_s, r):
    '''
    Initial compositional profile within a primordial layer

    Parameters
    ----------
    conc_l
        Mass fraction of light element in isentropic region
    dc
        Change in mass fraction of light element across layer
    dc_dr_s
        gradient in mass fraction at the base of the layer to satisfy
    r
        radial grid within layer to calculate mass fraction on

    Returns
    -------
    c
        Array with mass fraction values on the radial grid
    '''

    #Make compositional profile linear but ensure first 3 grid points satisfy stable chemical gradient above
    #This ensures numerical stability
    c = np.zeros(r.size)
    c[:3] = conc_l + dc_dr_s*(r[:3]-r[0])
    c[3:] = np.linspace(c[2], conc_l+dc, r.size-2)[1:]

    return c

def determine_strat_type(ub_type_c, ub_c, ADR, T, Ta, conc_l, E_T, prm):
    '''
    AI algorithm (lots of 'if' statements) to determine the dominant type of stratification: thermal or chemical

    Parameters
    ----------
    ub_type_T
        Type of upper boundary condition type for the thermal solution. Either 0 (fixed value) or 1 (fixed gradient)
    ub_T
        Value for the upper boundary condition for temperature
    ub_type_c
        Type of upper boundary condition type for the chemical solution. Either 0 (fixed value) or 1 (fixed gradient)
    ub_c
        Value for the upper boundary condition for composition
    ADR
        Adiabatic ratio, the ratio of CMB heat flux to the adiabatic heat flux
    T
        Temperature values on radial grid
    Ta
        Adiabatic temperature values on radial grid
    conc_l
        Mass fraction of light element in isentropic region
    E_T
        Entrainent parameters for thermal solution
    prm
        parameter class

    Returns
    -------
    strat
        String containing the type of dominant stratification ('thermal', 'chemical' or 'none')
    '''

    c_strat = prm.compositional_stratification
    T_strat = prm.thermal_stratification

    #Decide on dominant driver of stratification (thermal or chemical)
    if c_strat:

        #No thermal stratification
        if not T_strat:
            strat = 'chemical'

        #If super-adiabatic
        elif ADR >= 1:
            strat = 'chemical'

        #if thermal profile is still sub-adiabatic at interface
        elif T_strat and (T[0]-T[1]) < (Ta[0]-Ta[1]):
            strat = 'thermal'

        else:
            strat = 'chemical'

        #Make sure BC's on composition promote growth of a layer
        if strat == 'chemical':

            if (ub_type_c == 0 and ub_c <= conc_l) or (ub_type_c == 1 and ub_c <= 0):
                strat = 'none'


    elif T_strat:

        #Sub-adiabatic
        if ADR < 1-E_T:
            strat = 'thermal'

        #otherwise no stratification
        else:
            strat = 'none'

    else:
        strat = 'none'
        #assert T_strat or c_strat, 'No stratification enabled in parameters'

    if prm.primordial_layer:
        strat = 'chemical'

        if (ub_type_c == 0 and ub_c <= conc_l) or (ub_type_c == 1 and ub_c <= 0):
            strat = 'none'

    return strat


#Stable layer functions
def buffett_seagle_10_growth(r,c,conc_l,dc_dr):

    '''
    Calculate the layer interface movement using the method of Buffett and Seagle (2010)
    and regrid the solution.

    Parameters
    ----------
    r
        radial grid
    c
        Light element mass fraction on grid
    conc_l
        Mass fraction of isentropic region
    dc_dr
        Mass fraction gradient at interface in the layer

    Returns
    -------

    r,T,c
        Regridded solutions (radius, temperature, mass fraction)
    '''

    s = r[0]

    ds = (conc_l-c[0])/dc_dr

    s_new = float(s + ds)

    # T_rel = T - prof.adiabat(r,Tcen,adiabat_poly)
    # c_rel = c - conc_l


    # r, (T_rel,c_rel) = change_domain_size((T_rel,c_rel), r, s_new, resolution)
    # c = c_rel + conc_l
    # T = T_rel + prof.adiabat(r,Tcen,adiabat_poly)

    #T = adiabat(r,Tcen)


    return s_new



def mix_profile(x, I_prime, y_in):

    dx = x[1]-x[0]

    y = y_in.copy() #Make a copy as it is modified in place

    I = I_prime*y

    start=0
    for l in range(y.size-1):

        if y[l+1] < y[l]:

            start, end = l, l+1
            flag = True
            while flag and end < y.size:

                integral_before =  dx/2*(I[start] + 2*np.sum(I[start+1:end]) + I[end])

                y_cst = (integral_before*(2/dx))/(I_prime[start] + 2*np.sum(I_prime[start+1:end]) + I_prime[end])

                if start>0 and y_cst < y[start-1]:
                    start += -1

                elif end < y.size-1 and y_cst > y[end+1]:
                    end += 1

                else:
                    flag=False


            y[start:end+1] = y_cst

            y = mix_profile(x, I_prime, y)

    return y



def retreat_layer(r, T, c, Tcen, conc_l, adiabat_poly, resolution, alpha_T, alpha_c, density_grad_limit=0):
    '''
    Regrid the solution back to the radius that satisfies the stability conditions

    Parameters
    ----------
    r
        radial grid
    T
        Temperature on grid
    c
        Light element mass fraction on grid
    Tcen
        Temperature at r=0
    conc_l
        Mass fraction of isentropic region
    adiabat_poly
        Radial polynomials for adiabatic temperature
    resolution
        tuple of resolution parameters (number of grid points/meter, minimum number of points, maximum number of points)
    alpha_T
        Thermal expansivity
    alpha_c
        Chemical expansivity

    density_grad_limit
        Upper limit on the potential density gradient required for stability

    Returns
    -------

    r_new,T_new,c_new
        Regridded solutions (radius, temperature, mass fraction)
    '''

    d_rho_T = -alpha_T*(T - prof.adiabat(r,Tcen,adiabat_poly))   #Change in potential density due to temperature
    d_rho_c = -alpha_c*(c - conc_l)   #Change in potential density due to composition

    d_rho = d_rho_T + d_rho_c     #Total change in potential density

    if np.min(d_rho) >= 0:

        s_new = r[-1]

    elif np.max(d_rho) <= 0:

        s_new = r[0]

    else:

        for i in range(r.size):

            if d_rho[i] > 0:
                s_new = r[i]

            elif i < r.size-1 and (d_rho[i+1]-d_rho[i])/(r[i+1]-r[i]) > density_grad_limit:

                s_new = r[i]


    return s_new


def change_domain_size(X_rel,r,s_new,resolution):
    '''
    Change the domain size with linear interpolation

    Parameters
    ----------
    X_rel
        Property (array) or properties (tuple of arrays) to be regridded
    r
        radial grid of property/properties
    s_new
        new radius of the layer interface
    resolution
        tuple of resolution parameters (number of grid points/meter, minimum number of points, maximum number of points)

    Returns
    -------
    r_new, a
        New radial grid and tuple of arrays with regridded properties
    '''

    s_new = float(s_new)

    resolution, min_res, max_res = resolution

    if not type(X_rel) == tuple:
        X_rel = tuple([X_rel])

    r_cmb = r[-1]
    s = r[0]

    n_points = int(np.max([min_res,resolution*(r[-1]-s_new)]))
    n_points = int(np.min([max_res,n_points]))

    a = ()
    for x in X_rel:

        if s_new == r_cmb:
            r_new = np.ones(n_points)*r_cmb
            x_new = np.zeros(n_points)

        else:
            #Shrink layer
            if s_new > s:

                r_new = np.linspace(s_new,r_cmb,n_points)
                x_new = np.interp(r_new,r,x)

                if np.min(np.gradient(x_new,np.diff(r_new)[0]))>0:
                    x_new[x_new<0]=0

                x_new[0] = 0



            #Expand layer
            elif s_new < s:

                r_append = np.append(s_new,r)
                x_append = np.append(0,x)

                r_new = np.linspace(s_new,r_cmb,n_points)
                x_new = np.interp(r_new,r_append,x_append)

            #Layer has not moved
            else:

                r_new = r
                x_new = x


        a = a + tuple([x_new])

    if len(a) == 1:
        a = a[0]

    return r_new, a


    ###############################################################################
    ###############################################################################

#################################################################################
#################################################################################

def cubic_fit(x1,x2,f_x1,f_x2,f_x1_prime,f_x2_prime):
    '''
    Fit a cubic equation f(x) = m1*x^3 + m2*x^2 + m3*x + m4
    Defined conditions are f(x1), f(x2) and f'(x1), f'(x2)
    '''

    A = np.zeros([4,4])

    A[0,:] = [x1**3, x1**2, x1, 1]
    A[1,:] = [x2**3, x2**2, x2, 1]
    A[2,:] = [3*x1**2, 2*x1, 1, 0]
    A[3,:] = [3*x2**2, 2*x2, 1, 0]

    B = np.array([f_x1,f_x2,f_x1_prime,f_x2_prime])
    A = np.array(A)

    m = np.linalg.solve(A,B)

    return m





def pure_thermal_check(self):

    sl = self.stable_layer

    T_grad_cmb = eval_dT_dr_cmb(self)

    ADR = T_grad_cmb/self.core.profiles['dTa_dr'][-1]
    sl.ADR = ADR

    E_T = self.parameters.entrainment_T

    #If super adiabatic, layer is not present
    if ADR >= (1-E_T):
        sl.T_grad_s, sl.c_grad_s = 0, 0
        return False

    return True


def pure_chemical_check(self):

    prm = self.parameters

    ub_type, ub = self.mantle.chemical_bc_cmb
    conc_l = self.core.conc_l[0]

    if (ub_type == 1 and ub <= 0) or (ub_type == 0 and ub <= conc_l):
        self.stable_layer.T_grad_s, self.stable_layer.c_grad_s = 0, 0
        return False

    return True


def thermo_chemical_check(self):

    if pure_chemical_check(self) or pure_thermal_check(self):
        return True
    else:
        self.stable_layer.T_grad_s, self.stable_layer.c_grad_s = 0, 0
        return False




def eval_dT_dr_cmb(self):

    Q_cmb = self.mantle.Q_cmb
    k_cmb = self.core.profiles['k'][-1]
    r_cmb = self.parameters.r_cmb

    T_grad_cmb = Q_cmb / (4*np.pi*r_cmb**2 * -k_cmb)

    return T_grad_cmb


def setup_thermal_profile(self):

    prm = self.parameters

    r = self.stable_layer.profiles['r']
    Tcen = self.core.Tcen

    self.stable_layer.profiles['T'] = prof.adiabat(r, Tcen, prm.core_adiabat)
    self.stable_layer.T_grad_s = prof.adiabat_grad(r[0], Tcen, prm.core_adiabat)


def setup_chemical_profile(self):

    prm = self.parameters

    r = self.stable_layer.profiles['r']
    conc_l = self.core.conc_l[0]

    self.stable_layer.profiles['c'] = np.ones(r.size)*conc_l
    self.stable_layer.c_grad_s = 0

    if self.it == 1 and prm.primordial_layer:

        dc = prm.primordial_layer_excess_c

        #Calculate the stable chemical gradient at the base of the layer assuming CMB Temp gradient exists.
        #Ensures stability if initial profile can satisfy this in the first couple of grid points.
        T_grad_cmb = eval_dT_dr_cmb(self)
        dc_dr_s = -(prm.core_alpha_T/prm.alpha_c[0])*(T_grad_cmb-prof.adiabat_grad(r[0],self.core.Tcen,prm.core_adiabat))

        self.stable_layer.profiles['c'] = func.primordial_layer_init(conc_l, dc, dc_dr_s, r)
        self.stable_layer.c_grad_s = dc_dr_s

        # Ts = prof.adiabat(r[0], self.core.Tcen, prm.core_adiabat)

        # self.stable_layer.profiles['T'] = np.linspace(Ts, Ts+T_grad_cmb*self.stable_layer.layer_thickness, r.size)
        # self.stable_layer.T_grad_s = T_grad_cmb*1

    # self.stable_layer.c_grad_s = 0
