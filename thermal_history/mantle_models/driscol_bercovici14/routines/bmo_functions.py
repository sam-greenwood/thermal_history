
import numpy as np
import scipy


beta  = 1.0/3.0
k     = 8
delta = 100e3
nu    = 1e10
cp    = 1000
ds    = 300.0
cp_core = 860
rhom    = 5500
M_core  = 2e24
gamma_a = 5.0e-7
gamma_c = 1
gamma_m = 0.9e-7
eta     = 0.7
AO, AMg, AFe = 16, 24, 56


def latent_coeff(Tl,ds):
    L = Tl*ds
    return L
        
def latent_labrosse(r_l,r_t,L,dc):
    """Ql in eqn 1 of Labrosse is -4*pi*a^2*rho*L*da/dt and 
       da/dt = -(a^3 - b^3) /(a^2*dT*dc) * dTl/dt where dT < 0
       so the -ve sign goes away. 
    """
    dT = 3500-5500
    num = 4 * np.pi * rhom * L * (r_t**3 - r_l**3)
    den = 3 * dc * dT
    return num/den
    
def cmb_flux_mass(cc, cb, delta, rho, D):
    dcdr = -(cc-cb)/delta
    flux = -rho * D * dcdr
    return flux
   
def dcdt_fac_labrosse(r_l,r_t,dc):
    """eqn 2 of Labrosse    """
    num = 3 * r_t**2 * dc
    den = r_t**3 - r_l**3
    return -num/den
    
def dadt_fac_labrosse(r_l,r_t,dc):
    """Eqn 2 of Labrosse is correct. They then write dc/dt ~ dc/dTl * dTl/dt 
       and use the phase diagram in fig 2b to evaluate dc/dTl. This gives 
       dc = (1-0) > 0 and dT = (3500-5500) < 0. 
    """
    dT = 3500-5500
    num = r_t**3 - r_l**3
    den = 3 * r_t**2 * dc * dT
    return -num/den

def dcdt_fac(r_t, cl, Mbmo, drdt):
    fac = -4*np.pi*r_t**2 * rhom * cl / Mbmo * drdt
    return fac

def dadt_fac(dTadr, dTmdr): 
    """This is equation 25 of Gubbins et al 03, but in terms of dTdr rather than dT/dp. 
    Also I took the factor Ti/Tc = 1 as we are evaluating at the top, i.e. at Tc. """
    fac = 1 / (dTldr -dTadr)
    return fac

def getconc_bmo(Pc, dc):
    cFeO_bmo = dc / (1-Pc)
    cFeO_man = cFeO_bmo - dc
    return cFeO_bmo, cFeO_man

def mole2massconc_bmo(cbarFeO):

    AFeO = AO + AFe
    AMgO = AO + AMg
    cbarMgO = 1-cbarFeO
    Abar    = cbarFeO*AFeO + cbarMgO*AMgO
    cMgO    = cbarMgO*AMgO/Abar
    cFeO    = cbarFeO*AFeO/Abar    
    return cMgO, cFeO

def mass2moleconc_bmo(cMg, cFe):

    AFeO = AO + AFe
    AMgO = AO + AMg
    Abar   = 1 / (cMg/AMgO + cFe/AFeO)
    cbarMg = cMg / (AMgO/Abar)
    cbarFe = cFe / (AFeO/Abar)
    return cbarMg, cbarFe





