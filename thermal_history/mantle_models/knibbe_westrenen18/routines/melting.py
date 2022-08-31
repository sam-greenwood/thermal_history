import numpy as np
# from scipy.integrate import trapezoid
from thermal_history.utils.optimised_funcs import linspace, trapezoid
from . import profiles as prof
import matplotlib.pyplot as plt

def mantle_melt(D_ref, D_crust, r, rho, g, T, solidus_params, liquidus_params):

    #Increase grid resolution to accurately identify melting zones
    r_fine = linspace(r[0], r[-1], 100)
    rho_fine = np.interp(r_fine, r, rho)
    P_fine = prof.pressure(r_fine, rho_fine, g)
    T_fine = np.interp(r_fine, r, T)

    T_sol = prof.solidus(P_fine, solidus_params) + (D_crust/D_ref)* 150 
    T_liq = prof.liquidus(P_fine, liquidus_params) 

    mask = np.zeros(r_fine.size)
    mask[np.where(T_fine > T_sol)[0]] = 1

    #Keep melt fraction at <= 1
    T_fine[T_fine>T_liq] = T_liq[T_fine>T_liq]

    # V_melt = 4*np.pi*trapezoid(mask*r**2,x=r)
    V_melt = 4*np.pi*trapezoid(r_fine, mask*r_fine**2)[-1]

    if V_melt == 0:
        ma = 0
        dm_dT = 0
        melting_radii = np.array([np.nan, np.nan])
    else:
        # ma = 1/V_melt * 4*np.pi*trapezoid( mask*r**2 *(T-T_sol) / (T_liq-T_sol) , x=r)
        ma = 1/V_melt * 4*np.pi*trapezoid(r_fine, mask*r_fine**2 *(T_fine-T_sol) / (T_liq-T_sol))[-1]
        melting_radii = np.array([r_fine[T_fine>T_sol][0], r_fine[T_fine>T_sol][-1]])

        #Increase T by 1 degree and recalc ma to get dma_dT
        T2 = np.interp(r_fine, r, T+1)
        mask2 = np.zeros(r_fine.size)
        mask2[np.where(T2 > T_sol)[0]] = 1
        T2[T2>T_liq] = T_liq[T2>T_liq]
        V_melt2 = 4*np.pi*trapezoid(r_fine, mask*r_fine**2)[-1]
        ma2 = 1/V_melt2 * 4*np.pi*trapezoid(r_fine, mask*r_fine**2 *(T2-T_sol) / (T_liq-T_sol))[-1]
        dm_dT = ma2-ma

    #Save radial profiles
    profiles = {'r': r_fine,
                'rho': rho_fine,
                'P': P_fine,
                'T': T_fine,
                'T_liq': T_liq,
                'T_sol': T_sol}


    return V_melt, ma, melting_radii, dm_dT, profiles
