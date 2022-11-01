import numpy as np
from thermal_history.stable_layer_models.leeds_thermal.routines.diffusion import diffusion
from thermal_history.utils.optimised_funcs import trapezoid

def test_diffusion():

    #Check diffusion solution conserves energy to sufficient accuracy
    r = np.linspace(0,100, 100)
    k = np.ones(r.size)
    rho = np.ones(r.size)
    cp = np.ones(r.size)
    T = np.linspace(0,1,r.size)

    D_T = k/(rho*cp)
    dk_dr = 0

    lb_type, lb = 1, 0
    ub_type, ub = 1, 0

    dt = 1

    T_new = diffusion(T, r, dt, D_T, k, dk_dr, (lb_type,lb),(1,ub))

    dT = T_new-T
    dE = 4*np.pi*trapezoid(r, r**2 * rho * cp * dT)[-1]
    thermal_mass = 4*np.pi*trapezoid(r, r**2 * rho * cp * T)[-1]

    assert np.abs(dE/thermal_mass) < 1e-5, 'Secular cooling exceeds 0.001% of total thermal mass'

if __name__ == '__main__':
    test_diffusion()