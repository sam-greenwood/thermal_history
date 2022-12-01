
import numpy as np

#Run tests on melting curves
def test_alfe():

    from thermal_history.core_models.leeds.routines.chemistry import iron_melting, iron_melting_gradient

    melting_params = np.array([1,-1e-9,0,0]) #Pressure Polynomials
    P = np.array([0, 1e9]) #0 and 1 GPa

    Tm = iron_melting(P, melting_params)
    assert Tm[0] == 1 and Tm[-1] == 0, f"iron_melting should give [1,0], not: {Tm}"

    Tm_grad = iron_melting_gradient(P, melting_params)
    assert np.unique(Tm_grad) == -1e-9, f"iron_melting_gradient should give [-1e-9, -1e-9], not : {Tm_grad}"


def test_simon_glatzel():

    from thermal_history.core_models.leeds.routines.chemistry import simon_glatzel, simon_glatzel_gradient

    melting_params = np.array([1,2e9,2])

    P = np.array([6e9])

    Tm   = simon_glatzel(P, melting_params)
    assert Tm[0] == 2, f"iron_melting should give [2], not: {Tm}"

    Tm_grad = simon_glatzel_gradient(P, melting_params)
    assert Tm_grad*1e9 == 0.125, f"iron_melting_gradient should give [0.125e-9], not : {Tm_grad}"


if __name__ == '__main__':
    test_alfe()
    test_simon_glatzel()