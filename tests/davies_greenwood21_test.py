"""
Runs the combined core and mantle model of Davies and Greenwood (2021). Note that a small correction has been since fixed in the code,
resulting in the entropy of thermal conduction (Ek) to be fixed in time as it should be for an adiabatic temperature profile. The entropy
available for the magnetic field (Ej) is therefore slightly lower than published. 
"""

def test_davies_greenwood21():

    from thermal_history.model import Parameters, setup_model

    prm = Parameters('davies_greenwood21_params.py')

    model = setup_model(prm, core_method='leeds', stable_layer_method='leeds_chemical', mantle_method='driscol_bercovici14')

    while model.time < 4.5e9*prm.ys:
        #Need v. small timestep for stable solution at P=1. I found it was needed for first ~50,000 iterations before moving to a larger timestep.
        if prm.partition_FeO == 1 and model.it < 50000:
            dt = 0.01e6*prm.ys
        else:
            dt = 0.1e6*prm.ys

        model.evolve(dt)

    #Check final values are in the right ballpark.
    assert (model.core.ri > 1100e3 and model.core.ri < 1300e3), f"Inner core is too far from present radius, model.core.ri = {model.core.ri/1000: .2f} km"

    assert model.stable_layer.layer_thickness > 80e3, f"Stable layer is not large enough, only {model.stable_layer.layer_thickness/1000: .2f} km"

if __name__ == '__main__':
    test_davies_greenwood21()
