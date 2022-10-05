"""
Runs the mantle model of Driscol and Bercovici (2014) (minus mantle melting) coupled to the leeds core model.
Although implemented in a different way, the core model should exhibit equivalent behaviour.
Model parameters for figure 5 (Case E). Note since mantle melting is not included, this model does not reproduce
their results and still produces a thermal catastrophe. Only the last 1 Ga is modelled just to verify no errors are
produced by the code before the catastophe occurs.
"""

from tests.driscol_bercovici14_params import Q_surface

def test_driscol_bercovici14():

    from thermal_history.model import Parameters, setup_model

    prm = Parameters('driscol_bercovici14_params.py')

    model = setup_model(prm, core_method='leeds', mantle_method='driscol_bercovici14')
    
    model.time = 4.5e9*prm.ys #Start from present day
    dt = -1e6*prm.ys          #Iterate backwards in time.

    for i in range(1000):
        model.evolve(dt, print_freq=100)

    assert model.core.ri == 0, 'Inner core should have melted by 1 Ga.'

if __name__ == '__main__':
    test_driscol_bercovici14()