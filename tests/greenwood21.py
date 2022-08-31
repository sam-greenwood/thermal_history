"""
Test the model setup of Greenwood et al. (2021) Mars study.
Uses the 'thiriet19' mantle model coupled to 'leeds' and 'leeds_thermal' core/stable layer models.
greenwood21_params.py contains the parameters for the reference case in Figure 1 of Greenwood et al. (2021)
"""

from thermal_history.model import Parameters, setup_model
import matplotlib.pyplot as plt

prm = Parameters('greenwood21_params.py')

model = setup_model(prm, core_method='leeds', mantle_method='thiriet19', stable_layer_method='leeds_thermal')

dt = 1e6*prm.ys
for i in range(4500):
    model.evolve(dt)