import matplotlib.pyplot as plt
import sys
sys.path.append('../')
from thermal_history.model import Parameters, setup_model

prm = Parameters('simple_test_params.py')
model = setup_model(prm, core_method='simple_test')

dt = 1e6*prm.ys   #1 Myr time step
model.mantle.Q_cmb = 1e12 #1TW

for i in range(1000):
    model.evolve(dt)


data = model.save_dict_to_numpy_array()

time = data['core']['time']/(prm.ys*1e6)
T    = data['core']['T']

plt.plot(time, T)
plt.xlabel('Time [Myrs]')
plt.ylabel('Temperature [K]')
plt.show()
