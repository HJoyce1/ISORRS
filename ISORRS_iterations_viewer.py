# =============================================================================
#  MODULE TO PLOT CONTOUR PLOTS REGARDING WHAT MODEL ITERATIONS ARE DOING
# Author: H.S. Joyce
# =============================================================================

# import modules
import numpy as np
import matplotlib.pyplot as plt

# load in npz file with data
data = np.load('output_alliter_jupiter_test_0.1u_carley.npz')

# set up figure
plt.figure(figsize=(10,20))
plt.subplot(6,1,1)
# countour plot -  need grid, iterations and then the data itself
# * 100 as turning fractional into perctange
# / 0.01 to represent the timestep (dt)
plt.contourf(data['z']/1e3, data['iterations'][:],
    ((data['n_its_hplus'][2:-2,:].T*100)/0.01),
    25)
h=plt.colorbar()
# set labels
h.set_label(r'$\delta n (\rm{H}^+)$ [%]')
plt.ylabel('Iteration Number')

# next subplot
plt.subplot(6,1,2)
# plotting
plt.contourf(data['z']/1e3, data['iterations'][:],
    ((data['rho_its_hplus'][2:-2,:].T*100)/0.01),
    25)
h=plt.colorbar()
h.set_label(r'$\delta \rho (\rm{H}^+)$ [%]')
plt.ylabel('Iteration Number')

plt.subplot(6,1,3)
# plotting
plt.contourf(data['z']/1e3, data['iterations'][:],
    ((data['P_its_hplus'][2:-2,:].T*100)/0.01),
    25)
h=plt.colorbar()
h.set_label(r'$\delta P (\rm{H}^+)$ [%]')
plt.ylabel('Iteration Number')

plt.subplot(6,1,4)
# plotting
plt.contourf(data['z']/1e3, data['iterations'][:],
    ((data['T_its_hplus'][2:-2,:].T*100)/0.01),
    25)
h=plt.colorbar()
h.set_label(r'$\delta T (\rm{H}^+)$ [%]')
plt.ylabel('Iteration Number')

plt.subplot(6,1,5)
# plotting
plt.contourf(data['z']/1e3, data['iterations'][:],
    ((data['u_its_hplus'][2:-2,:].T*100)/0.01),
    25)
h=plt.colorbar()
h.set_label(r'$\delta u (\rm{H}^+)$ [%]')
plt.ylabel('Iteration Number')

plt.subplot(6,1,6)
# plotting
plt.contourf(data['z']/1e3, data['iterations'][:],
    ((data['kappa_its_hplus'][2:-2,:].T*100)/0.01),
    25)
h=plt.colorbar()
h.set_label(r'$\delta \kappa (\rm{H}^+)$ [%]')
plt.ylabel('Iteration Number')
plt.xlabel('Distance Along Field Line [km]')