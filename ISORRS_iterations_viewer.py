import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors

data = np.load('output_alliter_jupiter_test_stability_10000.npz')

plt.figure(figsize=(10,20))
plt.subplot(6,1,1)
plt.contourf(data['z']/1e3, data['iterations'][:],
    ((data['n_its_hplus'][2:-2,:])*100).T,
    25)
h=plt.colorbar()
h.set_label(r'$\delta n (\rm{H}^+)$ [%]')
plt.ylabel('Iteration Number')

plt.subplot(6,1,2)
plt.contourf(data['z']/1e3, data['iterations'][:],
    ((data['rho_its_hplus'][2:-2,:])*100).T,
    25)
h=plt.colorbar()
h.set_label(r'$\delta \rho (\rm{H}^+)$ [%]')
plt.ylabel('Iteration Number')

plt.subplot(6,1,3)
plt.contourf(data['z']/1e3, data['iterations'][:],
    ((data['P_its_hplus'][2:-2,:])*100).T,
    25)
h=plt.colorbar()
h.set_label(r'$\delta P (\rm{H}^+)$ [%]')
plt.ylabel('Iteration Number')

plt.subplot(6,1,4)
plt.contourf(data['z']/1e3, data['iterations'][:],
    ((data['T_its_hplus'][2:-2,:])*100).T,
    25)
h=plt.colorbar()
h.set_label(r'$\delta T (\rm{H}^+)$ [%]')
plt.ylabel('Iteration Number')

plt.subplot(6,1,5)
plt.contourf(data['z']/1e3, data['iterations'][:],
    ((data['u_its_hplus'][2:-2,:])*100).T,
    25)
h=plt.colorbar()
h.set_label(r'$\delta u (\rm{H}^+)$ [%]')
plt.ylabel('Iteration Number')

plt.subplot(6,1,6)
plt.contourf(data['z']/1e3, data['iterations'][:],
    ((data['kappa_its_hplus'][2:-2,:])*100).T,
    25)
h=plt.colorbar()
h.set_label(r'$\delta \kappa (\rm{H}^+)$ [%]')
plt.ylabel('Iteration Number')
plt.xlabel('Distance Along Field Line [km]')