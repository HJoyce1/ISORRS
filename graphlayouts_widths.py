"""
FUNCTION FOR TESTING PLOTTING FUNCTIONS FOR RESULTS
"""
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import FormatStrFormatter

TPSdawn = np.array([3.1892848095203414e+26, 6.378569619040683e+26, 9.567854428561023e+26, 1.2757139238081366e+27, 1.5946424047601706e+27, 1.9135708857122047e+27, 2.232499366664239e+27, 2.551427847616273e+27, 2.8703563285683074e+27, 3.189284809520341e+27, 3.508213290472375e+27, 3.8271417714244094e+27, 4.146070252376443e+27, 4.464998733328478e+27, 4.7839272142805113e+27, 5.102855695232546e+27, 5.42178417618458e+27, 5.740712657136615e+27, 6.059641138088649e+27, 6.378569619040682e+27])
TPSdusk = np.array([2.7798413585927282e+26, 5.5596827171854564e+26, 8.339524075778183e+26, 1.1119365434370913e+27, 1.389920679296364e+27, 1.6679048151556366e+27, 1.94588895101491e+27, 2.2238730868741826e+27, 2.5018572227334554e+27, 2.779841358592728e+27, 3.0578254944520005e+27, 3.3358096303112733e+27, 3.613793766170546e+27, 3.89177790202982e+27, 4.169762037889092e+27, 4.447746173748365e+27, 4.725730309607638e+27, 5.003714445466911e+27, 5.281698581326183e+27, 5.559682717185456e+27])

TMSdawn = np.array([0.5930585891954513, 1.1861171783909026, 1.7791757675863542, 2.372234356781805, 2.9652929459772572, 3.5583515351727084, 4.1514101243681605, 4.74446871356361, 5.337527302759062, 5.9305858919545145, 6.523644481149965, 7.116703070345417, 7.709761659540868, 8.302820248736321, 8.89587883793177, 9.48893742712722, 10.081996016322673, 10.675054605518124, 11.268113194713576, 11.861171783909029])
TMSdusk = np.array([0.517245060037595, 1.03449012007519, 1.5517351801127843, 2.06898024015038, 2.586225300187974, 3.1034703602255687, 3.6207154202631644, 4.13796048030076, 4.655205540338354, 5.172450600375948, 5.689695660413543, 6.206940720451137, 6.724185780488733, 7.241430840526329, 7.758675900563922, 8.27592096060152, 8.793166020639111, 9.310411080676708, 9.827656140714303, 10.344901200751895])

widths = np.array([0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,4.75,5.0])

TPSdawnn = (TPSdawn/1e27)
TPSduskk = (TPSdusk/1e27)


fig1,ax1=pl.subplots(2,1,sharex=True,figsize=(10,7))
pl.subplots_adjust(wspace=0,hspace=0)
# temperature
e,f=np.polyfit(widths, TPSdawnn, 1)
fitP=(e*widths+f)
ax1[0].plot(widths,fitP, color='lightseagreen',linestyle='dashed')
ax1[0].plot(widths,TPSdawnn,'s',label='Dawn')

c,d=np.polyfit(widths,TPSduskk,1)
fitP2 = (c*widths+d)
ax1[0].plot(widths, fitP2, color='lightseagreen',linestyle='dashed')
ax1[0].plot(widths,TPSduskk,'^',label='Dusk')

ax1[0].legend(loc='best',bbox_to_anchor=(0.15,0.97))
ax1[0].tick_params(direction='in',bottom=True, top=True, left=True, right=True)
ax1[0].set_ylabel('Total Particle Source / $10^{27}$ (s$^{-1}$)')


g,h=np.polyfit(widths,TMSdawn,1)
fitM=(g*widths+h)
ax1[1].plot(widths,fitM, color='lightseagreen',linestyle='dashed')
ax1[1].plot(widths,TMSdawn,'s')

a,b=np.polyfit(widths,TMSdusk,1)
fitM2=(a*widths+b)
ax1[1].plot(widths,fitM2, color='lightseagreen',linestyle='dashed')
ax1[1].plot(widths,TMSdusk,'^')

ax1[1].tick_params(direction='in',bottom=True, top=True, left=True, right=True)
ax1[1].set_xlabel('Auroral Width \N{DEGREE SIGN}')
ax1[1].set_ylabel('Total Mass Source (kgs$^{-1}$)')


folder = 'C:/Users/joyceh1/OneDrive - Lancaster University/pw_model/test_runs_asymmetries/tests/'

pl.savefig(folder+'widths_plot_neutrals.png')
