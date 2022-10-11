"""
FUNCTION FOR TESTING PLOTTING FUNCTIONS FOR RESULTS
"""
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import FormatStrFormatter

TPSdawn = np.array([9.378228582266759e+26, 9.574350206212505e+26, 9.765635224164507e+26, 9.952415185408572e+26, 1.013498524436655e+27, 1.0313609549535785e+27, 1.0488525685202299e+27, 1.0659948529577595e+27, 1.0828074124672676e+27, 1.0993084565324064e+27, 1.1155155330180039e+27, 1.1314466618048364e+27, 1.1471219925808536e+27, 1.1625660331201368e+27, 1.1778103940469108e+27, 1.1928969011409388e+27])
TPSdusk = np.array([2.162912649961894e+27, 2.2081381099199205e+27, 2.2522482531340358e+27, 2.2953195345543013e+27, 2.337420016509606e+27, 2.3786106113924076e+27, 2.4189461059130104e+27, 2.458476050782309e+27, 2.497245653100349e+27, 2.5352969037595167e+27, 2.572670267960169e+27, 2.609407299357317e+27, 2.645554462906072e+27, 2.681168273245558e+27, 2.7163216240926057e+27, 2.7511109651838644e+27])

TMSdawn = np.array([1.621847370352647, 1.6538046758588312, 1.68500964029057, 1.7155117518712037, 1.7453552778237917, 1.7745800084289014, 1.8032218735181362, 1.8313134734121013, 1.8588845867092878, 1.8859627543584985, 1.9125740757055283, 1.9387443617334084, 1.964500756586442, 1.9898738642122762, 2.0149003235306053, 2.0396256893016007])
TMSdusk = np.array([3.7399646472638195, 3.813657889913081, 3.885616241195935, 3.9559538204847193, 4.024772707574841, 4.0921646584753155, 4.15821252703583, 4.222991489285535, 4.286570214387573, 4.349012211502746, 4.410377665461638, 4.470726096119886, 4.530120097569001, 4.588630242144806, 4.646341018733112, 4.7033574760189385])
temps = np.array([ 500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000])

TPSdawnn = (TPSdawn/1e27)
TPSduskk = (TPSdusk/1e27)


fig1,ax1=pl.subplots(2,1,sharex=True,figsize=(10,7))
pl.subplots_adjust(wspace=0,hspace=0)
# temperature
e,f=np.polyfit(temps, TPSdawnn, 1)
fitP=(e*temps+f)
ax1[0].plot(temps,fitP, color='salmon',linestyle='dashed')
ax1[0].plot(temps,TPSdawnn,'s',label='Dawn')

c,d=np.polyfit(temps,TPSduskk,1)
fitP2 = (c*temps+d)
ax1[0].plot(temps, fitP2, color='salmon',linestyle='dashed')
ax1[0].plot(temps,TPSduskk,'^',label='Dusk')

ax1[0].legend(loc='best',bbox_to_anchor=(0.15,0.97))
ax1[0].tick_params(direction='in',bottom=True, top=True, left=True, right=True)
ax1[0].set_ylabel('Total Particle Source / $10^{27}$ (s$^{-1}$)')


g,h=np.polyfit(temps,TMSdawn,1)
fitM=(g*temps+h)
ax1[1].plot(temps,fitM, color='salmon',linestyle='dashed')
ax1[1].plot(temps,TMSdawn,'s')

a,b=np.polyfit(temps,TMSdusk,1)
fitM2=(a*temps+b)
ax1[1].plot(temps,fitM2, color='salmon',linestyle='dashed')
ax1[1].plot(temps,TMSdusk,'^')

ax1[1].tick_params(direction='in',bottom=True, top=True, left=True, right=True)
ax1[1].set_xlabel('Ionospheric Temperature (K)')
ax1[1].set_ylabel('Total Mass Source (kgs$^{-1}$)')

folder = '/Users/hannah/OneDrive - Lancaster University/pw_model/test_runs_asymmetries/tests/'
#folder = 'C:/Users/joyceh1/OneDrive - Lancaster University/pw_model/test_runs_asymmetries/tests/'

pl.savefig(folder+'temps_plot_density_test.png')
