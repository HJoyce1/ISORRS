"""
FUNCTION FOR TESTING PLOTTING FUNCTIONS FOR RESULTS
"""
import numpy as np
import matplotlib.pyplot as pl

TPSdawn = np.array([2.8042861962137288e+23, 5.6085723924274576e+23, 8.412858588641186e+23, 1.1217144784854915e+24, 1.4021430981068642e+24, 1.6825717177282372e+24, 1.96300033734961e+24, 2.243428956970983e+24, 2.5238575765923557e+24, 2.8042861962137284e+24, 3.084714815835101e+24, 3.3651434354564743e+24, 3.645572055077847e+24, 3.92600067469922e+24, 4.2064292943205923e+24, 4.486857913941966e+24, 4.767286533563339e+24, 5.047715153184711e+24, 5.328143772806085e+24, 5.608572392427457e+24])

TPSdusk = np.array([2.1559479901600134e+23, 4.311895980320027e+23, 6.467843970480039e+23, 8.623791960640054e+23, 1.0779739950800066e+24, 1.2935687940960078e+24, 1.509163593112009e+24, 1.7247583921280107e+24, 1.940353191144012e+24, 2.155947990160013e+24, 2.3715427891760144e+24, 2.5871375881920155e+24, 2.802732387208017e+24, 3.018327186224018e+24, 3.233921985240019e+24, 3.4495167842560214e+24, 3.665111583272023e+24, 3.880706382288024e+24, 4.0963011813040254e+24, 4.311895980320026e+24])

TMSdawn = np.array([0.000493334367169666, 0.000986668734339332, 0.0014800031015089977, 0.001973337468678664, 0.00246667183584833, 0.0029600062030179955, 0.0034533405701876623, 0.003946674937357328, 0.004440009304526994, 0.00493334367169666, 0.005426678038866326, 0.005920012406035991, 0.006413346773205657, 0.006906681140375325, 0.00740001550754499, 0.007893349874714656, 0.008386684241884322, 0.008880018609053988, 0.009373352976223655, 0.00986668734339332])

TMSdusk = np.array([0.0003791508183182188, 0.0007583016366364376, 0.0011374524549546564, 0.0015166032732728753, 0.0018957540915910938, 0.0022749049099093127, 0.002654055728227532, 0.0030332065465457506, 0.0034123573648639693, 0.0037915081831821876, 0.004170659001500406, 0.0045498098198186255, 0.004928960638136844, 0.005308111456455064, 0.005687262274773282, 0.006066413093091501, 0.006445563911409719, 0.006824714729727939, 0.007203865548046157, 0.007583016366364375])

widths = np.array([0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,4.75,5.0])

TPSdawnn = (TPSdawn/1e24)
TPSduskk = (TPSdusk/1e24)
#FACs = (FAC*1e)

TPSidawn = np.array([3.235889741584777e+23, 6.471779483169554e+23, 9.707669224754331e+23, 1.2943558966339108e+24, 1.6179448707923885e+24, 1.9415338449508662e+24, 2.2651228191093442e+24, 2.5887117932678216e+24, 2.9123007674262993e+24, 3.235889741584777e+24, 3.559478715743255e+24, 3.8830676899017324e+24, 4.20665666406021e+24, 4.5302456382186884e+24, 4.853834612377166e+24, 5.177423586535643e+24, 5.501012560694122e+24, 5.824601534852599e+24, 6.148190509011077e+24, 6.471779483169554e+24])
TPSidusk = np.array([2.2595600648886565e+23, 4.519120129777313e+23, 6.778680194665968e+23, 9.038240259554626e+23, 1.1297800324443281e+24, 1.3557360389331936e+24, 1.581692045422059e+24, 1.8076480519109252e+24, 2.0336040583997907e+24, 2.2595600648886562e+24, 2.4855160713775217e+24, 2.711472077866387e+24, 2.937428084355253e+24, 3.163384090844118e+24, 3.389340097332984e+24, 3.6152961038218504e+24, 3.8412521103107153e+24, 4.0672081167995814e+24, 4.293164123288447e+24, 4.5191201297773124e+24])

TPSidawnn = (TPSidawn/1e24)
TPSiduskk = (TPSidusk/1e24)


fig1,ax1=pl.subplots(2,1,sharex=True,figsize=(10,7))
pl.subplots_adjust(wspace=0,hspace=0)
# temperature
e,f=np.polyfit(widths, TPSidawnn, 1)
fitP=(e*widths+f)
#ax1[0].set_yscale('log')
ax1[0].plot(widths,fitP, color='lightseagreen',linestyle='dashed')
ax1[0].plot(widths,TPSidawnn,'s',color='orange',label='Dawn')

c,d=np.polyfit(widths,TPSiduskk,1)
fitP2 = (c*widths+d)
ax1[0].plot(widths, fitP2, color='lightseagreen',linestyle='dashed')
ax1[0].plot(widths,TPSiduskk,'^',color='blue',label='Dusk')

ax1[0].legend(loc='best',bbox_to_anchor=(0.15,0.97))
ax1[0].tick_params(direction='in',bottom=True, top=True, left=True, right=True)
ax1[0].set_ylabel('Total Particle Source / $10^{24}$ (s$^{-1}$)')
ax1[0].set_ylim(0, 7.5)

g,h=np.polyfit(widths,TMSdawn,1)
fitM=(g*widths+h)
#ax1[1].set_yscale('log')
ax1[1].plot(widths,fitM, color='lightseagreen',linestyle='dashed')
ax1[1].plot(widths,TMSdawn,'s',color='orange')

a,b=np.polyfit(widths,TMSdusk,1)
fitM2=(a*widths+b)
ax1[1].plot(widths,fitM2, color='lightseagreen',linestyle='dashed')
ax1[1].plot(widths,TMSdusk,'^',color='blue')

ax1[1].tick_params(direction='in',bottom=True, top=True, left=True, right=True)
ax1[1].set_xlabel('Auroral Width \N{DEGREE SIGN}')
ax1[1].set_ylabel('Total Mass Source (kgs$^{-1}$)')
ax1[1].set_ylim(0.000, 0.011)


folder = '/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/march/asym_runs/plots/'

pl.savefig(folder+'widthss_plot.png')


# ---------------

# import numpy as np
# import matplotlib.pyplot as pl
# from scipy.optimize import curve_fit

# def func(x, a, b, c):
#     return a * np.exp(b * x) + c

# TPSdawn = ([7.139646171608408e+23, 7.32642856441471e+23, 7.525260715141202e+23, 7.889869963232647e+23, 9.240980990286932e+23, 1.4548188105861335e+24, 3.159406449951686e+24, 7.555806518781051e+24, 1.7124261847511583e+25, 3.549132693769023e+25, 6.747410486983394e+25, 1.189862494379908e+26])
# TPSdusk = ([1.932790676407711e+24, 1.975861585696706e+24, 2.0216973959114256e+24, 2.1056480298788877e+24, 2.4165741050869248e+24, 3.638789696714353e+24, 7.568913015790242e+24, 1.770796145509749e+25, 3.977345048668259e+25, 8.212805465439497e+25, 1.558801858391675e+26, 2.746667779071461e+26])

# TMSdawn = np.array([0.0021455525059934843, 0.0022061612265193483, 0.002264729055185576, 0.0023214441242692214, 0.00237646723825049, 0.0024299373618188454, 0.0024819780612897182, 0.0025327171151694984, 0.002582365234794386, 0.002631476905850125, 0.0026816380537359893, 0.0027369567719293768])
# TMSdusk = np.array([0.0049477028232744955, 0.005087465922951312, 0.005222522747239571, 0.00535330712713752, 0.005480189876025662, 0.005603491444546924, 0.005723496770595885, 0.005840500478429421, 0.0059549883173403454, 0.006068238395357731, 0.006183906329357668, 0.006311461770583532])
# temps = np.array([ 500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600])

# TPSidawn = [1.2952339455793118e+24, 1.330760247196404e+24, 1.3651323147190636e+24, 1.3984526656901446e+24, 1.43080941444926e+24, 1.462279206554217e+24, 1.492931939704421e+24, 1.5228505332746796e+24, 1.5522165541179785e+24, 1.5815988254769597e+24, 1.612719818879989e+24, 1.6501273780605366e+24]
# TPSidusk = [2.9868092522625464e+24, 3.0687325416153937e+24, 3.1479941794076423e+24, 3.2248305699501933e+24, 3.299444903918329e+24, 3.372013924807197e+24, 3.442698812585161e+24, 3.511690755935488e+24, 3.579408340460897e+24, 3.6471628909687876e+24, 3.718925338529395e+24, 3.805180419611454e+24]

# TPSdawnn = np.array(TPSdawn)/1e24
# TPSduskk = np.array(TPSdusk)/1e24

# TPSidawnn = np.array(TPSidawn)/1e24
# TPSiduskk = np.array(TPSidusk)/1e24

# # popt, pcov = curve_fit(func, temps, TPSdawnn)
# # popt2, pcov2 = curve_fit(func, temps, TPSduskk)


# fig1,ax1=pl.subplots(2,1,sharex=True,figsize=(10,7))
# pl.subplots_adjust(wspace=0,hspace=0)

# # f=np.polyfit(temps, np.log(TPSdawnn), 1, w=np.sqrt(TPSdawnn))
# # a=np.exp(f[1])
# # b=f[0]
# # fit_x = np.linspace(np.min(temps),np.max(temps),12)
# # fit = a * np.exp(b*fit_x)


# # g = np.polyfit(temps,np.log(TPSduskk),1)
# # c=np.exp(g[1])
# # d=g[0]
# # fit1 = a*np.exp(b*temps)

# # def polyfit(y):
# #     fit = np.polyfit(temps,y,4)
# #     a = fit[0]
# #     b = fit[1]
# #     c = fit[2]
# #     d= fit[3]
# #     e=fit[4]
# #     fit_eq = a * temps**4 + b*temps**3 +c*temps**2 + d**temps + e
# #     return(fit_eq)

# # temperature
# e,f=np.polyfit(temps,TPSidawnn, 1)
# fitP=(e*temps+f)
# #ax1[0].set_yscale('log')
# #ax1[0].plot(temps,func(temps, *popt),color='green',linestyle='dashed')
# ax1[0].plot(temps,fitP, color='salmon',linestyle='dashed')
# ax1[0].plot(temps,TPSidawnn,'s',color='orange',label='Dawn')
# #ax1.plot(temps,polyfit(TPSidawnn),color='salmon',linestyle='dashed')

# c,d=np.polyfit(temps,TPSiduskk,1)
# fitP2 = (c*temps+d)
# ax1[0].plot(temps, fitP2, color='salmon',linestyle='dashed')
# #ax1[0].plot(temps,fit1,color='salmon',linestyle='dashed')
# ax1[0].plot(temps,TPSiduskk,'^',color='blue',label='Dusk')
# #ax1.plot(temps,polyfit(TPSiduskk),color='salmon',linestyle='dashed')

# #ax1[0].plot(temps,func(temps, *popt2),color='salmon',linestyle='dashed')
# ax1[0].legend(loc='best',bbox_to_anchor=(0.15,0.97))
# #ax1.set_xlim(650, 1050)
# #ax1[0].set_ylim(0,5)
# ax1[0].tick_params(direction='in',bottom=True, top=True, left=True, right=True)
# ax1[0].set_ylim(0.6, 4.5)
# ax1[0].set_ylabel('Total Particle Source / $10^{24}$ (s$^{-1}$)')
# ax1[0].set_xlabel('Ionospheric Temperature (K)')


# g,h=np.polyfit(temps,TMSdawn,1)
# fitM=(g*temps+h)
# ax1[1].plot(temps,fitM, color='salmon',linestyle='dashed')
# ax1[1].plot(temps,TMSdawn,'s',color='orange')
# #ax1[1].plot(temps,polyfit(TMSdawn),color='salmon',linestyle='dashed')


# a,b=np.polyfit(temps,TMSdusk,1)
# fitM2=(a*temps+b)
# ax1[1].plot(temps,fitM2, color='salmon',linestyle='dashed')
# ax1[1].plot(temps,TMSdusk,'^',color='blue')
# #ax1[1].plot(temps,polyfit(TMSdusk),color='salmon',linestyle='dashed')
# ax1[1].tick_params(direction='in',bottom=True, top=True, left=True, right=True)
# ax1[1].set_xlabel('Ionospheric Temperature (K)')
# ax1[1].set_ylabel('Total Mass Source (kgs$^{-1}$)')
# ax1[1].set_ylim(0.001, 0.0075)

# folder = '/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/march/asym_runs/plots/'

# pl.savefig(folder+'temps_plot_ions_50k.pdf')
