"""
FUNCTION FOR TESTING PLOTTING FUNCTIONS FOR RESULTS
"""
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import FormatStrFormatter

TPSdawn = np.array([1.3269169076616805e+27, 1.3268926817031633e+27, 1.3268684557446538e+27, 1.326844229786148e+27, 1.3268200038276453e+27, 1.3267957778691457e+27, 1.326771551910641e+27, 1.3267473259521603e+27, 1.3267230999936793e+27, 1.3266988740351894e+27, 1.3266746480767214e+27, 1.326650422118252e+27, 1.326626196159796e+27, 1.3266019702013427e+27, 1.326577744242901e+27, 1.3265535182844447e+27, 1.3265292923260067e+27, 1.326505066367555e+27, 1.326480840409135e+27, 1.3264566144507143e+27, 1.3264323884923013e+27])
TPSdusk = np.array([3.002329229514746e+27, 3.002273361606615e+27, 3.002217493698474e+27, 3.0021616257904025e+27, 3.002105757882263e+27, 3.002049889974125e+27, 3.0019940220660934e+27, 3.001938154158002e+27, 3.001882286249987e+27, 3.00182641834191e+27, 3.001770550433884e+27, 3.0017146825258864e+27, 3.0016588146178404e+27, 3.0016029467098427e+27, 3.001547078801826e+27, 3.00149121089384e+27, 3.0014353429858635e+27, 3.0013794750779247e+27, 3.001323607169997e+27, 3.001267739262044e+27, 3.0012118713541164e+27])

TMSdawn = np.array([2.467123388125489, 2.4671234357543432, 2.4671234833831583, 2.467123531011966, 2.4671235786407975, 2.4671236262696477, 2.467123673898458, 2.4671237215273374, 2.4671237691562005, 2.467123816784994, 2.4671238644139044, 2.4671239120428, 2.467123959671686, 2.4671240073005607, 2.467124054929477, 2.467124102558342, 2.467124150187279, 2.467124197816206, 2.467124245445116, 2.4671242930740798, 2.4671243407030294])
TMSdusk = np.array([5.58624645947698, 5.586246564437535, 5.586246669398034, 5.586246774358605, 5.586246879319105, 5.586246984279577, 5.586247089240341, 5.586247194200866, 5.586247299161577, 5.586247404122174, 5.586247509082873, 5.586247614043522, 5.58624771900414, 5.5862478239648405, 5.586247928925501, 5.586248033886162, 5.5862481388469325, 5.586248243807699, 5.586248348768541, 5.586248453729289, 5.586248558690098])

j = np.arange(0.1e-11,10.6e-11,0.5e-11)

TPSidawn =np.array([6.63460880142663e+26, 6.634608987225174e+26, 6.634609173023756e+26, 6.63460935882236e+26, 6.634609544620978e+26, 6.63460973041961e+26, 6.63460991621822e+26, 6.634610102016947e+26, 6.634610287815675e+26, 6.634610473614355e+26, 6.634610659413146e+26, 6.63461084521193e+26, 6.634611031010781e+26, 6.634611216809648e+26, 6.63461140260857e+26, 6.63461158840742e+26, 6.634611774206361e+26, 6.634611960005235e+26, 6.634612145804264e+26, 6.634612331603293e+26, 6.634612517402359e+26])
TPSidusk =np.array([1.5011702098094084e+27, 1.5011702511155196e+27, 1.5011702924216258e+27, 1.501170333727767e+27, 1.501170375033874e+27, 1.5011704163399815e+27, 1.5011704576461424e+27, 1.501170498952273e+27, 1.501170540258443e+27, 1.501170581564581e+27, 1.5011706228707448e+27, 1.5011706641769228e+27, 1.501170705483076e+27, 1.501170746789254e+27, 1.5011707880954226e+27, 1.5011708294016066e+27, 1.5011708707077947e+27, 1.5011709120140023e+27, 1.5011709533202152e+27, 1.5011709946264154e+27, 1.5011710359326283e+27])

TPSdawnn = (TPSdawn/1e27)
TPSduskk = (TPSdusk/1e27)
FAC = (j)/1e-4
FACs = np.array(FAC)

TPSidawnn = (TPSidawn/1e27)
TPSiduskk = (TPSidusk/1e27)


fig1,ax1=pl.subplots(2,1,sharex=True,figsize=(10,7))
pl.subplots_adjust(wspace=0,hspace=0)
e,f=np.polyfit(FACs, TPSdawnn, 1)
fitP=(e*FACs+f)
ax1[0].plot(FACs,fitP, color='goldenrod',linestyle='dashed')
ax1[0].plot(FACs,TPSdawnn,'s',label='Dawn')

c,d=np.polyfit(FACs,TPSduskk,1)
fitP2 = (c*FACs+d)
ax1[0].plot(FACs, fitP2, color='goldenrod',linestyle='dashed')
ax1[0].plot(FACs,TPSduskk,'^',label='Dusk')

ax1[0].legend(loc='best',bbox_to_anchor=(0.15,0.9))
ax1[0].tick_params(direction='in',bottom=True, top=True, left=True, right=True)
#ax1[0].yaxis.set_major_formatter(FormatStrFormatter('%.10f'))
ax1[0].set_xscale('log')
ax1[0].set_ylabel('Total Particle Source / $10^{27}$ (s$^{-1}$)')

# k,l=np.polyfit(FACs, TPSdawnn, 1)
# fitPi=(k*FACs+l)
# ax1[1].plot(FACs,fitPi, color='goldenrod',linestyle='dashed')
# ax1[1].plot(FACs,TPSdawnn,'s')

# m,n=np.polyfit(FACs,TPSduskk,1)
# fitPi2 = (m*FACs+n)
# ax1[1].plot(FACs, fitPi2, color='goldenrod',linestyle='dashed')
# ax1[1].plot(FACs,TPSduskk,'^')
# ax1[1].tick_params(direction='in',bottom=True, top=True, left=True, right=True)
# #ax1[1].yaxis.set_major_formatter(FormatStrFormatter('%.10f'))
# ax1[1].set_xscale('log')
# ax1[1].set_ylabel('Total Ion Particle Source / $10^{27}$ (s$^{-1}$)')


g,h=np.polyfit(FACs,TMSdawn,1)
fitM=(g*FACs+h)
ax1[1].plot(FACs,fitM, color='goldenrod',linestyle='dashed')
ax1[1].plot(FACs,TMSdawn,'s')

a,b=np.polyfit(FACs,TMSdusk,1)
fitM2=(a*FACs+b)
ax1[1].plot(FACs,fitM2, color='goldenrod',linestyle='dashed')
ax1[1].plot(FACs,TMSdusk,'^')

ax1[1].tick_params(direction='in',bottom=True, top=True, left=True, right=True)
#ax1[2].yaxis.set_major_formatter(FormatStrFormatter('%.10f'))
ax1[1].set_xscale('log')
ax1[1].set_xlabel('Current Strength / $(Am^{-2})$')
ax1[1].set_ylabel('Total Mass Source (kgs$^{-1}$)')


#folder = 'C:/Users/joyceh1/OneDrive - Lancaster University/pw_model/test_runs_asymmetries/tests/'
folder = '/Users/hannah/OneDrive - Lancaster University/pw_model/test_runs_asymmetries/tests/'

pl.savefig(folder+'FACs_logs_smol.pdf')
