#Module for plotting functions to tidy code up
# I think you will need to change some parameters depending on which planet
# you want to plot out. I think these can replace the plotting in the
# main functions if you want to


#Personalised plotting function by me for me just to quickly check on things
def plot_me_quick(x,y,xlabel,ylabel):
    import matplotlib.pyplot as pl
    pl.plot(x,y)
    pl.ylabel(ylabel)
    pl.xlabel(xlabel)
    #pl.grid(grid) 
    
def plot_me_log(x,y,xlabel,ylabel):
    import matplotlib.pyplot as pl
    pl.plot(x,y)
    pl.ylabel(ylabel)
    pl.xlabel(xlabel)
    pl.yscale('log')
    #pl.grid(grid)
 
    #plot specifically for results
def results_plot(z,z_ext,radius,num_ion,e_charge,E,ion_dict,electron_dict,ac,ag,e_flux,ion_flux):
    import matplotlib.pyplot as pl
    pl.figure(figsize=(8,10))
    
    ax1 = pl.subplot(4,1,1)
    pl.plot(z/1000,(E)/1e-7,linestyle='-',color='black')
    pl.ylabel('a) $E_{\parallel}$ $(V$ $m^{-1}) $x$ 10^{-7}$')
    pl.xlim([0,z_ext[-1]/1000])
    ax1.set_xlabel('')
    ax1.xaxis.set_ticklabels([])
    ax1.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    #pl.grid('on') 
    ax2 = ax1.twiny()
    ax2.plot(z/radius,(E)/1e-7,linestyle= '-',color=[0,0,0])
    ax2.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax2.set_xlim([0,z_ext[-1]/radius])
    ax2.tick_params(direction='in')
    
    ax3 = pl.subplot(4,1,2)
    pl.plot(z/1000,((e_charge/ion_dict[1]["mass"])*(E)),linestyle= '-',color='blue')
    pl.plot(z/1000,((e_charge/ion_dict[2]["mass"])*(E)),linestyle= '-',color='lightskyblue')
    pl.plot(z/1000,ac,linestyle= '--',color='peru')
    pl.plot(z/1000,ag,linestyle= '-.',color='saddlebrown')
    pl.ylabel('b) Acceleration Terms \n $(m$ $s^{-2})$')
    pl.xlim([0,z_ext[-1]/1000])
    ax3.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    #pl.grid('on') 
    ax4 = ax3.twiny()
    ax4.plot(z/radius,(e_charge/ion_dict[1]["mass"])*(E),linestyle='-',color='blue')
    ax4.set_xlim([0,z_ext[-1]/radius])
    ax3.set_xlabel('')
    ax3.xaxis.set_ticklabels([])
    ax4.set_xlabel('')
    ax4.xaxis.set_ticklabels([])
    ax4.tick_params(direction='in')
    ax3.legend(['$H^+$ Electric Field','$H_3^+$ Electric Field','Centrifugal', 'Gravitational'],loc='upper right',bbox_to_anchor=(0.99,0.98))
            
        
    ax5 = pl.subplot(4,1,3)
    pl.yscale('log')
    pl.plot(z/1000,e_flux[2:-2,-1],color='orange')
    pl.ylabel('c) Electron Flux \n $(m^{-2}$ $s^{-1})$ x A $(m^2)$')
    pl.xlim([0,z_ext[-1]/1000])
    ax5.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    #pl.grid('on') 
    ax6 = ax5.twiny()
    ax6.plot(z/radius,e_flux[2:-2,-1],color='orange')
    ax6.set_xlim([0,z_ext[-1]/radius])
    ax5.set_xlabel('')
    ax5.xaxis.set_ticklabels([])
    ax6.set_xlabel('')
    ax6.tick_params(direction='in')
    ax6.xaxis.set_ticklabels([])
    ax5.legend(['$e^-$'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    ax7 = pl.subplot(4,1,4)
    pl.yscale('log')
    pl.plot(z/1000,ion_flux[2:-2,-1,0],color='blue')
    pl.plot(z/1000,ion_flux[2:-2,-1,1],color='lightskyblue')    
    pl.ylabel('d) Ion Flux \n $(m^{-2}$ $s^{-1})$ x A $(m^2)$')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    ax7.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    #pl.grid('on') 
    ax8 = ax7.twiny()
    ax8.plot(z/radius,ion_flux[2:-2,-1,0],color='blue')
    ax8.set_xlabel('')
    ax8.set_xlim([0,z_ext[-1]/radius])
    ax8.xaxis.set_ticklabels([])
    ax8.tick_params(direction='in')
    ax7.legend(['$H^+$','$H_3^+$'],loc='upper right',bbox_to_anchor=(0.99,0.95))

        
    pl.subplots_adjust(hspace=0.0,wspace=0.5)
    pl.suptitle("Results")

def species_plot(z,z_ext,ion_dict,radius):
    import matplotlib.pyplot as pl
    pl.figure(figsize=(15,8))
    ax1 = pl.subplot(2,2,1)
    pl.yscale('log')
    pl.plot(z/1000,ion_dict["n"][2:-2,-1])
    pl.ylabel('Density (kg/m^3)')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    #pl.grid('on') 
    ax2 = ax1.twiny()
    ax2.plot(z/radius,ion_dict["n"][2:-2,-1])
    ax2.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax2.set_xlim([0,z_ext[-1]/radius])
    
    ax3 =pl.subplot(2,2,2)
    pl.plot(z/1000,ion_dict["u"][2:-2,-1])
    pl.ylabel('Velocity (m/s)')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    #pl.grid('on') 
    ax6 = ax3.twiny()
    ax6.plot(z/radius,ion_dict["u"][2:-2,-1])
    ax6.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax6.set_xlim([0,z_ext[-1]/radius])

    
    ax4=pl.subplot(2,2,3)
    pl.plot(z/1000,ion_dict["P"][2:-2,-1])
    pl.ylabel('Pressure (N/m^2)')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    #pl.grid('on') 
    ax7 = ax4.twiny()
    ax7.plot(z/radius,ion_dict["P"][2:-2,-1])
    ax7.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax7.set_xlim([0,z_ext[-1]/radius])
    pl.yscale('log')
    
    ax5=pl.subplot(2,2,4)
    pl.plot(z/1000,ion_dict["T"][2:-2,-1])
    pl.ylabel('Temperature (K)')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    #pl.grid('on') 
    ax8 = ax5.twiny()
    ax8.plot(z/radius,ion_dict["T"][2:-2,-1])
    ax8.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax8.set_xlim([0,z_ext[-1]/radius])
 
    pl.suptitle(ion_dict["name"])
    pl.subplots_adjust(wspace=0.5, hspace=0.5)  
    
def input_plot(ions,electrons,neutrals,z,z_ext,A,radius):
    import matplotlib.pyplot as pl
    
    #length of arrays
    num_ion = len(ions)
    num_neut = len(neutrals)
    
    # Blue = ions
    #bl = ['blue'],['lightskyblue']
    bl = [0,0.3,0.6],[0,0.4,0.8],[0,0.5,1],[0.2,0.6,1],[0.4,0.7,1],[0.6,0.8,1],[0.8,0.9,1]
    # Green = electrons
    #gr = ['goldenrod']
    gr = [0.2,1,0.2]
    # Red = neutrals
    #rd = ['maroon'],['red'],['salmon']
    rd = [0.6,0,0],[0.8,0,0],[1,0,0],[1,0.2,0.2],[1,0.4,0.4],[1,0.6,0.6],[1,0.8,0.8]
    
    fig = pl.figure(figsize=(8,10))
    ax1 = fig.add_subplot(4,2,1)
    ax1.plot(z/1000,A[2:-2],color='black')
    ax1.set_xlim([0,z_ext[-1]/1000])
    ax1.set_xlabel([])
    ax1.xaxis.set_ticklabels([])
    ax1.set_ylabel('Cross-sectional Area \n A $(m^{2})$')
    #ax1.grid('on')
    ax1.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax2 = ax1.twiny()
    ax2.plot(z/radius,A[2:-2],color='black')
    ax2.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax2.set_xlim([0,z_ext[-1]/radius])
    ax2.tick_params(direction='in')

    
    ax3 = fig.add_subplot(4,2,2)
    ax3.plot(z/1000,ions[1]["u"][2:-2,0],color='blue') 
    ax3.plot(z/1000,ions[2]["u"][2:-2,0],color='lightskyblue')
    ax3.plot(z/1000,electrons["u"][2:-2,0],color='orange')
    ax3.set_ylabel('Initial Velocity $(ms^{-1})$')
    ax3.set_xlabel([])
    ax3.xaxis.set_ticklabels([])
    ax3.set_xlim([0,z_ext[-1]/1000])
    #ax3.grid('on')
    ax3.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax4 = ax3.twiny()
    ax4.plot(z/1000,electrons["u"][2:-2,0],color='orange')
    
    ax4.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax4.set_xlim([0,z_ext[-1]/radius]) 
    
    
    ax5 = fig.add_subplot(4,2,3)
    ax5.plot(z/1000,ions[1]["n"][2:-2,0],color='blue')    
    ax5.plot(z/1000,ions[2]["n"][2:-2,0],color='lightskyblue') 
    ax5.plot(z/1000,electrons["n"][2:-2,0],color='orange')
    #neutrals
    ax5.plot(z/1000,neutrals[1]["n"][2:-2],color='maroon')
    ax5.plot(z/1000,neutrals[2]["n"][2:-2],color='red')
    ax5.plot(z/1000,neutrals[3]["n"][2:-2],color='salmon')
    ax5.set_ylabel('Initial Number \n Density $(m^{-3})$')
    ax5.set_xlabel([])
    ax5.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax5.xaxis.set_ticklabels([])
    ax5.set_xlim([0,z_ext[-1]/1000])
    ax6=ax5.twiny()
    ax6.plot(z/radius,ions[1]["n"][2:-2,0],color='blue') 
    ax6.set_xlim([0,z_ext[-1]/radius])
    ax6.xaxis.set_ticklabels([])
    ax6.tick_params(direction='in')
    #pl.grid('on')  
    pl.yscale('log')
    
    ax7 = pl.subplot(4,2,4)
    ax7.plot(z/1000,ions[1]["rho"][2:-2,0],color='blue') 
    ax7.plot(z/1000,ions[2]["rho"][2:-2,0],color='lightskyblue')
    ax7.plot(z/1000,electrons["rho"][2:-2,0],color='orange')
    #neutrals
    ax7.plot(z/1000,neutrals[1]["rho"][2:-2],color='maroon') 
    ax7.plot(z/1000,neutrals[2]["rho"][2:-2],color='red') 
    ax7.plot(z/1000,neutrals[3]["rho"][2:-2],color='salmon') 
    ax7.set_ylabel('Initial Mass \n Density $(m^{-3})$')
    ax7.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax7.set_xlabel([])
    ax7.xaxis.set_ticklabels([])
    ax7.set_xlim([0,z_ext[-1]/1000])
    ax8 = ax7.twiny()
    ax8.plot(z/radius,ions[1]["rho"][2:-2,0],color='blue') 
    ax8.set_xlim([0,z_ext[-1]/radius])
    ax8.xaxis.set_ticklabels([])
    ax8.tick_params(direction='in')
    #pl.grid('on') 
    pl.yscale('log')
    
    ax9 = pl.subplot(4,2,5)
    #ions
    ax9.plot(z/1000,ions[1]["S"][2:-2],color='blue')
    ax9.plot(z/1000,ions[2]["S"][2:-2],color='lightskyblue')
    ax9.plot(z/1000,electrons["S"][2:-2],color='orange')    
    ax9.set_ylabel('Mass Production Rate \n $(kg m^{-3} s^{-1})$')
    ax9.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    pl.xlabel([])
    ax7.xaxis.set_ticklabels([])
    pl.xlim([0,z_ext[-1]/1000])
    ax10=ax9.twiny()
    ax10.plot(z/radius,ions[1]["S"][2:-2],color='blue')
    ax10.set_xlim([0,z_ext[-1]/radius])
    ax10.tick_params(direction='in')
    ax10.xaxis.set_ticklabels([])
    #pl.grid('on')  
    pl.yscale('log')
    
    
    ax11=pl.subplot(4,2,6)
    #ions
    ax11.plot(z/1000,ions[1]["T"][2:-2,0],color='blue')
    ax11.plot(z/1000,ions[2]["T"][2:-2,0],color='lightskyblue')
    ax11.plot(z/1000,electrons["T"][2:-2,0],color='orange')
    ax11.plot(z/1000,neutrals[1]["T"][2:-2],color='red') #all neutrals have same temp
    ax11.set_ylabel('Temperature (K)')
    ax11.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax11.set_xlabel([])
    ax11.xaxis.set_ticklabels([])
    pl.xlim([0,z_ext[-1]/1000])
    ax12=ax11.twiny()
    ax12.plot(z/radius,ions[1]["T"][2:-2,0],color='blue')
    ax12.set_xlim([0,z_ext[-1]/radius])
    ax12.tick_params(direction='in')
    ax12.xaxis.set_ticklabels([])
    #pl.grid('on') 
    
    ax13=pl.subplot(4,2,7)
    #ions
    ax13.plot(z/1000,ions[1]["P"][2:-2,0],color='blue')
    ax13.plot(z/1000,ions[2]["P"][2:-2,0],color='lightskyblue')
    ax13.plot(z/1000,electrons["P"][2:-2,0],color='orange')
    ax13.plot(z/1000,neutrals[1]["P"][2:-2],color='red') # same temp same pressure
    ax13.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    pl.ylabel('Pressure $(Nm^{-2})$')
    pl.xlabel('Distance Along Field Line (km)')
    pl.yscale('log')
    pl.xlim([0,z_ext[-1]/1000])
    ax14 =ax13.twiny()
    ax14.plot(z/1000,ions[1]["P"][2:-2,0],color='blue')
    ax14.set_xlim([0,z_ext[-1]/radius])
    ax14.xaxis.set_ticklabels([])
    ax14.tick_params(direction='in')
    #pl.grid('on')    
    
    
    ax15=pl.subplot(4,2,8)
    #ions
    ax15.plot(z/1000,ions[1]["kappa"][2:-2,0],color='blue')
    ax15.plot(z/1000,ions[2]["kappa"][2:-2,0],color='lightskyblue')
    ax15.plot(z/1000,electrons["kappa"][2:-2,0],color='orange')
    ax15.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    pl.ylabel('Heat Conductivity \n $(J m^{-1}s^{-1}K^{-1})$')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    ax16=ax15.twiny()
    ax16.plot(z/1000,ions[1]["kappa"][2:-2,0],color='blue')
    ax16.set_xlim([0,z_ext[-1]/radius])
    ax16.xaxis.set_ticklabels([])
    ax16.tick_params(direction='in')
    #pl.grid('on')   
    
    
    pl.subplots_adjust(wspace=0.5, hspace=0.0)  
    pl.suptitle("Input")  
    
def output_plots(ions,electrons,neutrals,z,z_ext,A,radius):
    import matplotlib.pyplot as pl
    
    #length of arrays
    num_ion = len(ions)
    num_neut = len(neutrals)
    
    fig = pl.figure(figsize=(8,10))
    ax1 = fig.add_subplot(4,2,1)
    ax1.plot(z/1000,A[2:-2],color='black')
    ax1.set_xlim([0,z_ext[-1]/1000])
    ax1.set_xlabel([])
    ax1.xaxis.set_ticklabels([])
    ax1.set_ylabel('Cross-sectional Area \n A $(m^{2})$')
    #ax1.grid('on')
    ax1.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax2 = ax1.twiny()
    ax2.plot(z/radius,A[2:-2],color='black')
    ax2.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax2.set_xlim([0,z_ext[-1]/radius])
    ax2.tick_params(direction='in')

    
    ax3 = fig.add_subplot(4,2,2)
    ax3.plot(z/1000,ions[1]["u"][2:-2,-1],color='blue') 
    ax3.plot(z/1000,ions[2]["u"][2:-2,-1],color='lightskyblue')
    ax3.plot(z/1000,electrons["u"][2:-2,-1],color='orange')
    ax3.set_ylabel('Initial Velocity $(ms^{-1})$')
    ax3.set_xlabel([])
    ax3.xaxis.set_ticklabels([])
    ax3.set_xlim([0,z_ext[-1]/1000])
    #ax3.grid('on')
    ax3.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax4 = ax3.twiny()
    ax4.plot(z/1000,electrons["u"][2:-2,-1],color='orange')
    
    ax4.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax4.set_xlim([0,z_ext[-1]/radius]) 
    
    
    ax5 = fig.add_subplot(4,2,3)
    ax5.plot(z/1000,ions[1]["n"][2:-2,-1],color='blue')    
    ax5.plot(z/1000,ions[2]["n"][2:-2,-1],color='lightskyblue') 
    ax5.plot(z/1000,electrons["n"][2:-2,-1],color='orange')
    #neutrals
    # ax5.plot(z/1000,neutrals[1]["n"][2:-2],color='maroon')
    # ax5.plot(z/1000,neutrals[2]["n"][2:-2],color='red')
    # ax5.plot(z/1000,neutrals[3]["n"][2:-2],color='salmon')
    ax5.set_ylabel('Initial Number \n Density $(m^{-3})$')
    ax5.set_xlabel([])
    ax5.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax5.xaxis.set_ticklabels([])
    ax5.set_xlim([0,z_ext[-1]/1000])
    ax6=ax5.twiny()
    ax6.plot(z/radius,ions[1]["n"][2:-2,-1],color='blue') 
    ax6.set_xlim([0,z_ext[-1]/radius])
    ax6.xaxis.set_ticklabels([])
    ax6.tick_params(direction='in')
    #pl.grid('on')  
    pl.yscale('log')
    
    ax7 = pl.subplot(4,2,4)
    ax7.plot(z/1000,ions[1]["rho"][2:-2,-1],color='blue') 
    ax7.plot(z/1000,ions[2]["rho"][2:-2,-1],color='lightskyblue')
    ax7.plot(z/1000,electrons["rho"][2:-2,-1],color='orange')
    #neutrals
    # ax7.plot(z/1000,neutrals[1]["rho"][2:-2],color='maroon') 
    # ax7.plot(z/1000,neutrals[2]["rho"][2:-2],color='red') 
    # ax7.plot(z/1000,neutrals[3]["rho"][2:-2],color='salmon') 
    ax7.set_ylabel('Initial Mass \n Density $(m^{-3})$')
    ax7.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax7.set_xlabel([])
    ax7.xaxis.set_ticklabels([])
    ax7.set_xlim([0,z_ext[-1]/1000])
    ax8 = ax7.twiny()
    ax8.plot(z/radius,ions[1]["rho"][2:-2,-1],color='blue') 
    ax8.set_xlim([0,z_ext[-1]/radius])
    ax8.xaxis.set_ticklabels([])
    ax8.tick_params(direction='in')
    #pl.grid('on') 
    pl.yscale('log')
    
    ax9 = pl.subplot(4,2,5)
    #ions
    ax9.plot(z/1000,ions[1]["S"][2:-2],color='blue')
    ax9.plot(z/1000,ions[2]["S"][2:-2],color='lightskyblue')
    ax9.plot(z/1000,electrons["S"][2:-2],color='orange')    
    ax9.set_ylabel('Mass Production Rate \n $(kg m^{-3} s^{-1})$')
    ax9.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    pl.xlabel([])
    ax7.xaxis.set_ticklabels([])
    pl.xlim([0,z_ext[-1]/1000])
    ax10=ax9.twiny()
    ax10.plot(z/radius,ions[1]["S"][2:-2],color='blue')
    ax10.set_xlim([0,z_ext[-1]/radius])
    ax10.tick_params(direction='in')
    ax10.xaxis.set_ticklabels([])
    #pl.grid('on')  
    pl.yscale('log')
    
    
    ax11=pl.subplot(4,2,6)
    #ions
    ax11.plot(z/1000,ions[1]["T"][2:-2,-1],color='blue')
    ax11.plot(z/1000,ions[2]["T"][2:-2,-1],color='lightskyblue')
    ax11.plot(z/1000,electrons["T"][2:-2,-1],color='orange')
    #ax11.plot(z/1000,neutrals[1]["T"][2:-2],color='red') #all neutrals have same temp
    ax11.set_ylabel('Temperature (K)')
    ax11.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax11.set_xlabel([])
    ax11.xaxis.set_ticklabels([])
    pl.xlim([0,z_ext[-1]/1000])
    ax12=ax11.twiny()
    ax12.plot(z/radius,ions[1]["T"][2:-2,-1],color='blue')
    ax12.set_xlim([0,z_ext[-1]/radius])
    ax12.tick_params(direction='in')
    ax12.xaxis.set_ticklabels([])
    #pl.grid('on') 
    
    ax13=pl.subplot(4,2,7)
    #ions
    ax13.plot(z/1000,ions[1]["P"][2:-2,-1],color='blue')
    ax13.plot(z/1000,ions[2]["P"][2:-2,-1],color='lightskyblue')
    ax13.plot(z/1000,electrons["P"][2:-2,-1],color='orange')
    #ax13.plot(z/1000,neutrals[1]["P"][2:-2],color='red') # same temp same pressure
    ax13.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    pl.ylabel('Pressure $(Nm^{-2})$')
    pl.xlabel('Distance Along Field Line (km)')
    pl.yscale('log')
    pl.xlim([0,z_ext[-1]/1000])
    ax14 =ax13.twiny()
    ax14.plot(z/1000,ions[1]["P"][2:-2,-1],color='blue')
    ax14.set_xlim([0,z_ext[-1]/radius])
    ax14.xaxis.set_ticklabels([])
    ax14.tick_params(direction='in')
    #pl.grid('on')    
    
    
    ax15=pl.subplot(4,2,8)
    #ions
    ax15.plot(z/1000,ions[1]["kappa"][2:-2,-1],color='blue')
    ax15.plot(z/1000,ions[2]["kappa"][2:-2,-1],color='lightskyblue')
    ax15.plot(z/1000,electrons["kappa"][2:-2,-1],color='orange')
    ax15.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    pl.ylabel('Heat Conductivity \n $(J m^{-1}s^{-1}K^{-1})$')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    ax16=ax15.twiny()
    ax16.plot(z/1000,ions[1]["kappa"][2:-2,-1],color='blue')
    ax16.set_xlim([0,z_ext[-1]/radius])
    ax16.xaxis.set_ticklabels([])
    ax16.tick_params(direction='in')
    #pl.grid('on')   
    
    
    pl.subplots_adjust(wspace=0.5, hspace=0.0)  
    pl.suptitle("Output")
    
def plot_single(z,z_ext,item,radius,item_str):
    import matplotlib.pyplot as pl
    ax1 = pl.subplot(1,1,1)
    ax1.plot(z/1000,item,color=[0,0,0])
    ax1.set_xlim([0,z_ext[-1]/1000])
    ax1.set_xlabel('Distance Along Field Line (km)')
    ax1.set_ylabel(item_str)
    #ax1.grid('on')
    ax2 = ax1.twiny()
    ax2.plot(z/radius,item,color=[0,0,0])
    ax2.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax2.set_xlim([0,z_ext[-1]/radius])
