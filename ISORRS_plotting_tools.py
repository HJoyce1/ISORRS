# Module for plotting functions to tidy code up
# Written by Carley Martin


# quick plotting function just to quickly check on things
def plot_me_quick(x,y,xlabel,ylabel):
    import matplotlib.pyplot as pl
    # plot x vs y
    pl.plot(x,y)
    # add labels
    pl.ylabel(ylabel)
    pl.xlabel(xlabel)

# quick plotting function with log y-scale
def plot_me_log(x,y,xlabel,ylabel):
    import matplotlib.pyplot as pl
    # plot x vs y
    pl.plot(x,y)
    # set labels
    pl.ylabel(ylabel)
    pl.xlabel(xlabel)
    # set y-scale as log
    pl.yscale('log')

# results plot showing the ambipolar electric field, acceleration terms, ion flux, electron flux and total flux
def results_plot(z,z_ext,radius,num_ion,e_charge,E,ion_dict,electron_dict,ac,ag,e_flux,ion_flux,ion_flux_tot):
    import matplotlib.pyplot as pl
    # set figure size
    pl.figure(figsize=(8,10))
    # first plot, ambipolar electric field
    ax1 = pl.subplot(5,1,1)
    # plot grid vs E field
    pl.plot(z/1000,(E)/1e-7,linestyle='-',color='black')
    # label y axis
    pl.ylabel('a) $E_{\parallel}$ $(V$ $m^{-1}) $x$ 10^{-7}$')
    # set limit on x-axis based on grid
    pl.xlim([0,z_ext[-1]/1000])
    # set no label or tick numbers on first axes as want to twin this axis
    ax1.set_xlabel('')
    ax1.xaxis.set_ticklabels([])
    # aet ticks
    ax1.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    # new axis on top of above axis
    ax2 = ax1.twiny()
    # replot electric field but with jupiter radii
    ax2.plot(z/radius,(E)/1e-7,linestyle= '-',color='black')
    # add x axis label to top of whole figure
    ax2.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    # set limit based on jupiter radii
    ax2.set_xlim([0,z_ext[-1]/radius])
    # set ticks
    ax2.tick_params(direction='in')

    # second plot, acceleration terms
    ax3 = pl.subplot(5,1,2)
    # plot H+ electric field vs grid
    pl.plot(z/1000,((e_charge/ion_dict[1]["mass"])*(E)),linestyle= '-',color='blue')
    # plot H3+ electric field vs grid
    pl.plot(z/1000,((e_charge/ion_dict[2]["mass"])*(E)),linestyle= '-',color='lightskyblue')
    # plot centrifugal acceleration vs grid
    pl.plot(z/1000,ac,linestyle= '--',color='peru')
    # plot graviatational acceleration vs grid
    pl.plot(z/1000,ag,linestyle= '-.',color='saddlebrown')
    # add y label
    pl.ylabel('b) Acceleration Terms \n $(m$ $s^{-2})$')
    # set limit based on grid
    pl.xlim([0,z_ext[-1]/1000])
    # set ticks
    ax3.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    # new axis on top of above axis
    ax4 = ax3.twiny()
    # replot H+ electric field vs jupiter radii this time
    ax4.plot(z/radius,(e_charge/ion_dict[1]["mass"])*(E),linestyle='-',color='blue')
    # set limit based on jupiter radii
    ax4.set_xlim([0,z_ext[-1]/radius])
    # remove label and numbers from both axes
    ax3.set_xlabel('')
    ax3.xaxis.set_ticklabels([])
    ax4.set_xlabel('')
    ax4.xaxis.set_ticklabels([])
    # set ticks for second axes
    ax4.tick_params(direction='in')
    # add legend to visualise different lines
    ax3.legend(['$H^+$ Electric Field','$H_3^+$ Electric Field','Centrifugal', 'Gravitational'],loc='upper right',bbox_to_anchor=(0.99,0.98))

    # third plot, electron flux
    ax5 = pl.subplot(5,1,3)
    # plot on log scale
    pl.yscale('log')
    # plot electron flux vs grid
    pl.plot(z/1000,e_flux[2:-2,-1],color='orange')
    # label y axis
    pl.ylabel('c) Electron Flux \n $(m^{-2}$ $s^{-1})$ x A $(m^2)$')
    # set limit based on grid
    pl.xlim([0,z_ext[-1]/1000])
    # set ticks
    ax5.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    # new axis on top of above axis
    ax6 = ax5.twiny()
    # replot electron flux but vs jupiter radii now
    ax6.plot(z/radius,e_flux[2:-2,-1],color='orange')
    # set limits based on RJ
    ax6.set_xlim([0,z_ext[-1]/radius])
    # remove labels and numbers from both axes
    ax5.set_xlabel('')
    ax5.xaxis.set_ticklabels([])
    ax6.set_xlabel('')
    ax6.xaxis.set_ticklabels([])
    # set ticks for second axes
    ax6.tick_params(direction='in')
    # set legend
    ax5.legend(['$e^-$'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # fourth plot, ion fluxes
    ax7 = pl.subplot(5,1,4)
    # plot using log scale
    pl.yscale('log')
    # plot H+ flux vs grid
    pl.plot(z/1000,ion_flux[2:-2,-1,0],color='blue')
    # plot H3+ flux vs grid
    pl.plot(z/1000,ion_flux[2:-2,-1,1],color='lightskyblue')
    # add label to y-scale
    pl.ylabel('d) Ion Flux \n $(m^{-2}$ $s^{-1})$ x A $(m^2)$')
    # set limit based on grid
    pl.xlim([0,z_ext[-1]/1000])
    # set ticks
    ax7.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    # set new axis on top of above axis
    ax8 = ax7.twiny()
    # replot H+ flux vs now jupiter radii
    ax8.plot(z/radius,ion_flux[2:-2,-1,0],color='blue')
    # set limits based on jupiter radii
    ax8.set_xlim([0,z_ext[-1]/radius])
    # remove labels and numbers for these axes
    ax8.set_xlabel('')
    ax7.set_xlabel('')
    ax8.xaxis.set_ticklabels([])
    ax7.xaxis.set_ticklabels([])
    # set ticks for top axes
    ax8.tick_params(direction='in')
    # set legend
    ax7.legend(['$H^+$','$H_3^+$'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # fifth plot, total ion flux
    ax9 = pl.subplot(5,1,5)
    # plot using log scale
    pl.yscale('log')
    # plot electron flux vs grid
    pl.plot(z/1000,e_flux[2:-2,-1],color='orange')
    # plot H+ flux vs grid
    pl.plot(z/1000,ion_flux[2:-2,-1,0],color='blue')
    # plot H3+ flux vs grid
    pl.plot(z/1000,ion_flux[2:-2,-1,1],color='lightskyblue')
    # plot H+ and H3+ combined flux vs grid
    pl.plot(z/1000,ion_flux_tot[2:-2,-1],color='black')
    # add y label
    pl.ylabel('e) Total Flux \n $(m^{-2}$ $s^{-1})$ x A $(m^2)$')
    # set limit based on grid
    pl.xlim([0,z_ext[-1]/1000])
    # set ticks
    ax9.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    # add axes on top of current axes
    ax10 = ax9.twiny()
    # replot combined flux vs jupiter radii now
    ax10.plot(z/1000,ion_flux_tot[2:-2,-1],color='black')
    # set limit based on RJ
    ax10.set_xlim([0,z_ext[-1]/radius])
    # remove label and numbers on RJ axes
    ax10.set_xlabel('')
    ax10.xaxis.set_ticklabels([])
    # set ticks on top axes
    ax10.tick_params(direction='in')
    # set x label based on grid
    ax9.set_xlabel('Distance Along Field Line (km)')
    # set legend
    ax9.legend(['$e^-$','$H^+$','$H_3^+$','$H^+$ & $H_3^+$'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # add title to plots and adjust so plots stack
    pl.subplots_adjust(hspace=0.0,wspace=0.5)
    pl.suptitle("Results")

# individual plots looking at specific species
def species_plot(z,z_ext,ion_dict,radius):
    import matplotlib.pyplot as pl
    # create figure
    pl.figure(figsize=(15,8))
    # number density plot
    ax1 = pl.subplot(2,2,1)
    # log scale
    pl.yscale('log')
    # grid vs num density
    pl.plot(z/1000,ion_dict["n"][2:-2,-1])
    pl.ylabel('Density (kg/m^3)')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    # redo plot for jupiter radii
    ax2 = ax1.twiny()
    ax2.plot(z/radius,ion_dict["n"][2:-2,-1])
    ax2.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax2.set_xlim([0,z_ext[-1]/radius])

    # plot for velocity
    ax3 =pl.subplot(2,2,2)
    # plot vs grid
    pl.plot(z/1000,ion_dict["u"][2:-2,-1])
    pl.ylabel('Velocity (m/s)')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    # plot vs jupiter radii
    ax6 = ax3.twiny()
    ax6.plot(z/radius,ion_dict["u"][2:-2,-1])
    ax6.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax6.set_xlim([0,z_ext[-1]/radius])

    # plot for pressure
    ax4=pl.subplot(2,2,3)
    # plot vs grid
    pl.plot(z/1000,ion_dict["P"][2:-2,-1])
    pl.ylabel('Pressure (N/m^2)')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    # plot vs RJ, tinned axes
    ax7 = ax4.twiny()
    ax7.plot(z/radius,ion_dict["P"][2:-2,-1])
    ax7.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax7.set_xlim([0,z_ext[-1]/radius])
    pl.yscale('log')

    # plot for temperature
    ax5=pl.subplot(2,2,4)
    # plot vs grid
    pl.plot(z/1000,ion_dict["T"][2:-2,-1])
    pl.ylabel('Temperature (K)')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    # plot vs radii, twinned axes
    ax8 = ax5.twiny()
    ax8.plot(z/radius,ion_dict["T"][2:-2,-1])
    ax8.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax8.set_xlim([0,z_ext[-1]/radius])

    # add title and adjust for stacking
    pl.suptitle(ion_dict["name"])
    pl.subplots_adjust(wspace=0.5, hspace=0.5)


# plots of the input data
def input_plot(ions,electrons,neutrals,z,z_ext,A,radius):
    import matplotlib.pyplot as pl
    # create figure
    fig = pl.figure(figsize=(8,10))

    #plot for cross sectional area
    ax1 = fig.add_subplot(4,2,1)
    # plot vs grid
    ax1.plot(z/1000,A[2:-2],color='black')
    ax1.set_xlim([0,z_ext[-1]/1000])
    ax1.set_xlabel([])
    ax1.xaxis.set_ticklabels([])
    ax1.set_ylabel('Cross-sectional Area \n A $(m^{2})$')
    ax1.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    # plot vs jupiter radii, twinned axes
    ax2 = ax1.twiny()
    ax2.plot(z/radius,A[2:-2],color='black')
    ax2.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax2.set_xlim([0,z_ext[-1]/radius])
    ax2.tick_params(direction='in')

    # plot for velocity
    ax3 = fig.add_subplot(4,2,2)
    # H+ velocity vs grid
    ax3.plot(z/1000,ions[1]["u"][2:-2,0],color='blue')
    # H3+ velocity vs grid
    ax3.plot(z/1000,ions[2]["u"][2:-2,0],color='lightskyblue')
    # electron velocity vs grid
    ax3.plot(z/1000,electrons["u"][2:-2,0],color='orange')
    ax3.set_ylabel('Initial Velocity $(ms^{-1})$')
    ax3.set_xlabel([])
    ax3.xaxis.set_ticklabels([])
    ax3.set_xlim([0,z_ext[-1]/1000])
    ax3.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    # twin axes for jupiter radii
    ax4 = ax3.twiny()
    # replot electron velocity so twinned axes 'have data'
    ax4.plot(z/1000,electrons["u"][2:-2,0],color='orange')
    ax4.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax4.set_xlim([0,z_ext[-1]/radius])
    # add legend
    ax3.legend(['$H^+$','$H_3^+$','$e^-$'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # plot for number density
    ax5 = fig.add_subplot(4,2,3)
    # H+ vs grid
    ax5.plot(z/1000,ions[1]["n"][2:-2,0],color='blue')
    # H3+ vs grid
    ax5.plot(z/1000,ions[2]["n"][2:-2,0],color='lightskyblue')
    # electrons vs grid
    ax5.plot(z/1000,electrons["n"][2:-2,0],color='orange')
    # neutrals
    # H2 vs grid
    ax5.plot(z/1000,neutrals[1]["n"][2:-2],color='maroon')
    # He vs grid
    ax5.plot(z/1000,neutrals[2]["n"][2:-2],color='red')
    # H vs grid
    ax5.plot(z/1000,neutrals[3]["n"][2:-2],color='salmon')
    ax5.set_ylabel('Initial Number \n Density $(m^{-3})$')
    ax5.set_xlabel([])
    ax5.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax5.xaxis.set_ticklabels([])
    ax5.set_xlim([0,z_ext[-1]/1000])
    # twinned axes for RJ
    ax6=ax5.twiny()
    # plot for twinned axes
    ax6.plot(z/radius,ions[1]["n"][2:-2,0],color='blue')
    ax6.set_xlim([0,z_ext[-1]/radius])
    ax6.xaxis.set_ticklabels([])
    ax6.tick_params(direction='in')
    pl.yscale('log')
    # add legend
    ax5.legend(['$H^+$','$H_3^+$','$e^-$','$H_2$','He','H'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # plot for mass density
    ax7 = pl.subplot(4,2,4)
    # H+ vs grid
    ax7.plot(z/1000,ions[1]["rho"][2:-2,0],color='blue')
    # H3+ vs grid
    ax7.plot(z/1000,ions[2]["rho"][2:-2,0],color='lightskyblue')
    # electrons vs grid
    ax7.plot(z/1000,electrons["rho"][2:-2,0],color='orange')
    # neutrals
    # H2 vs grid
    ax7.plot(z/1000,neutrals[1]["rho"][2:-2],color='maroon')
    # He vs grid
    ax7.plot(z/1000,neutrals[2]["rho"][2:-2],color='red')
    # H vs grid
    ax7.plot(z/1000,neutrals[3]["rho"][2:-2],color='salmon')
    ax7.set_ylabel('Initial Mass \n Density $(m^{-3})$')
    ax7.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax7.set_xlabel([])
    ax7.xaxis.set_ticklabels([])
    ax7.set_xlim([0,z_ext[-1]/1000])
    # twin axes so can show RJ
    ax8 = ax7.twiny()
    ax8.plot(z/radius,ions[1]["rho"][2:-2,0],color='blue')
    ax8.set_xlim([0,z_ext[-1]/radius])
    ax8.xaxis.set_ticklabels([])
    ax8.tick_params(direction='in')
    pl.yscale('log')
    # legend
    ax7.legend(['$H^+$','$H_3^+$','$e^-$','$H_2$','He','H'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # plot for mass prodeuction rate
    ax9 = pl.subplot(4,2,5)
    # H+ vs grid
    ax9.plot(z/1000,ions[1]["S"][2:-2],color='blue')
    # H3+ vs grid
    ax9.plot(z/1000,ions[2]["S"][2:-2],color='lightskyblue')
    # electrons vs grid
    ax9.plot(z/1000,electrons["S"][2:-2],color='orange')
    ax9.set_ylabel('Mass Production Rate \n $(kg m^{-3} s^{-1})$')
    ax9.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    pl.xlabel([])
    ax7.xaxis.set_ticklabels([])
    pl.xlim([0,z_ext[-1]/1000])
    # twin axes so can show RJ
    ax10=ax9.twiny()
    ax10.plot(z/radius,ions[1]["S"][2:-2],color='blue')
    ax10.set_xlim([0,z_ext[-1]/radius])
    ax10.tick_params(direction='in')
    ax10.xaxis.set_ticklabels([])
    pl.yscale('log')
    # legend
    ax9.legend(['$H^+$','$H_3^+$','$e^-$'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # plot for temperature
    ax11=pl.subplot(4,2,6)
    # H+ vs grid
    ax11.plot(z/1000,ions[1]["T"][2:-2,0],color='blue')
    # H3+ vs grid
    ax11.plot(z/1000,ions[2]["T"][2:-2,0],color='lightskyblue')
    # electrons vs grid
    ax11.plot(z/1000,electrons["T"][2:-2,0],color='orange')
    # neutrals vs grid - all neutrals have same temp, using H2 here
    ax11.plot(z/1000,neutrals[1]["T"][2:-2],color='red')
    ax11.set_ylabel('Temperature (K)')
    ax11.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax11.set_xlabel([])
    ax11.xaxis.set_ticklabels([])
    pl.xlim([0,z_ext[-1]/1000])
    # twin axes for RJ
    ax12=ax11.twiny()
    ax12.plot(z/radius,ions[1]["T"][2:-2,0],color='blue')
    ax12.set_xlim([0,z_ext[-1]/radius])
    ax12.tick_params(direction='in')
    ax12.xaxis.set_ticklabels([])
    # legend
    ax11.legend(['$H^+$','$H_3^+$','$e^-$','Neutrals'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # plot for pressure
    ax13=pl.subplot(4,2,7)
    # H+ vs grid
    ax13.plot(z/1000,ions[1]["P"][2:-2,0],color='blue')
    # H3+ vs grid
    ax13.plot(z/1000,ions[2]["P"][2:-2,0],color='lightskyblue')
    # electrons vs grid
    ax13.plot(z/1000,electrons["P"][2:-2,0],color='orange')
    # neutrals vs grid - same temp same pressure
    ax13.plot(z/1000,neutrals[1]["P"][2:-2],color='red')
    ax13.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    pl.ylabel('Pressure $(Nm^{-2})$')
    pl.xlabel('Distance Along Field Line (km)')
    pl.yscale('log')
    pl.xlim([0,z_ext[-1]/1000])
    # twin axes for RJ
    ax14 =ax13.twiny()
    ax14.plot(z/1000,ions[1]["P"][2:-2,0],color='blue')
    ax14.set_xlim([0,z_ext[-1]/radius])
    ax14.xaxis.set_ticklabels([])
    ax14.tick_params(direction='in')
    # legend
    ax13.legend(['$H^+$','$H_3^+$','$e^-$','Neutrals'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # plot for heat conductivity
    ax15=pl.subplot(4,2,8)
    # H+ vs grid
    ax15.plot(z/1000,ions[1]["kappa"][2:-2,0],color='blue')
    # H3+ vs grid
    ax15.plot(z/1000,ions[2]["kappa"][2:-2,0],color='lightskyblue')
    # electrons vs grid
    ax15.plot(z/1000,electrons["kappa"][2:-2,0],color='orange')
    ax15.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    pl.ylabel('Heat Conductivity \n $(J m^{-1}s^{-1}K^{-1})$')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    # twin axes for RJ
    ax16=ax15.twiny()
    ax16.plot(z/1000,ions[1]["kappa"][2:-2,0],color='blue')
    ax16.set_xlim([0,z_ext[-1]/radius])
    ax16.xaxis.set_ticklabels([])
    ax16.tick_params(direction='in')
    # legend
    ax15.legend(['$H^+$','$H_3^+$','$e^-$'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # add title and adjust for stacking
    pl.subplots_adjust(wspace=0.5, hspace=0.0)
    pl.suptitle("Input")

# plots of the input variables after they have been interated through (outputs)
def output_plots(ions,electrons,neutrals,z,z_ext,A,radius):
    import matplotlib.pyplot as pl

    # plot for cross sectional area
    fig = pl.figure(figsize=(8,10))
    ax1 = fig.add_subplot(4,2,1)
    # plot vs grid
    ax1.plot(z/1000,A[2:-2],color='black')
    ax1.set_xlim([0,z_ext[-1]/1000])
    ax1.set_xlabel([])
    ax1.xaxis.set_ticklabels([])
    ax1.set_ylabel('Cross-sectional Area \n A $(m^{2})$')
    ax1.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    # twinned axes for RJ
    ax2 = ax1.twiny()
    ax2.plot(z/radius,A[2:-2],color='black')
    ax2.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax2.set_xlim([0,z_ext[-1]/radius])
    ax2.tick_params(direction='in')

    # plot for velocity
    ax3 = fig.add_subplot(4,2,2)
    # H+ vs grid
    ax3.plot(z/1000,ions[1]["u"][2:-2,-1],color='blue')
    # H3+ vs grid
    ax3.plot(z/1000,ions[2]["u"][2:-2,-1],color='lightskyblue')
    # electrons vs grid
    ax3.plot(z/1000,electrons["u"][2:-2,-1],color='orange')
    ax3.set_ylabel('Initial Velocity $(ms^{-1})$')
    ax3.set_xlabel([])
    ax3.xaxis.set_ticklabels([])
    ax3.set_xlim([0,z_ext[-1]/1000])
    ax3.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    # twinned axes for RJ
    ax4 = ax3.twiny()
    ax4.plot(z/1000,electrons["u"][2:-2,-1],color='orange')
    ax4.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    ax4.set_xlim([0,z_ext[-1]/radius])
    # legend
    ax3.legend(['$H^+$','$H_3^+$','$e^-$'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # plot for number density
    ax5 = fig.add_subplot(4,2,3)
    # H+ vs grid
    ax5.plot(z/1000,ions[1]["n"][2:-2,-1],color='blue')
    # H3+ vs grid
    ax5.plot(z/1000,ions[2]["n"][2:-2,-1],color='lightskyblue')
    # electrons vs grid
    ax5.plot(z/1000,electrons["n"][2:-2,-1],color='orange')
    # neutrals
    # H2 vs grid
    ax5.plot(z/1000,neutrals[1]["n"][2:-2],color='maroon')
    # He va grid
    ax5.plot(z/1000,neutrals[2]["n"][2:-2],color='red')
    # H vs grid
    ax5.plot(z/1000,neutrals[3]["n"][2:-2],color='salmon')
    ax5.set_ylabel('Initial Number \n Density $(m^{-3})$')
    ax5.set_xlabel([])
    ax5.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax5.xaxis.set_ticklabels([])
    ax5.set_xlim([0,z_ext[-1]/1000])
    # twinned axes for RJ
    ax6=ax5.twiny()
    ax6.plot(z/radius,ions[1]["n"][2:-2,-1],color='blue')
    ax6.set_xlim([0,z_ext[-1]/radius])
    ax6.xaxis.set_ticklabels([])
    ax6.tick_params(direction='in')
    pl.yscale('log')
    # legend
    ax5.legend(['$H^+$','$H_3^+$','$e^-$','$H_2$','He','H'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # plot for mass density
    ax7 = pl.subplot(4,2,4)
    # H+ vs grid
    ax7.plot(z/1000,ions[1]["rho"][2:-2,-1],color='blue')
    # H3+ vs grid
    ax7.plot(z/1000,ions[2]["rho"][2:-2,-1],color='lightskyblue')
    # electrons vs grid
    ax7.plot(z/1000,electrons["rho"][2:-2,-1],color='orange')
    # neutrals
    # H2 vs grid
    ax7.plot(z/1000,neutrals[1]["rho"][2:-2],color='maroon')
    # He vs grid
    ax7.plot(z/1000,neutrals[2]["rho"][2:-2],color='red')
    # H vs grid
    ax7.plot(z/1000,neutrals[3]["rho"][2:-2],color='salmon')
    ax7.set_ylabel('Initial Mass \n Density $(m^{-3})$')
    ax7.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax7.set_xlabel([])
    ax7.xaxis.set_ticklabels([])
    ax7.set_xlim([0,z_ext[-1]/1000])
    # twinned axes for RJ
    ax8 = ax7.twiny()
    ax8.plot(z/radius,ions[1]["rho"][2:-2,-1],color='blue')
    ax8.set_xlim([0,z_ext[-1]/radius])
    ax8.xaxis.set_ticklabels([])
    ax8.tick_params(direction='in')
    pl.yscale('log')
    # legend
    ax7.legend(['$H^+$','$H_3^+$','$e^-$','$H_2$','He','H'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # plot for mass production rate
    ax9 = pl.subplot(4,2,5)
    # H+ vs grid
    ax9.plot(z/1000,ions[1]["S"][2:-2],color='blue')
    # H3+ vs grid
    ax9.plot(z/1000,ions[2]["S"][2:-2],color='lightskyblue')
    # electrons vs grid
    ax9.plot(z/1000,electrons["S"][2:-2],color='orange')
    ax9.set_ylabel('Mass Production Rate \n $(kg m^{-3} s^{-1})$')
    ax9.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    pl.xlabel([])
    ax7.xaxis.set_ticklabels([])
    pl.xlim([0,z_ext[-1]/1000])
    # twinned axes for RJ
    ax10=ax9.twiny()
    ax10.plot(z/radius,ions[1]["S"][2:-2],color='blue')
    ax10.set_xlim([0,z_ext[-1]/radius])
    ax10.tick_params(direction='in')
    ax10.xaxis.set_ticklabels([])
    pl.yscale('log')
    # legend
    ax9.legend(['$H^+$','$H_3^+$','$e^-$'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # plot for temperature
    ax11=pl.subplot(4,2,6)
    # H+ vs grid
    ax11.plot(z/1000,ions[1]["T"][2:-2,-1],color='blue')
    # H3+ vs grid
    ax11.plot(z/1000,ions[2]["T"][2:-2,-1],color='lightskyblue')
    # electrons vs grid
    ax11.plot(z/1000,electrons["T"][2:-2,-1],color='orange')
    # neutrals vs grid - all neutrals have same temp
    ax11.plot(z/1000,neutrals[1]["T"][2:-2],color='red')
    ax11.set_ylabel('Temperature (K)')
    ax11.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    ax11.set_xlabel([])
    ax11.xaxis.set_ticklabels([])
    pl.xlim([0,z_ext[-1]/1000])
    # twinned axes to show RJ
    ax12=ax11.twiny()
    ax12.plot(z/radius,ions[1]["T"][2:-2,-1],color='blue')
    ax12.set_xlim([0,z_ext[-1]/radius])
    ax12.tick_params(direction='in')
    ax12.xaxis.set_ticklabels([])
    # legend
    ax11.legend(['$H^+$','$H_3^+$','$e^-$','$H_2$','He','H'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # plot for pressure
    ax13=pl.subplot(4,2,7)
    # H+ vs grid
    ax13.plot(z/1000,ions[1]["P"][2:-2,-1],color='blue')
    # H3+ vs grid
    ax13.plot(z/1000,ions[2]["P"][2:-2,-1],color='lightskyblue')
    # electrons vs grid
    ax13.plot(z/1000,electrons["P"][2:-2,-1],color='orange')
    # neutrals vs grid - same temp = same pressure
    ax13.plot(z/1000,neutrals[1]["P"][2:-2],color='red')
    ax13.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    pl.ylabel('Pressure $(Nm^{-2})$')
    pl.xlabel('Distance Along Field Line (km)')
    pl.yscale('log')
    pl.xlim([0,z_ext[-1]/1000])
    # twinned axes to show RJ
    ax14 =ax13.twiny()
    ax14.plot(z/1000,ions[1]["P"][2:-2,-1],color='blue')
    ax14.set_xlim([0,z_ext[-1]/radius])
    ax14.xaxis.set_ticklabels([])
    ax14.tick_params(direction='in')
    # legend
    ax13.legend(['$H^+$','$H_3^+$','$e^-$','$H_2$','He','H'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # plot for heat conductvity
    ax15=pl.subplot(4,2,8)
    # H+ vs grid
    ax15.plot(z/1000,ions[1]["kappa"][2:-2,-1],color='blue')
    # H3+ vs grid
    ax15.plot(z/1000,ions[2]["kappa"][2:-2,-1],color='lightskyblue')
    # electrons vs grid
    ax15.plot(z/1000,electrons["kappa"][2:-2,-1],color='orange')
    ax15.tick_params(direction='in',bottom=True, top=False, left=True, right=True)
    pl.ylabel('Heat Conductivity \n $(J m^{-1}s^{-1}K^{-1})$')
    pl.xlabel('Distance Along Field Line (km)')
    pl.xlim([0,z_ext[-1]/1000])
    # twinned axes to show jupiter radii
    ax16=ax15.twiny()
    ax16.plot(z/1000,ions[1]["kappa"][2:-2,-1],color='blue')
    ax16.set_xlim([0,z_ext[-1]/radius])
    ax16.xaxis.set_ticklabels([])
    ax16.tick_params(direction='in')
    # legend
    ax15.legend(['$H^+$','$H_3^+$','$e^-$'],loc='upper right',bbox_to_anchor=(0.99,0.95))

    # add title and adjust for plot stacking
    pl.subplots_adjust(wspace=0.5, hspace=0.0)
    pl.suptitle("Output")

# function template for plotting
def plot_single(z,z_ext,item,radius,item_str):
    import matplotlib.pyplot as pl
    # set subplot
    ax1 = pl.subplot(1,1,1)
    # plot vs grid
    ax1.plot(z/1000,item,color='black')
    # set limit based on grid
    ax1.set_xlim([0,z_ext[-1]/1000])
    # set labels
    ax1.set_xlabel('Distance Along Field Line (km)')
    ax1.set_ylabel(item_str)
    # axes for RJ set on top of current axis
    ax2 = ax1.twiny()
    # plot vs RJ
    ax2.plot(z/radius,item,color='black')
    # set x label
    ax2.set_xlabel('Distance Along Field Line \n (Planetary Radii)')
    # set limit based on RJ
    ax2.set_xlim([0,z_ext[-1]/radius])
