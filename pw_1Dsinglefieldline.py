# =============================================================================
#             now fully generalised for any number of ionic or neutral species
# AUTHOR: C.J.Martin
# INSTITUTION: Lancaster University

# INPUTS: dipolefield.py, pw.py, pw_plotting_tools.py and planet.py modules required
#         planet_name = string with planet name (saturn or jupiter)
#         dt = time step (0.01)
#         its = iterations (100)
#         rs = rough estimate of how far in planetary radii you want to go (add a few more fo jupiter)
#         lshell = lshell of field line (33 is around 80 deg lat for both)
#         FAC_flag = flag for field aligned currents (0=no change, 1=removed)
#         CF_flag = flag for centrifugal force (0=no change, 1=removed)
#         plots = do you want plotting (1=yes, 0=no)
#         saves = do you want to save the output to file (1=yes, 0=no)

# OUTPUT: Plots for now

# NOTES: 
# All in SI units unless specified
# Part of the lancaster polar wind model
# 2019

#==============================================================================

#EDITED by Hannah Joyce, Lancaster University, 2022

# =============================================================================
#EXAMPLE: polar_wind_outflow('jupiter',0.01,10000,3,16,0,0,1,0,'anything'):
def polar_wind_outflow(planet_name,dt,its,rs,lshell,FAC_flag,CF_flag,plots,saves,run_name): 
    # ---------------------------------Imports-------------------------------------
    import numpy as np
    import matplotlib.pyplot as pl
    import dipolefield
    import pw
    import planet_edit as planet
    import pw_plotting_tools_editt as pwpl # already been incorporated?
    # ---------------------------------Start Main-------------------------------------
    # planet_name = 'jupiter'
    # dt = 0.01
    # its=10000
    # rs = 3
    # lshell=16 #y_eq[h-1]
    # FAC_flag =0
    # CF_flag =0
    # plots=1
    # saves=0
    # run_name='testB'
    #    Main folder to save plots and data to (obv change this)
    #folder = 'H:/pw_model/test_runs/more_tests_feb/'
    folder = 'C:/Users/joyceh1/OneDrive - Lancaster University/pw_model/test_runs/desktop/'
    #folder = 'C:/Users/Knowhow/OneDrive - Lancaster University/pw_model/test_runs/laptop/'
    # ---------------------------------Grid set up-------------------------------------    
    #800 spatial grids per Rs roughly - a bit more to round the numbers
    inner = 1400000  #lower boundary (1400km)
    numpoints = int(rs * 800) 
    dz = 75000.0 # grid spacing, 75km
    
    #set up grid
    z=np.zeros([numpoints,])
    z_ext=np.zeros([numpoints+4,]) #for ghost points
    z[0] = inner #set lower boundary
    z_ext[0] = inner - 2*dz # set ghost points
    z_ext[1] = inner - dz
    for k in range(1,numpoints):
        z[k] = z[k-1] + dz #fill z array with values 
    for l in range(1,numpoints+3):    
        z_ext[l+1] = z_ext[l] + dz #fill z_ext array with values
    x = np.linspace(0,numpoints-1,numpoints) # linear array 0-end for making functions
    # global zz
    # zz = z/1000
    
    
    #ghost points only
    gb_z = np.linspace(z_ext[0],z_ext[1],2)
    ge_z = np.linspace(z_ext[-2],z_ext[-1],2)
    gb_x = np.linspace(-1,0,2)
    ge_x = np.linspace(numpoints,numpoints+1,2)
    
    ghosts = [gb_z,ge_z,gb_x,ge_x] #combine to read into input module
    
    # -------------------------Physical Constants----------------------------------
    # Constants set-up
    e_charge = 1.60217662*10**-19  #C - electron charge
    m_p = 1.6726*10**-27 # kg - mass of a proton
    m_e = 9.10938356 * 10**-31 # kg - mass of an electron
    k_b = 1.38064852*10**-23 #m2 kg s-2 K-1  - boltzmann constant
    gamma = 5/3 #specific heat ratio
    
    phys_consts = [m_e,m_p,k_b,e_charge,gamma] #combine so can be used in planet module
    
    # -----------------------READ IN INITIAL CONDITIONS----------------------------
    #choose planet for input module - name of planet no capitals
    if planet_name == 'jupiter':
        consts,A,FAC,ions,electrons,neutrals = planet.jupiter(its,phys_consts,z,z_ext,x,ghosts)
        print('-------Running for Jupiter--------')
    elif planet_name == 'saturn':
        consts,A,FAC,ions,electrons,neutrals = planet.saturn(its,phys_consts,z,z_ext,x,ghosts)
        print('-------Running for Saturn--------')       
    else:
        print('Incorrect planet_name input, please use ''jupiter'' or ''saturn''')
        exit()
        
    #optional field aligned currents    
    if FAC_flag != 0:
        FAC = np.zeros(np.size(z_ext)) #if testing without field aligned currents
        print('Field aligned currents removed')
    else:
        print('Field aligned currents included')
        
    # preallocate arrays not used in initial conditions - electric field and flux
    E = np.empty([len(z)+4,its])
    E[:,0] = np.nan # no initial values for electric field
    e_flux = np.empty([len(z)+4,its])
    e_flux[:,0] = np.nan # no initial values for electron flux
    ion1_flux = np.empty([len(z)+4,its])
    ion1_flux[:,0] = np.nan # no initial values for ion1 flux
    ion2_flux = np.empty([len(z)+4,its])
    ion2_flux[:,0] = np.nan # no initial values for ion2 flux
    
    # unpack constants
    radius = consts[0]
    mass_planet= consts[1]
    b0= consts[2]
    rot_period= consts[3]
    dipole_offset= consts[4]
    #    g= consts[5]
    #    lshell=33 #equivalent to 80 degrees latitude
    
    #Determine number of ion and neutral species
    num_ionic_species = len(ions)
    num_neutral_species = len(neutrals)
    
    # plot input values - can take up to 7-neutral 7-ion species for diff colours
    if plots ==1:
    #        pl.figure(1)
        pwpl.input_plot(ions,electrons,neutrals,z,z_ext,A,radius)  
        pl.savefig(folder+'inputs_plot_%s_%s.png' %(planet_name,run_name))
        
    #calculation for centrifugal and gravitational accelleration - Dave's function
    phi,ag,ac = dipolefield.dipolefield_carley(radius,inner,z[-1],mass_planet,b0,lshell,rot_period,dipole_offset,numpoints)   
       
    
    #optional centrifugal force    
    if CF_flag != 0:
        ac=np.zeros(np.size(ag)) #if testing without centrifugal acceleration
        print('Centrifugal force removed')
    else:
        print('Centrifugal force included')
    
    
    #preallocate arrays that are filled but then updated at next step 
    dMdt = np.empty([len(z_ext),num_ionic_species])
    dMdt_tmp = np.empty([len(z_ext),num_neutral_species])
    dEdt = np.empty([len(z_ext),num_ionic_species])
    dEdt_tmp = np.empty([len(z_ext),num_neutral_species])
    dArhou = np.empty([len(z_ext),num_ionic_species])
    dArhou2= np.empty([len(z_ext),num_ionic_species])
    dPdr= np.empty([len(z_ext),num_ionic_species])
    dTdr= np.empty([len(z_ext),num_ionic_species])
    dEngdr= np.empty([len(z_ext),num_ionic_species])
    dkdr= np.empty([len(z_ext),num_ionic_species])
    ion_flux = np.empty([len(z_ext),its,num_ionic_species,])
    # ----------------------------------------------------------------------------   
    for i in range(1,its):
    #        print(i)
        for n in range(1,num_ionic_species+1):
            for p in range(1,num_neutral_species+1):
                # calculate momentum exchange rate for each neutral species
                dMdt_tmp[:,p-1] = pw.momentum_rate_1species(neutrals[p]["rho"],ions[n]["mass"],neutrals[p]["mass"],neutrals[p]["lambda"], e_charge,ions[n]["rho"][:,i-1],ions[n]["u"][:,i-1])
                dEdt_tmp[:,p-1] = pw.energy_rate_1species(ions[n]["rho"][:,i-1],neutrals[p]["rho"],ions[n]["mass"],neutrals[p]["mass"], neutrals[p]["lambda"], e_charge, neutrals[p]["T"],ions[n]["T"][:,i-1],ions[n]["u"][:,i-1],k_b)
            # sum each neutral species for each ionic species to get momentum exchange rate for each ionic species
            dMdt[:,n-1] = - np.sum(dMdt_tmp, axis=1)
            dMdt_tmp = np.empty([len(z_ext),num_neutral_species])    
            dEdt[:,n-1] = np.sum(dEdt_tmp, axis=1)
            dEdt_tmp = np.empty([len(z_ext),num_neutral_species]) 
                
    
    #numerically calculate differentials- central difference using roll function
    # TODO: Better numerical algorithm than the central difference!
    # ions
    
        for m in range(1,num_ionic_species+1):
            dArhou[:,m-1] = (np.roll(A * ions[m]["u"][:,i-1] *ions[m]["rho"][:,i-1],-1) - np.roll(A * ions[m]["u"][:,i-1] *ions[m]["rho"][:,i-1],1))/(2*dz)
            dArhou[0,m-1] = np.nan
            dArhou[-1,m-1] = np.nan
            
            dArhou2[:,m-1]= (np.roll(A * ions[m]["u"][:,i-1]**2 * ions[m]["rho"][:,i-1],-1)-np.roll(A * ions[m]["u"][:,i-1]**2 * ions[m]["rho"][:,i-1],1))/(2*dz)
            dArhou2[0,m-1] = np.nan
            dArhou2[-1,m-1] = np.nan
         
            dPdr[:,m-1] = (np.roll(ions[m]["P"][:,i-1],-1)-np.roll(ions[m]["P"][:,i-1],1))/(2*dz)
            dPdr[0,m-1] = np.nan
            dPdr[-1,m-1] = np.nan 
             
            dTdr[:,m-1] = (np.roll(ions[m]["T"][:,i-1],-1)-np.roll(ions[m]["T"][:,i-1],1))/(2*dz)  
            dTdr[0,m-1] = np.nan
            dTdr[-1,m-1] = np.nan 
            
            dEngdr[:,m-1]= (np.roll((0.5 * A * ions[m]["u"][:,i-1]**3 * ions[m]["rho"][:,i-1] + gamma/(gamma-1) * A * ions[m]["u"][:,i-1] * ions[m]["P"][:,i-1]),-1)-np.roll((0.5 * A * ions[m]["u"][:,i-1]**3 * ions[m]["rho"][:,i-1] + gamma/(gamma-1) * A * ions[m]["u"][:,i-1] * ions[m]["P"][:,i-1]),1))/(2*dz)
            dEngdr[0,m-1] = np.nan
            dEngdr[-1,m-1] = np.nan 
                    
            dkdr[:,m-1] = (np.roll(ions[m]["kappa"][:,i-1],-1)-np.roll(ions[m]["kappa"][:,i-1],1))/(2*dz)
            dkdr[0,m-1] = np.nan
            dkdr[-1,m-1] = np.nan
        
        #uses both ions for one term - can take as many ion species as it wants
        dEdr = (np.roll(pw.E_second_term(electrons,ions,dMdt,num_ionic_species,i),-1)-np.roll(pw.E_second_term(electrons,ions,dMdt,num_ionic_species,i),1))/(2*dz)
        dEdr[0] = np.nan
        dEdr[-1] = np.nan
        
        #electrons
        dAudr = (np.roll(A*electrons["u"][:,i-1],-1)-np.roll(A*electrons["u"][:,i-1],1))/(2*dz)#np.roll(a, shift, axis=None)
        dAudr[0] = np.nan
        dAudr[-1] = np.nan
    
        dTedr = (np.roll(electrons["T"][:,i-1],-1)-np.roll(electrons["T"][:,i-1],1))/(2*dz)  
        dTedr[0] = np.nan
        dTedr[-1] = np.nan    
        
        dPrhou2 = (np.roll(electrons["P"][:,i-1]+electrons["rho"][:,i-1]*electrons["u"][:,i-1]**2,-1)-np.roll(electrons["P"][:,i-1]+electrons["rho"][:,i-1]*electrons["u"][:,i-1]**2,1))/(2*dz)
        dPrhou2[0] = np.nan
        dPrhou2[-1] = np.nan
        
        dkdr_e = (np.roll(electrons["kappa"][:,i-1],-1)-np.roll(electrons["kappa"][:,i-1],1))/(2*dz) 
        dkdr_e[0] = np.nan
        dkdr_e[-1] = np.nan 
        
        #others    
        dAdr = (np.roll(A,-1)-np.roll(A,1))/(2*dz)
        dAdr[0] = np.nan
        dAdr[-1] = np.nan
      
        dakTdr =np.zeros(np.size(z_ext))
        dakTedr = np.zeros(np.size(z_ext))
        
        #parallel electric field - initial (short assumption)
        E[2:-2,i] = pw.E_parallel_short(e_charge, electrons["n"][2:-2,i-1].T, dPrhou2[2:-2], A[2:-2].T, dAdr[2:-2], electrons["rho"][2:-2,i-1].T, electrons["u"][2:-2,i-1].T) + 1/(e_charge*electrons["n"][2:-2,i-1]) * dEdr[2:-2].T
        E[0:2,i]=pw.extrap_start(E[2:-2,i])
        E[-2:,i]=pw.extrap_end(E[2:-2,i])
    
        # updated values for next step
        #ions
        for l in range(1,num_ionic_species+1):
            #Mass conservation equation
            ions[l]["rho"][2:-2,i] = pw.density_dt_ion(dt,A[2:-2],ions[l]["S"][2:-2],ions[l]["rho"][2:-2,i-1],dArhou[2:-2,l-1].T)
            ions[l]["rho"][0:2,i]= pw.extrap_start(ions[l]["rho"][2:-2,i])
            ions[l]["rho"][-2:,i]= pw.extrap_end(ions[l]["rho"][2:-2,i]) 
            ions[l]["n"][2:-2,i] = ions[l]["rho"][2:-2,i] / ions[l]["mass"]
            ions[l]["n"][0:2,i]=pw.extrap_start(ions[l]["n"][2:-2,i])
            ions[l]["n"][-2:,i]=pw.extrap_end(ions[l]["n"][2:-2,i]) 
            #Momentum conservation equation
            ions[l]["u"][2:-2,i] = pw.velocity_dt_ion(dt,A[2:-2],ions[l]["rho"][2:-2,i-1],ions[l]["rho"][2:-2,i],dArhou2[2:-2,l-1].T,dPdr[2:-2,l-1].T,ions[l]["mass"],E[2:-2,i],e_charge,-ag,dMdt[2:-2,l-1],ions[l]["u"][2:-2,i-1],ions[l]["S"][2:-2],ac)
            ions[l]["u"][0:2,i]= pw.extrap_start(ions[l]["u"][2:-2,i])
            ions[l]["u"][-2:,i]= pw.extrap_end(ions[l]["u"][2:-2,i])  
            #Energy conservation equation
            ions[l]["P"][2:-2,i] = pw.pressure_dt_ion(dt,A[2:-2],gamma,ions[l]["rho"][2:-2,i-1],ions[l]["rho"][2:-2,i],ions[l]["u"][2:-2,i-1],ions[l]["u"][2:-2,i],ions[l]["mass"], e_charge,E[2:-2,i],-ag,dEdt[2:-2,l-1],dMdt[2:-2,l-1],ions[l]["P"][2:-2,i-1],dEngdr[2:-2,l-1].T, dakTdr[2:-2].T,ions[l]["S"][2:-2],ac) 
            ions[l]["P"][0:2,i]= pw.extrap_start(ions[l]["P"][2:-2,i])
            ions[l]["P"][-2:,i]= pw.extrap_end(ions[l]["P"][2:-2,i]) 
            ions[l]["T"][2:-2,i] = pw.plasma_temperature(ions[l]["n"][2:-2,i],k_b,ions[l]["P"][2:-2,i])
            ions[l]["T"][0,i]=400
            ions[l]["T"][1,i]=400
            ions[l]["T"][-2:,i]= pw.plasma_temperature(ions[l]["n"][-2:,i],k_b,ions[l]["P"][-2:,i])
            #Heat conductivities
            ions[l]["kappa"][2:-2,i] = pw.heat_conductivity(ions[l]["T"][2:-2,i],e_charge,ions[l]["mass"],m_p)
            ions[l]["kappa"][0:2,i]= pw.extrap_start(ions[l]["kappa"][2:-2,i])
            ions[l]["kappa"][-2:,i]= pw.extrap_end(ions[l]["kappa"][2:-2,i]) 
         
        #Electrons
        #mass conservation
        electrons["rho"][2:-2,i] = pw.density_dt_electron(electrons,ions,num_ionic_species,i)
        electrons["rho"][0:2,i]=pw.extrap_start(electrons["rho"][2:-2,i])
        electrons["rho"][-2:,i]=pw.extrap_end(electrons["rho"][2:-2,i])      
        electrons["n"][2:-2,i] = electrons["rho"][2:-2,i] / electrons["mass"]
        electrons["n"][0:2,i]= pw.extrap_start(electrons["n"][2:-2,i])
        electrons["n"][-2:,i]= pw.extrap_end(electrons["n"][2:-2,i])   
        #momentum conservation  
        electrons["u"][2:-2,i] = pw.velocity_dt_electron(electrons,ions,num_ionic_species,i,FAC,e_charge)
        electrons["u"][0:2,i]= pw.extrap_start(electrons["u"][2:-2,i])
        electrons["u"][-2:,i]= pw.extrap_end(electrons["u"][2:-2,i]) 
        #energy conservation
        electrons["T"][2:-2,i] =  pw.temperature_dt_electron(dt,gamma,electrons["mass"],k_b,A[2:-2],electrons["rho"][2:-2,i-1],electrons["u"][2:-2,i-1],electrons["T"][2:-2,i-1],electrons["S"][2:-2],dTedr[2:-2].T,dAudr[2:-2].T,dakTedr[2:-2].T)
        electrons["T"][0:2,i]=pw.extrap_start(electrons["T"][2:-2,i])
        electrons["T"][-2:,i]=pw.extrap_end(electrons["T"][2:-2,i])  
        electrons["P"][2:-2,i] = pw.plasma_pressure(electrons["rho"][2:-2,i]/electrons["mass"], k_b,electrons["T"][2:-2,i])
        electrons["P"][0:2,i]=pw.extrap_start(electrons["P"][2:-2,i])
        electrons["P"][-2:,i]=pw.extrap_end(electrons["P"][2:-2,i])      
        #Heat conductivities
        electrons["kappa"][2:-2,i] = pw.heat_conductivity_electrons((electrons["T"][2:-2,i]),e_charge,gamma)
        electrons["kappa"][0:2,i]= pw.extrap_start(electrons["kappa"][2:-2,i])
        electrons["kappa"][-2:,i]= pw.extrap_end(electrons["kappa"][2:-2,i]) 
        
        #Electron and ion Flux    
        for w in range(1,num_ionic_species+1):
            ion_flux[2:-2,i,w-1] = ions[w]["n"][2:-2,i]*ions[w]["u"][2:-2,i] * A[2:-2] 
        e_flux[2:-2,i] = pw.electron_flux_e(electrons,i)*A[2:-2]
        
    
    
    #plotting 
    if plots ==1:
        print('Plots on screen')
    #        pl.figure(2)    
        pwpl.results_plot(z,num_ionic_species,e_charge,E[2:-2,-1],ions,electrons,ac,ag,e_flux,ion_flux) # z_ext,radius removed from second and third place
        pl.savefig(folder+'overview_results_plot_%s_%s.png' %(planet_name,run_name))
    #        pl.figure(3)
        pwpl.species_plot(z,z_ext,electrons,radius)
        pl.savefig(folder+'species_plot_%s_electrons_%s.png' %(planet_name,run_name))
        for q in range(1,num_ionic_species+1):
    #            pl.figure(q+3)
            pwpl.species_plot(z,z_ext,ions[q],radius)
            pl.savefig(folder+'species_plot_%s_%s_%s.png' %(planet_name,ions[q]['name'],run_name))
        print('Plots saved to: %s' %folder)
              
    
    #Calculating total particle source
    ind = np.max(np.where(z<10000000)) +1 #2# at altitude of 10000km 
    print(ind)
    elef = e_flux[ind,-1]/A[ind] # electron flux at altitude of 10000km
    ionf=np.empty([num_ionic_species])
    for v in range(1,num_ionic_species+1):
        ionf[v-1] = ion_flux[ind,-1,v-1] / A[ind]# ion flux at altitude of 10000km
    
    arc = 10/360 * 2*np.pi*(radius+10000000) #auroral arc of 2 deg width #2
    circ = 2*np.pi*(radius+10000000)*np.sin(np.deg2rad(14)) #auroral arc centred on 14deg #2
    
    elecTPS = elef*arc*circ*2
    ion1TPS = ionf[0]*arc*circ*2
    ion2TPS = ionf[1]*arc*circ*2 
    print('Total particle source [/s]')
    print(2*elecTPS)
    
    elecTMS = elecTPS * electrons['mass']
    ion1TMS = ion1TPS * ions[1]["mass"]
    ion2TMS = ion2TPS * ions[2]["mass"]
    
    print('Total Mass Source [kg/s]')
    print(elecTMS+ion1TMS+ion2TMS)
    # breakpoint()
    
    pl.figure(6)
    # for k in range(1,num_ionic_species+1):
    #     pl.plot(z/1000,ions[k]["S"][2:-2])
    pl.plot(z/1000,electrons["S"][2:-2])
    pl.xlabel('Distance (km)')
    pl.ylabel('Mass Production Rate')
    pl.yscale('log')
    pl.title('Test Plot for Electron Mass Production Rate')
    pl.savefig(folder+'test_plot_%s_%s_%s.png' %(planet_name,ions[q]['name'],run_name))
    pl.show()
    
    
    #pwpl.plot_me_quick(z/1000,E[2:-2],'Distance Along Field Line (km)','E Parallel (V/m)','on')
    #pl.savefig(folder+'test_plot_%s_%s.png' %(planet_name,run_name))
    #for some reason this is just making another H3+ graph??
    # try manually writing out own plotting function
    
    if saves ==1:
        
        #create and open file
        fid = folder+"polar_wind_output_%s_%s.txt" %(planet_name,run_name)
        f= open(fid,"w+")
        print('Data saved to file: %s ' %fid)
        f.write("-------SETUP------- \n")
        f.write("Planet=%s\n" %(planet_name))
        f.write("Time Step=%s\n" %(dt))
        f.write("Spatial Step=%s\n" %(dz))
        f.write("Iterations=%s\n" %(its))
        f.write("Outer limit (radii)=%s\n" %(rs))
        f.write("L-Shell=%s\n" %(lshell))
        f.write("Field Aligned Currents removed:1,included:0 =%s\n" %(FAC_flag))
        f.write("Centrifugal Stress  removed:1,included:0 =%s\n" %(CF_flag))
        f.write("Number of Ionic species=%s\n" %(num_ionic_species))
        f.write("Number of Neutral species=%s\n" %(num_neutral_species))
        for s in range(1,num_ionic_species+1):
            f.write("Ion Species %d: %s\n" %(s,ions[s]["name"]))
        for w in range(1,num_neutral_species+1):
            f.write("Neutral Species %d: %s\n" %(w,neutrals[w]["name"]))     
        f.write("-------VECTORS------- \n")    
        
        f.write("Spatial Grid: \n")
        for cc in range(0,len(z)):
            f.write('%s,' %(z[cc]))
        f.write("\n")
        
        f.write("Gravitational Acceleration: \n")
        for aa in range(0,len(ag)):
            f.write('%s,' %(ag[aa]))
        f.write("\n")  
        
        f.write("Centrifugal Acceleration: \n")
        for bb in range(0,len(ac)):
            f.write('%s,' %(ac[bb]))
        f.write("\n")
        
        f.write("Electric Field: \n")
        for dd in range(0,len(E[2:-2,-1])):
            f.write('%s,' %(E[2:-2,-1][dd]))
        f.write("\n")
        
        f.write("-------ELECTRONS------- \n") 
        
        f.write("rho (kg/m^3): \n")
        for ee in range(0,len(electrons["rho"][2:-2,-1])):
            f.write('%s,' %(electrons["rho"][2:-2,-1][ee]))
        f.write("\n")
        
        f.write("n (/m^3): \n")
        for ff in range(0,len(electrons["n"][2:-2,-1])):
            f.write('%s,' %(electrons["n"][2:-2,-1][ff]))
        f.write("\n")
        
        f.write("u (m/s): \n")
        for gg in range(0,len(electrons["u"][2:-2,-1])):
            f.write('%s,' %(electrons["u"][2:-2,-1][gg]))
        f.write("\n")
        
        f.write("P (Pa): \n")
        for hh in range(0,len(electrons["P"][2:-2,-1])):
            f.write('%s,' %(electrons["P"][2:-2,-1][hh]))
        f.write("\n")  
        
        f.write("T (K): \n")
        for ii in range(0,len(electrons["T"][2:-2,-1])):
            f.write('%s,' %(electrons["T"][2:-2,-1][ii]))
        f.write("\n")
        
        f.write("kappa: \n")
        for jj in range(0,len(electrons["kappa"][2:-2,-1])):
            f.write('%s,' %(electrons["kappa"][2:-2,-1][jj]))
        f.write("\n")
        
        
        for kk in range(1,num_ionic_species+1):
            f.write("-------%s------- \n" %ions[kk]["name"]) 
        
            f.write("rho (kg/m^3): \n")
            for ee in range(0,len(ions[kk]["rho"][2:-2,-1])):
                f.write('%s,' %(ions[kk]["rho"][2:-2,-1][ee]))
            f.write("\n")
            
            f.write("n (/m^3): \n")
            for ff in range(0,len(ions[kk]["n"][2:-2,-1])):
                f.write('%s,' %(ions[kk]["n"][2:-2,-1][ff]))
            f.write("\n")
            
            f.write("u (m/s): \n")
            for gg in range(0,len(ions[kk]["u"][2:-2,-1])):
                f.write('%s,' %(ions[kk]["u"][2:-2,-1][gg]))
            f.write("\n")
            
            f.write("P (Pa): \n")
            for hh in range(0,len(ions[kk]["P"][2:-2,-1])):
                f.write('%s,' %(ions[kk]["P"][2:-2,-1][hh]))
            f.write("\n")  
            
            f.write("T (K): \n")
            for ii in range(0,len(ions[kk]["T"][2:-2,-1])):
                f.write('%s,' %(ions[kk]["T"][2:-2,-1][ii]))
            f.write("\n")
            
            f.write("kappa: \n")
            for jj in range(0,len(ions[kk]["kappa"][2:-2,-1])):
                f.write('%s,' %(ions[kk]["kappa"][2:-2,-1][jj]))
            f.write("\n")
        
        f.close()

# if __name__ == "__main__":
# 	import sys
# 	polar_wind_outflow(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],
# 					   sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10])