# DISCRIPTION:Analytical model of polar wind at Saturn from equations in Glocer+ 2007
#			 now fully generalised for any number of ionic or neutral species
# AUTHOR: C.J.Martin
# INSTITUTION: Lancaster University

# INPUTS: dipolefield.py, pw.py, pw_plotting_tools.py and planet.py modules required
#		 planet_name = string with planet name (saturn or jupiter)
#		 dt = time step (0.01)
#		 its = iterations (100)
#		 rs = rough estimate of how far in planetary radii you want to go (add a few more fo jupiter)
#		 lshell = lshell of field line (33 is around 80 deg lat for both)
#		 FAC_flag = flag for field aligned currents (0=no change, 1=removed)
#		 CF_flag = flag for centrifugal force (0=no change, 1=removed)
#		 plots = do you want plotting (1=yes, 0=no)
#		 saves = do you want to save the output to file (1=yes, 0=no)

# OUTPUT: Plots for now

# NOTES: 
# All in SI units unless specified
# Part of the lancaster polar wind model
# 2019

#polar_wind_outflow('jupiter',0.01,100,3,33,0,0,1,0,'anything')
def polar_wind_outflow(planet_name,dt,its,rs,lshell,FAC_flag,CF_flag,plots,saves,run_name):
#def polar_wind_outflow(planet_name,**kwargs):	
	# ---------------------------------Imports-------------------------------------
	import numpy as np
	import matplotlib.pyplot as pl
	import dipolefield
	import pw
	import planet
	import pw_plotting_tools as pwpl
	# ---------------------------------Start Main-------------------------------------
	# SECTION FOR CALCULATING RADIAL DIFFERENCES- JUPITER  
	x = [75,75.5,76,76.5,77] #- Auroral oval 
	y = [19.773,27.073,40.863,63.170,89.031]#Vogt mapping to 
	p = np.polyfit(x,y,2) # fits a polynomial to above x and y
	x2 = np.linspace(75,77,101) # for from 75 degrees to 77 degrees in 100 field lines
	y_eq = np.polyval(p,x2) # 100 unequally spaced lengths at equator - find where the x2 lines hit the equator
	x_is = np.linspace(73.01,76.99,100) # equally spaced in ionosphere
	area_eq = np.empty([100]) # initialise the equatorial area array
	# Calculate the area of the equator that the area in the ionsphere maps to:
	for m in range(1,101):
		area_eq[m-1] = np.pi*(y_eq[m]*7.1492e7)**2 *((1)/(2* np.pi*(7.1492e7+20000000)))- np.pi*(y_eq[m-1]*7.1492e7)**2 *((1)/(2* np.pi*(7.1492e7+20000000)))
	area_is = 25000 
	# If we include the current density, get the values from pw
	j = pw.current_density_ray2015()*1e-6
	print(area_eq)
	# TODO: make a if planets = saturn do below/ if planet = Jupiter do above section
	# SECTION FOR CALCULATING RADIAL DIFFERENCES- Saturn
	#y=[10.891181988742963,11.172607879924952,11.454033771106939,11.791744840525325,12.101313320825513,12.326454033771105,12.607879924953094,12.833020637898684,13.114446529080674,13.395872420262663,13.649155722326451,13.90243902439024,14.127579737335832,14.380863039399625,14.577861163227016,14.831144465290805,15.112570356472794,15.309568480300186,15.647279549718572,15.872420262664164, 16.097560975609753,16.37898686679174,16.575984990619133,16.74484052532833,17.082551594746715,17.335834896810507,17.504690431519695,17.701688555347086,17.926829268292682,18.15196998123827,18.348968105065662]
	#x=[72.95902883156297, 73.0804248861912, 73.20182094081942, 73.32321699544764, 73.44461305007587, 73.53566009104703, 73.59635811836115, 73.65705614567526, 73.77845220030349, 73.86949924127465, 73.96054628224582, 74.05159332321699, 74.14264036418817, 74.20333839150227, 74.32473444613049, 74.38543247344461, 74.47647951441579, 74.597875569044, 74.65857359635811, 74.77996965098635, 74.87101669195751, 74.96206373292867, 75.05311077389985, 75.11380880121396, 75.23520485584218, 75.3566009104704, 75.47799696509864, 75.5690440060698, 75.66009104704096, 75.7814871016692, 75.90288315629742]
	#p = np.polyfit(x,y,2)
	#x2 = np.linspace(75,77,101)
	#y_eq = np.polyval(p,x2) # 100 unequally spaced lengths at equator
	#x_is = np.linspace(75.01,76.99,100) # equally spaced in ionosphere
	#area_is = 23000 
	#area_eq = np.empty([100])
	#for m in range(1,100):
	#	area_eq[m-1] = np.pi*(y_eq[m]*5.8232e7)**2 *((1)/(2* np.pi*(5.8232e7+20000000)))- np.pi*(y_eq[m-1]*5.8232e7)**2 *((1)/(2* np.pi*(5.8232e7+20000000)))
	# 

	# Initialize the arrays for equatorial fluxes
	eflux_eq = np.empty([len(x_is)])
	ionfH_eq = np.empty([len(x_is)])
	ionfH3_eq = np.empty([len(x_is)])
	elecEM = np.empty([len(x_is)])
	ion1EM = np.empty([len(x_is)])
	ion2EM = np.empty([len(x_is)])
	tot_mass_eq = np.empty([len(x_is)])
	
	# For each of the ionospheric field lines
	for h in range(1,len(x_is)+1):
		# Make sure the read in values are actually in the correct format
		dt = float(dt)
		its = int(its)
		rs = int(rs)
		lshell = int(lshell) #y_eq[h-1]

		# TODO: Make folder as a read in or make automatically in same repo?
		#	Main folder to save plots and data to
		folder = '/Users/carleymartin/Desktop/polar_wind_output/'
		# ---------------------------------Grid set up-------------------------------------	
		#800 spatial grids per Rs roughly - a bit more to round the numbers
		#	rs = 3 #how far out in rs do you want to go - roughly?
		inner = 1400000  #lower boundary
		numpoints = rs * 800 
		dz = 75000.0 # grid spacing
	
		# TODO: Explain ghost points, I have a ghost of a memory of why I used them...
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
	
		#ghost points only
		gb_z = np.linspace(z_ext[0],z_ext[1],2)
		ge_z = np.linspace(z_ext[-2],z_ext[-1],2)
		gb_x = np.linspace(-1,0,2)
		ge_x = np.linspace(numpoints,numpoints+1,2)
	
		ghosts = [gb_z,ge_z,gb_x,ge_x] #combine to read into input module

		# -------------------------Physical Constants----------------------------------
		# Constants set-up: don't change so just hard coded in
		e_charge = 1.60217662*10**-19  #C - electron charge
		m_p = 1.6726*10**-27 # kg - mass of a proton
		m_e = 9.10938356 * 10**-31 # kg - mass of an electron
		k_b = 1.38064852*10**-23 #m2 kg s-2 K-1  - boltzmann constant
		gamma = 5/3 #specific heat ratio
	
		phys_consts = [m_e,m_p,k_b,e_charge,gamma] #combine to read into input module

		# -----------------------READ IN INITIAL CONDITIONS----------------------------
		# choose planet for input module - name of planet no capitals
		# initial conditions are dictionaries for each species
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
			FAC = np.zeros(np.size(z_ext)) #if testing with/without field aligned currents
			print('Field aligned currents removed')
		else:
			print('Field aligned currents included')
	  
		#option to include rudimental upward and downward current (very basic)
		#if x_is[h-1] >= 76.0:	
		#	FAC = 1e-11 * (A[2]/A) # 1-7 microAm-1 from Ray+ 2009 1e-11 5e-9
		#else:
		#	FAC = - 0.3 * 1e-11  * (A[2]/A)
		#	
	
		# TODO: I'm not sure if this is overwriting the above if removed is chosen
		# option for RAY2015 current density profile
		# FAC = j[h] * (A[2]/A)
		
		# preallocate arrays not used in initial conditions - electric field and flux
		E = np.empty([len(z)+4,its])
		E[:,0] = np.nan # no initial values for electric field
		e_flux = np.empty([len(z)+4,its])
		e_flux[:,0] = np.nan # no initial values for electron flux
		ion1_flux = np.empty([len(z)+4,its])
		ion1_flux[:,0] = np.nan # no initial values for electron flux
		ion2_flux = np.empty([len(z)+4,its])
		ion2_flux[:,0] = np.nan # no initial values for electron flux
	
		# unpack constants from planet
		radius = consts[0]
		mass_planet= consts[1]
		b0= consts[2]
		rot_period= consts[3]
		dipole_offset= consts[4]
		#	g= consts[5]
		#	lshell=33 #equivalent to 80 degrees latitude
	
		#Determine number of ion and neutral species
		num_ionic_species = len(ions)
		num_neutral_species = len(neutrals)
	
		# plot input values - can take up to 7-neutral 7-ion species for diff colours
		if plots ==1:
		#		pl.figure(1)
			pwpl.input_plot(ions,electrons,neutrals,z,z_ext,A,radius)  
			pl.savefig(folder+'inputs_plot_%s_%s.png' %(planet_name,run_name))
		
		#calculation for centrifugal and gravitational accelleration - Dave's function
		phi,ag,ac = dipolefield.dipolefield_carley(radius,inner,z[-1],mass_planet,b0,lshell,rot_period,dipole_offset,numpoints)   
	   
		#optional centrifugal force	
		if CF_flag != 0:
			ac=np.zeros(np.size(ag)) #if testing with/without centrifugal acceleration
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
		
		# This is where the heavy lifting happens
		# ----------------------------------------------------------------------------   
		for i in range(1,its):
		#		print(i)
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
			dAudr = (np.roll(A*electrons["u"][:,i-1],-1)-np.roll(A*electrons["u"][:,i-1],1))/(2*dz)
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
		#		pl.figure(2)	
			pwpl.results_plot(z,z_ext,radius,num_ionic_species,e_charge,E[2:-2,-1],ions,electrons,ac,ag,e_flux,ion_flux)
			pl.savefig(folder+'overview_results_plot_%s_%s.png' %(planet_name,run_name))
		#		pl.figure(3)
			pwpl.species_plot(z,z_ext,electrons,radius)
			pl.savefig(folder+'species_plot_%s_electrons_%s.png' %(planet_name,run_name))
			for q in range(1,num_ionic_species+1):
		#			pl.figure(q+3)
				pwpl.species_plot(z,z_ext,ions[q],radius)
				pl.savefig(folder+'species_plot_%s_%s_%s.png' %(planet_name,ions[q]['name'],run_name))
			print('Plots saved to: %s' %folder)
		# plot single value with r if needed
		#pwpl.plot_single(z,z_ext,item,radius)
		#return E,z,z_ext,radius
	
		#Calculating total particle source
		ind = np.max(np.where(z<20000000)) +1 # at altitude of 10000km
		elef = e_flux[ind,-1]/A[ind] # electron flux at altitude of 10000km
		ionf=np.empty([num_ionic_species])
		for v in range(1,num_ionic_species+1):
			ionf[v-1] = ion_flux[ind,-1,v-1] / A[ind]# ion flux at altitude of 10000km
	
		arc = 2/360 * 2*np.pi*(radius+20000000) #auroral arc of 2 deg width
		circ = 2*np.pi*(radius+20000000)*np.sin(np.deg2rad(14)) #auroral arc centred on 14deg
	
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
			
		
		# SECTION FOR CALCULATING RADIAL DIFFERENCES  
		eflux_eq[h-1] = (elef/area_eq[h-1]) *area_is
		ionfH_eq[h-1] = (ionf[0]/area_eq[h-1]) *area_is
		ionfH3_eq[h-1] = (ionf[1]/area_eq[h-1]) *area_is
		
	
		elecEM[h-1] = eflux_eq[h-1] * electrons['mass']
		ion1EM[h-1] = ionfH_eq[h-1] * ions[1]["mass"]
		ion2EM[h-1] = ionfH3_eq[h-1] * ions[2]["mass"]
		tot_mass_eq[h-1] = elecEM[h-1]+ion1EM[h-1]+ion2EM[h-1]
	
	pl.figure(figsize=(8,10))

	ax1 = pl.subplot(3,1,1)
	pl.yscale('log')
	pl.plot(y_eq[0:100],eflux_eq,linestyle='-',color=[0.2,1,0.2])
	pl.ylabel('a) Electron Flux (equator)\n $(m^{-2}$ $s^{-1})$')
	pl.xlim([19,90])
	ax1.set_xlabel('')
	ax1.xaxis.set_ticklabels([])
	pl.grid('on') 
	ax2 = ax1.twiny()
	ax2.plot(y_eq[0:100],eflux_eq,linestyle= '-',color=[0.2,1,0.2])
	ax2.set_xlabel('Ionospheric Latitude (degrees)')
	ax2.set_xlim([19,90])
	pl.xticks([19.77,27.07,40.86,63.17,89.03], ['75.0','75.5','76.0','76.5','77.0'])
	ax1.set_xlabel('')
	ax1.legend(['$e^-$'])

	ax7 = pl.subplot(3,1,2)
	pl.yscale('log')
	pl.plot(y_eq[0:100],ionfH_eq,color=[0,0.2,0.5])
	pl.plot(y_eq[0:100],ionfH3_eq,color=[0.2,0.6,1])	
	pl.ylabel('b) Ion Fluxes (equator) \n $(m^{-2}$ $s^{-1})$')
	pl.xlabel('')
	pl.xlim([19,90])
	pl.grid('on') 
	ax8 = ax7.twiny()
	ax8.plot(y_eq[0:100],ionfH_eq,color=[0,0.2,0.5])
	ax8.set_xlabel('')
	ax8.set_xlim([19,90])
	ax8.xaxis.set_ticklabels([])
	ax7.xaxis.set_ticklabels([])
	ax7.legend(['$H^+$','$H_3^+$'])
		
	
	ax5 = pl.subplot(3,1,3)
	pl.yscale('log')
	pl.plot(y_eq[0:100],tot_mass_eq,color=[0.0,0,0.0])
	pl.ylabel('c) Mass Flux \n $(kg$ $m^{-2}$ $s^{-1})$')
	pl.xlim([19,90])
	pl.grid('on') 
	ax6 = ax5.twiny()
	ax6.plot(y_eq[0:100],tot_mass_eq,color=[0.0,0,0.0])
	ax6.set_xlim([19,90])
	ax5.set_xlabel('Radial Distance at the Equator $(R_J)$')
	ax6.set_xlabel('')
	ax6.xaxis.set_ticklabels([])


	pl.subplots_adjust(hspace=0.0,wspace=0.5)

	# TODO: Automate from the planets value (if planet = saturn/ plot this)
	#SATURN PLOT
	#
	#pl.figure(figsize=(8,10))
	#
	#ax1 = pl.subplot(3,1,1)
	#pl.yscale('log')
	#pl.plot(y_eq[0:100],eflux_eq,linestyle='-',color=[0.2,1,0.2])
	#pl.ylabel('a) Electron Flux (equator)\n $(m^{-2}$ $s^{-1})$')
	#pl.xlim([16,21])
	#ax1.set_xlabel('')
	#ax1.xaxis.set_ticklabels([])
	#pl.grid('on') 
	#ax2 = ax1.twiny()
	#ax2.plot(y_eq[0:100],eflux_eq,linestyle= '-',color=[0.2,1,0.2])
	#ax2.set_xlabel('Ionospheric Latitude (degrees)')
	#ax2.set_xlim([16,21])
	#pl.xticks([16.4,17.6,18.7,19.7,20.6], ['75.0','75.5','76.0','76.5','77.0'])
	#ax1.set_xlabel('')
	#ax1.legend(['$e^-$'])
	#
	#ax7 = pl.subplot(3,1,2)
	#pl.yscale('log')
	#pl.plot(y_eq[0:100],ionfH_eq,color=[0,0.2,0.5])
	#pl.plot(y_eq[0:100],ionfH3_eq,color=[0.2,0.6,1])	
	#pl.ylabel('b) Ion Fluxes (equator) \n $(m^{-2}$ $s^{-1})$')
	#pl.xlabel('')
	#pl.xlim([16,21])
	#pl.grid('on') 
	#ax8 = ax7.twiny()
	#ax8.plot(y_eq[0:100],ionfH_eq,color=[0,0.2,0.5])
	#ax8.set_xlabel('')
	#ax8.set_xlim([16,21])
	#ax8.xaxis.set_ticklabels([])
	#ax7.xaxis.set_ticklabels([])
	#ax7.legend(['$H^+$','$H_3^+$'])
	#		
	#	
	#ax5 = pl.subplot(3,1,3)
	#pl.yscale('log')
	#pl.plot(y_eq[0:100],tot_mass_eq,color=[0.0,0,0.0])
	#pl.ylabel('c) Mass Flux \n $(kg$ $m^{-2}$ $s^{-1})$')
	#pl.xlim([16,21])
	#pl.grid('on') 
	#ax6 = ax5.twiny()
	#ax6.plot(y_eq[0:100],tot_mass_eq,color=[0.0,0,0.0])
	#ax6.set_xlim([16,21])
	#ax5.set_xlabel('Radial Distance at the Equator $(R_S)$')
	#ax6.set_xlabel('')
	#ax6.xaxis.set_ticklabels([])
	#
	#
	#pl.subplots_adjust(hspace=0.0,wspace=0.5)


if __name__ == "__main__":
	import sys
	polar_wind_outflow(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],
					   sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10])
