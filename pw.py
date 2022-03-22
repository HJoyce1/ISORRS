# Functions for polar wind model
# All in SI units unless specified
# Part of the lancaster polar wind model
# Author: C.J. Martin @ Lancaster University
# 2019

#EDITED BY H.S. Joyce, Lancaster University
#==============================================================================
def red_mass(m_i,m_j): 
#Returns: reduced mass
#Requires: mass of two species - ionic-i and neutral-j
#Ref: SCHUNK + NAGY 2000
    return (m_i*m_j)/(m_i+m_j) #ion mass x neutral mass / ion mass + neutral mass
#==============================================================================    
def collision_freq(rho_n,m_i,m_j,lambda_n, e_charge):
#Returns: collision frequency between two species - ionic-i and neutral-j
#Requires: neutral mass density, ionic mass, neutral mass, neutral polarisability
#            and electron charge
#Ref: SCHUNK + NAGY 2000
    import numpy as np
    return 2.21*np.pi*(rho_n/(m_i+m_j))*((lambda_n*e_charge**2)/\
                       red_mass(m_i,m_j))**0.5 
        # 2.21 x pi x (mass density / ion mass + neutral mass) x ((neutral polarisability x charge of electron^2) / 
    #     reduced mass)^0.5
#==============================================================================
def momentum_rate_1species(rho_n,m_i,m_j,lambda_n, e_charge,rho_i,u_i):
#Returns:momentum exhange rate for one ionic species and one neutral species
#Requires:neutral mass density, ionic mass, neutral mass, neutral polarisability 
#         electron charge, ionic mass density and ion speed (assuming neutral 
#         velocity is 0)
#Ref: SCHUNK + NAGY 2000
    return rho_i*collision_freq(rho_n,m_i,m_j,lambda_n, e_charge)*u_i
# ion mass density x collision frequency x speed of ions 
#==============================================================================
def energy_rate_1species(rho_i,rho_n,m_i,m_j,lambda_n, e_charge,T_j,T_i,u_i,k_b):
#Returns:energy exhange rate for one ionic species and one neutral species
#Requires:ionic mass density, neutral mass density, ionic mass, neutral mass,  
#         neutral polarisability, electron charge, neutral temp, ion temp,
#         ion speed and boltzmann constant (assuming neutral velocity is 0)
#Ref: SCHUNK + NAGY 2000
    return ((rho_i * collision_freq(rho_n,m_i,m_j,lambda_n, e_charge))/\
            (m_i + m_j)) * (3* k_b * (T_j-T_i) + m_j*u_i**2)
        # ion mass density x collision freq / (ion mass + neutral mass) x 3xboltzmann x (temp neutrals - temp ions) +
        # neutral mass x neutral velocity^2
#==============================================================================
def heat_conductivity(T,e_charge,m_i,m_p):
#Returns: ion heat conductivity
#Requires: electron temp. electron charge, mass ion, mass proton
#Ref: Banks&Kockarts1973, Moore+ 2008
    return (4.6*10**4 *(m_i/m_p)**-0.5 * T**2.5)*e_charge*100 #J m-1 s-1 K-1
#( 4.6x10^4 x (ion mass/ proton mass)^-0.5 x temp^2.5) x electron charge x100
#==============================================================================
def heat_conductivity_electrons(T,e_charge,gamma):
#Returns: electron heat conductivity
#Requires: electron temp. electron charge, specific heat
#Ref: Banks&Kockarts1973, Raitt+ 1975: 7.7 x 10^5
    return ((1.8*10**6)*(T**(5/2)))*e_charge*100 #J m-1 s-1 K-1
# ((1.8x10^6) x temp^5/2) x elecltron charge
#==============================================================================
#Function to calculate plasma pressure
def plasma_pressure(n,k_b,T):
#Returns: plasma pressure
#Requires: number density, boltzmann constant, temperature
#Ref: P = nkT  
    return (n * k_b * T)
# number density x boltzmanns x temp
#==============================================================================
def plasma_temperature(n,k_b,P):
#Returns: plasma temperature
#Requires: number density, boltzmann constant, pressure
#Ref: P = nkT  
    return (P/(n*k_b))
# pressure / number density x boltzmanns (pressure is from plasma pressure function)
#==============================================================================   
def E_second_term(electrons,ions,dMdt,num_ionic_species,j):
#Returns: extended electric field calculation 
#Requires: electron dictionary, ion dictionary, momentum exchange rate, number 
#       ionic species, iteration number
#Ref: Gombozi+1985, Glocer+ 2007 etc.
    import numpy as np
    E_tmp = np.empty([len(dMdt),num_ionic_species])
    for i in range(1,num_ionic_species+1):
        E_tmp[:,i-1] = (electrons["mass"]/ions[i]["mass"] * ((electrons["u"][:,j-1]-ions[i]["u"][:,j-1])*ions[i]["S"]-dMdt[:,i-1]))       
    return np.sum(E_tmp,axis=1)
# create empty array of length dMdt, number of ionic species
# fill as loop, electron mass / ion mass x( electron velocity - ion velocity) * ion mass production rate - dMdt
# sum it all up on -1 axis
# CHECK THIS ONE
#==============================================================================
def E_parallel_short(e_charge, n_e, dFdt,A, dAdr, rho_e,u_e):
#Returns: electric field, short term

#Requires: electron charge, electron number density, dFdt = d/dr(P_e - rho_e*u_e^2)
#          areaA, spatial area differential, electron mass density, electron velocity
#Ref: Gombozi+1985, Glocer+ 2007 etc.
    return (-1/(e_charge*n_e) * (dFdt + (dAdr/A)*rho_e*u_e**2))
# (-1 / electron charge x electron number density) x (dFdt + (dAdr/cross secion area)*mass density electronsxelectron speed^2
#CHECK THIS ONE
#==============================================================================
def plot_me_quick(x,y,xlabel,ylabel,grid):
#Returns: plot - also found in pw_plotting tools
#Requires: x,y, axis labels and grid 'on' or 'off'    
    import matplotlib.pyplot as pl
    pl.plot(x,y)
    pl.ylabel(ylabel)
    pl.xlabel(xlabel)
    pl.grid(grid) 
#==============================================================================
def density_dt_ion(dt,A,S,rho,dFdt):
#returns: ion density
#requires: time step, area vector, mass production, mass density
#           dFdt = d/dr(A rho u)
#ref: Gombozi+1985, Glocer+ 2007 etc.
    return (((dt/A)*(A*S-dFdt)) + rho)
# dt/A x (area of flux tube *mass production -dFdt) +mass density
#==============================================================================    
def density_dt_electron(electrons,ions,num_ionic_species,j):
#returns: electron density
#requires: electron dictionary, ion dictionary, number of ionic species, iteration
#ref: Gombozi+1985, Glocer+ 2007 etc.    
    import numpy as np
    term1 = np.empty([len(ions[1]["n"])-4,num_ionic_species])
    for i in range(1,num_ionic_species+1):
        term1[:,i-1] = ions[i]["rho"][2:-2,j]/ions[i]["mass"]
    return np.sum(term1,axis=1)*electrons["mass"]
# empty array of length ions number density -4, number of ions
# in loop add ion mass density / ion mass
# then sum this * electron mass
#==============================================================================
def velocity_dt_ion(dt,A,rho,rho_1,dFdr,dPdr,m_i,E,e_charge,ag,dMdt,u,S,ac):
#returns: ion velocity
#requires: time step, area vector, mass density previous step, mass density
#           updated step, dFdr = d/dr(A rho u^2), dPdr- pressure diff
#           ions mass, electric field, electron charge, grav acceleration
#           momentum exchange rate, velocity, mass production, centrifugal 
#           accelleration
#ref: Gombozi+1985, Glocer+ 2007 etc.
    return ((dt/(A*rho_1))*(-dFdr - (A*dPdr) +(A*rho*((e_charge/m_i)*E - ag + ac)) + \
           A*dMdt + A*u*S)+((rho*u)/rho_1))
        # ((dt/ cross section areax rho_1? )) * (-dFdr - (area*dPdr) + (area* mass density*((electron charge/ion mass)* E field
        # - gravitational acceleration  + centrifugal acceleration)) + area x dMdt + area x velocity x mass production) +
        # ((mass density x velocity/ rho_1
        # CHECK THIS
#==============================================================================
def velocity_dt_electron(electrons,ions,num_ionic_species,j,cd,e_charge):
#returns: electron velocity
#requires:electron dictionary, ion dictionary, number of ionic species,
#         iteration, current density, electron charge
#ref: Gombozi+1985, Glocer+ 2007 etc.
    import numpy as np
    term1 = np.empty([len(ions[1]["n"])-4,num_ionic_species])
    for i in range(1,num_ionic_species+1):
        term1[:,i-1] = ions[i]["n"][2:-2,j]*ions[i]["u"][2:-2,j]
    return 1/electrons["n"][2:-2,j] * (np.sum(term1,axis=1) - cd[2:-2]/e_charge) 
# empty array made of length ion number density - 4, number of ionic species 
# in loop: number density of ions * ion velocity (3rd to 3rd to last term, axis j (to remove ghost points)
# then 1/ electron number density * sum of loop
#==============================================================================    
def pressure_dt_ion(dt,A,gamma,rho,rho_1,u,u_1,m_i, e_charge,E,ag,dEdt,dMdt,P,dEngdr,dakTdr,S,ac):   
#returns: ion pressure
#requires:time step,area vector,specific heat,mass density previous step,mass 
#         density updated step, velocity previous step, velocity updated step
#         ion mass, electron charge, electric field, gravitational acceleration
#         dEdt = energy exchange rate, momentum exchange rate, pressure 
#         previous step, dEngdr = d/dr(0.5*A*rho*u^3 - gamma/(1-gamma) *A*u*P)
#         dakTdr = heat conduction differential, mass production, centrifugal
#         acceleration
#ref: Gombozi+1985, Glocer+ 2007 etc.
    return (((dt/A)*(gamma-1))*(-dEngdr+A*rho*u*((e_charge/m_i)*E-ag+ac)+dakTdr+A*dEdt+\
     A*u*dMdt+0.5*A*u**2*S)+(0.5*rho*u**2 *(gamma-1))+P-(0.5*rho_1*u_1**2*(gamma-1)))
        # (((dt/area) x gamma -1 x (-dEngdr + areaxmass density x velocity x (electron charge / ion mass) x electric field - 
        #gravitational field + centrifugal acceleration ) + dakTdr + area x dEdt + area x velocity x dMdt + 0.5 x area x velocity^2
        # x mass production ) + (0.5x mass density x velocity^2 x (gamma -1)) + pressure - (0.5 x rho_1 x u_1^2 x gamma -1))
        # CHECK THIS
#==============================================================================     
def temperature_dt_electron(dt,gamma,m_e,k_b,A,rho,u,T,S,dTedr,dAudr,dakTedr):
#returns: electron temperature
#requires:time step,specific heat, electron mass, boltzmann constant, area, 
#         electronmass density, electron velocity, electron temperature, 
#         electron mass production, dTedr- temp differential, dAudr - d/dr(au)
#         dakTedr - electron heat conduction equation
#ref: Gombozi+1985, Glocer+ 2007 etc   
    return (dt*((gamma-1)*(m_e/(k_b*A*rho))*dakTedr-u*dTedr-(T/rho)*(S+((gamma-1)/A)*rho*dAudr))+T)
# dt *( gamma -1 x (electron mass / boltzmann x area x mass density)) x dakTedr - velocity x dTedr - (temp/mass density) x 
# (mass production +( gamma -1)/area) x mass density x dAudr) + temp
# CHECK THIS
#==============================================================================
def electron_flux2(ions,num_ionic_species,j):
#returns: electron flux
#requires:ion dictionary, number of ionic species, iteration
#ref: Gombozi+1985, Glocer+ 2007 etc
    import numpy as np
    term1 = np.empty([len(ions[1]["n"])-4,num_ionic_species])
    for i in range(1,num_ionic_species+1):
        term1[:,i-1] = ions[i]["n"][2:-2,j]*ions[i]["u"][2:-2,j]   
    return np.sum(term1,axis=1)
# empty array of length ion number density -4, number of ion species
# loop: number density ions x ion velocity
# sum up terms
#==============================================================================
def electron_flux_e(electrons,j):
#returns: electron flux
#requires:electron dictionary, iteration
    return electrons["n"][2:-2,j]*electrons["u"][2:-2,j]
# electron number density x electron velocity
#==============================================================================    
def eV2vel(eV, m):
#returns: velocity
#requires:electron volts, mass
    import numpy as np
    return np.sqrt((eV*2*1.6*10**-19)/(m))
# square root of: eV x 2 x 1.6e-19 / m
# CHECK THIS
#==============================================================================
def vel2eV(vel, m):
#returns: electron volts
#requires:velocity, mass
    return (m*vel**2)/(2*1.6*10**-19)
# m x vel^2 / 2x1.6e19
# CHECK THIS
#==============================================================================
def v2T(v,m):
#returns: temperature
#requires: velocity, mass
    return (m*v**2)/(2*1.38*10**-23)
# m x v^2 / 2 x 1.38^-23
# CHECK THIS
#==============================================================================
def T2v(T,m):
#returns: velocity
#requires:temperature, mass
    import numpy as np
    return np.sqrt((2*1.38*10**-23*T)/m)
# square root of:( 2x1.38^-23 x temp) / m
# CHECK THIS
#==============================================================================
def extrap_end(a):
#returns: two values to fill ghost points by extrapolating the function a
#requires:vector a
    return [2*a[-1]-a[-2],2*(2*a[-1]-a[-2])-a[-1]]
#==============================================================================
def extrap_start(a):
#returns: two values to fill ghost points by extrapolating the function a
#requires:vector a
    return [2*(2*a[0]-a[1])-a[0],2*a[0]-a[1]]
#==============================================================================
def current_density_ray2015():
# returns a vector of values of the upward and downward currnt region with size
# 100 taken form Ray+ 2015  
    import numpy as np   
    x = np.linspace(1,100,100)
    #1 degree of upward current from 0-0.2 microA/m2
    y1 = 0.004*x[0:50]
    #1/2 degree of constant upward current at 0.2 microA/m2
    y2 = np.ones(25)*0.2
    # 1/4 degree to downward current
    x2 = np.linspace(0,12,13)
    y3 = -0.038*x2+0.2
    # 1/4 degree return to 0
    x3 = np.linspace(0,11,12)
    y4 = 0.025*x3+y3[-1]
    j = np.concatenate([y1,y2,y3,y4],axis=0)   
    return j
#==============================================================================
# End of file
