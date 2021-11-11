# Functions for polar wind model
# All in SI units unless specified
# Part of the lancaster polar wind model
# Author: C.J. Martin @ Lancaster University
# 2019
#==============================================================================
def red_mass(m_i,m_j): 
#Returns: reduced mass
#Requires: mass of two species - ionic-i and neutral-j
#Ref: SCHUNK + NAGY 2000
    return (m_i*m_j)/(m_i+m_j)
#==============================================================================    
def collision_freq(rho_n,m_i,m_j,lambda_n, e_charge):
#Returns: collision frequency between two species - ionic-i and neutral-j
#Requires: neutral mass density, ionic mass, neutral mass, neutral polarisability
#            and electron charge
#Ref: SCHUNK + NAGY 2000
    import numpy as np
    return 2.21*np.pi*(rho_n/(m_i+m_j))*((lambda_n*e_charge**2)/\
                       red_mass(m_i,m_j))**0.5
#==============================================================================
def momentum_rate_1species(rho_n,m_i,m_j,lambda_n, e_charge,rho_i,u_i):
#Returns:momentum exhange rate for one ionic species and one neutral species
#Requires:neutral mass density, ionic mass, neutral mass, neutral polarisability 
#         electron charge, ionic mass density and ion speed (assuming neutral 
#         velocity is 0)
#Ref: SCHUNK + NAGY 2000
    return rho_i*collision_freq(rho_n,m_i,m_j,lambda_n, e_charge)*u_i
#==============================================================================
def energy_rate_1species(rho_i,rho_n,m_i,m_j,lambda_n, e_charge,T_j,T_i,u_i,k_b):
#Returns:energy exhange rate for one ionic species and one neutral species
#Requires:ionic mass density, neutral mass density, ionic mass, neutral mass,  
#         neutral polarisability, electron charge, neutral temp, ion temp,
#         ion speed and boltzmann constant (assuming neutral velocity is 0)
#Ref: SCHUNK + NAGY 2000
    return ((rho_i * collision_freq(rho_n,m_i,m_j,lambda_n, e_charge))/\
            (m_i + m_j)) * (3* k_b * (T_j-T_i) + m_j*u_i**2)
#==============================================================================
def heat_conductivity(T,e_charge,m_i,m_p):
#Returns: ion heat conductivity
#Requires: electron temp. electron charge, mass ion, mass neutral
#Ref: Banks&Kockarts1973, Moore+ 2008
    return (4.6*10**4 *(m_i/m_p)**-0.5 * T**2.5)*e_charge*100 #J m-1 s-1 K-1
#==============================================================================
def heat_conductivity_electrons(T,e_charge,gamma):
#Returns: electron heat conductivity
#Requires: electron temp. electron charge, specific heat
#Ref: Banks&Kockarts1973, Raitt+ 1975: 7.7 x 10^5
    return ((1.8*10**6)*(T**(5/2)))*e_charge*100 #J m-1 s-1 K-1
#==============================================================================
#Function to calculate plasma pressure
def plasma_pressure(n,k_b,T):
#Returns: plasma pressure
#Requires: number density, boltzmann constant, temperature
#Ref: P = nkT  
    return (n * k_b * T)
#==============================================================================
def plasma_temperature(n,k_b,P):
#Returns: plasma temperature
#Requires: number density, boltzmann constant, pressure
#Ref: P = nkT  
    return (P/(n*k_b))
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
#==============================================================================
def E_parallel_short(e_charge, n_e, dFdt,A, dAdr, rho_e,u_e):
#Returns: electric field, short term
#Requires: electron charge, electron number density, dFdt = d/dr(P_e - rho_e*u_e^2)
#          areaA, spatial area differential, electron mass density, electron velocity
#Ref: Gombozi+1985, Glocer+ 2007 etc.
    return (-1/(e_charge*n_e) * (dFdt + (dAdr/A)*rho_e*u_e**2))
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
    return ((dt/A)*(A*S-dFdt)+rho)
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
#==============================================================================     
def temperature_dt_electron(dt,gamma,m_e,k_b,A,rho,u,T,S,dTedr,dAudr,dakTedr):
#returns: electron temperature
#requires:time step,specific heat, electron mass, boltzmann constant, area, 
#         electronmass density, electron velocity, electron temperature, 
#         electron mass production, dTedr- temp differential, dAudr - d/dr(au)
#         dakTedr - electron heat conduction equation
#ref: Gombozi+1985, Glocer+ 2007 etc   
    return (dt*((gamma-1)*(m_e/(k_b*A*rho))*dakTedr-u*dTedr-(T/rho)*(S+((gamma-1)/A)*rho*dAudr))+T)
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
#==============================================================================
def electron_flux_e(electrons,j):
#returns: electron flux
#requires:electron dictionary, iteration
    return electrons["n"][2:-2,j]*electrons["u"][2:-2,j]
#==============================================================================    
def eV2vel(eV, m):
#returns: velocity
#requires:electron volts, mass
    import numpy as np
    return np.sqrt((eV*2*1.6*10**-19)/(m))
#==============================================================================
def vel2eV(vel, m):
#returns: electron volts
#requires:velocity, mass
    return (m*vel**2)/(2*1.6*10**-19)
#==============================================================================
def v2T(v,m):
#returns: temperature
#requires: velocity, mass
    return (m*v**2)/(2*1.38*10**-23)
#==============================================================================
def T2v(T,m):
#returns: velocity
#requires:temperature, mass
    import numpy as np
    return np.sqrt((2*1.38*10**-23*T)/m)
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
# End of file