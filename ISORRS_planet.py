# Module to read in input arrays to model for specified planet currently
# Saturn and Jupiter
# all functions are not physical or from a model, just functions to give a 'shape'
# TODO: Actual physical functions would be a good place to start (but at time of writing 
# thre were no better, possibly better now with Juno data)

# Returns: empty array with first row being the input values in dictionaries
#          Dictionaries are electrons, ions and neutrals
#          A - area function
#          FAC - field aligned currents specific for each planetfrom Ray+ 2009,2013
#          consts- array of various constants about each planet
# Requires: pw.py function module, its - iterations, 
#          phys_constants = electron mass, proton mass, boltzmann, electron charge, specific heat
#          z- length of grid, z_ext - length of grid with ghost points
#           x- linear array of index for z - mainly for ease of developing functions
#          ghosts - ghost points
# Ref: discussed in functions

# All in SI units unless specified
# Part of the lancaster polar wind model
# Author: C.J. Martin @ Lancaster University
# 2019
#==============================================================================
def saturn(its,phys_consts,z,z_ext,x,ghosts):
    import numpy as np
    import ISORRS_equations as pw
    
    # constants for Saturn
    radius= 5.8232e7
    mass_planet = 5.683e26 
    b0=21e-6 #magnetic field at surface
    rot_period=(10+(42/60))*3600
    dipole_offset=0
    g = 10.44 #graviational acceleration
    
    # combine to read out
    consts=[radius, mass_planet,b0,rot_period,dipole_offset,g]
    
    # physical constants (read in)
    m_e = phys_consts[0]
    m_p = phys_consts[1]
    k_b = phys_consts[2]
    e_charge = phys_consts[3]
    gamma = phys_consts[4]
    
    # masses of neutral and ionic species at Saturn
    m_H2 = 2*m_p
    m_He = 4*m_p
    m_H = m_p
    m_H2O = 18* m_p
    m_H_plus = m_p
    m_H3_plus = 3*m_p

    # ghost points (read in)
    gb_z = ghosts[0]
    ge_z = ghosts[1]
    gb_x = ghosts[2]
    ge_x = ghosts[3]
    
    # prefill fixed arrays
    A = np.empty([len(z)+4,])
    n_H2 = np.empty([len(z)+4,])
    rho_H2 = np.empty([len(z)+4,])
    n_He = np.empty([len(z)+4,])
    rho_He = np.empty([len(z)+4,])
    n_H = np.empty([len(z)+4,])
    rho_H = np.empty([len(z)+4,])
    n_H2O = np.empty([len(z)+4,])
    rho_H2O = np.empty([len(z)+4,])
    S_H_plus = np.empty([len(z)+4,])
    S_H3_plus = np.empty([len(z)+4,])
    S_e = np.empty([len(z)+4,])
    T_neut = np.empty([len(z)+4,])
    P_neut = np.empty([len(z)+4,])
    # prefill iterative arrays
    rho_H_plus = np.empty([len(z)+4,its])
    rho_H3_plus = np.empty([len(z)+4,its])
    rho_e = np.empty([len(z)+4,its])
    n_H_plus = np.empty([len(z)+4,its])
    n_H3_plus = np.empty([len(z)+4,its])
    n_e = np.empty([len(z)+4,its])
    u_H_plus = np.empty([len(z)+4,its])
    u_H3_plus = np.empty([len(z)+4,its])
    u_e = np.empty([len(z)+4,its])
    P_H_plus = np.empty([len(z)+4,its])
    P_H3_plus = np.empty([len(z)+4,its])
    P_e = np.empty([len(z)+4,its])
    T_H_plus = np.empty([len(z)+4,its])
    T_H3_plus = np.empty([len(z)+4,its])
    T_e = np.empty([len(z)+4,its])
    kappa_H_plus = np.empty([len(z)+4,its])
    kappa_H3_plus = np.empty([len(z)+4,its])
    kappa_e = np.empty([len(z)+4,its])
    
    #==========================================================================
    #     INITIAL CONDITIONS
    #==========================================================================
    
    # DIPOLE A
    # cross sectional area A
    alp= (1/(z[0]**3)) *0.00000001 
    A[2:-2] = alp*(z**3) +0.0001
    A[0:2]=alp*(gb_z**3)+0.0001
    A[-2:]=alp*(ge_z**3)+0.0001
    
    # field aligned currents
#    idx = (np.abs(z - radius)).argmin() #find index of closest point to 1Rs
    FAC = 5e-10 * (A[2]/A) # 54.5-572 nAm-2 from Ray+ 2013
    
    btemp=700
    # initial H+ temperature profile and hence pressure (Log values require testing and justifying (it was a test))
    T_H_plus[2:-2,0] = 25* (x/20)**2*np.exp(-0.1*(x/20)) +btemp +200*np.log(0.001*x+1)
    T_H_plus[0,0]=btemp
    T_H_plus[1,0]=btemp
    T_H_plus[-2:,0]=25* (ge_x/20)**2*np.exp(-0.1*(ge_x/20)) +btemp +200*np.log(0.001*ge_x+1)
    
    # initial H3+ temperature profile and hence pressure
    T_H3_plus[2:-2,0] = 30* (x/15)**2*np.exp(-0.1*(x/15)) +btemp +200*np.log(0.001*x+1)
    T_H3_plus[0,0]=btemp
    T_H3_plus[1,0]=btemp
    T_H3_plus[-2:,0]=30* (ge_x/15)**2*np.exp(-0.1*(ge_x/15)) +btemp +200*np.log(0.001*ge_x+1)
    
    # initial H+ ion velocity - 1eV proton eV2vel(1,m_H_plus)
    u_H_plus[2:-2,0] = pw.T2v(T_H_plus[2:-2,0],m_H_plus)
    u_H_plus[0,0]=pw.T2v(T_H_plus[0,0],m_H_plus)
    u_H_plus[1,0]=pw.T2v(T_H_plus[1,0],m_H_plus)
    u_H_plus[-1,0]=pw.T2v(T_H_plus[-1,0],m_H_plus)
    u_H_plus[-2,0]=pw.T2v(T_H_plus[-2,0],m_H_plus)
        
    # initial H3+ ion velocity - 1eV ion
    u_H3_plus[2:-2,0] = pw.T2v(T_H3_plus[2:-2,0],m_H3_plus)
    u_H3_plus[0,0]=pw.T2v(T_H3_plus[0,0],m_H3_plus)
    u_H3_plus[1,0]=pw.T2v(T_H3_plus[1,0],m_H3_plus)
    u_H3_plus[-1,0]=pw.T2v(T_H3_plus[-1,0],m_H3_plus)
    u_H3_plus[-2,0]=pw.T2v(T_H3_plus[-2,0],m_H3_plus)
        
    #initial H+ ion density - exponential decrease place holder
    n_H_plus[2:-2,0] = 5*10**10*np.exp(-0.5*(z/z[0])) + 6*10**5
    n_H_plus[0:2,0]=5*10**10*np.exp(-0.5*(gb_z/z[0])) + 6*10**5
    n_H_plus[-2:,0]=5*10**10*np.exp(-0.5*(ge_z/z[0])) + 6*10**5
    rho_H_plus[2:-2,0] = n_H_plus[2:-2,0] * m_H_plus
    rho_H_plus[0:2,0]=n_H_plus[0:2,0] * m_H_plus
    rho_H_plus[-2:,0]=n_H_plus[-2:,0] * m_H_plus

    #initial H3+ ion density - exponential decrease place holder
    n_H3_plus[2:-2,0] = 2*10**10*np.exp(-0.4*(z/z[0])) + 10**5
    n_H3_plus[0:2,0]=2*10**10*np.exp(-0.4*(gb_z/z[0])) + 10**5
    n_H3_plus[-2:,0]=2*10**10*np.exp(-0.4*(ge_z/z[0])) + 10**5
    rho_H3_plus[2:-2,0] = n_H3_plus[2:-2,0] * m_H3_plus
    rho_H3_plus[0:2,0]=n_H3_plus[0:2,0] * m_H3_plus
    rho_H3_plus[-2:,0]=n_H3_plus[-2:,0] * m_H3_plus
     
    #initial electron density - quasi-neutrality
    n_e[2:-2,0] = n_H3_plus[2:-2,0]+n_H_plus[2:-2,0]
    n_e[0:2,0]=n_H3_plus[0:2,0]+n_H_plus[0:2,0]
    n_e[-2:,0]=n_H3_plus[-2:,0]+n_H_plus[-2:,0]
    rho_e[2:-2,0] = n_e[2:-2,0] * m_e
    rho_e[0:2,0]=n_e[0:2,0]* m_e
    rho_e[-2:,0]=n_e[-2:,0]* m_e
        
    # initial electron velocity - calcualted from ion velocities and densities - quasinauetrality
    u_e[2:-2,0] = (1/n_e[2:-2,0]) * (n_H3_plus[2:-2,0] * u_H3_plus[2:-2,0] + n_H_plus[2:-2,0]*u_H_plus[2:-2,0])
    u_e[0:2,0]=(1/n_e[0:2,0]) * (n_H3_plus[0:2,0] * u_H3_plus[0:2,0] + n_H_plus[0:2,0]*u_H_plus[0:2,0])
    u_e[-2:,0]=(1/n_e[-2:,0]) * (n_H3_plus[-2:,0] * u_H3_plus[-2:,0] + n_H_plus[-2:,0]*u_H_plus[-2:,0])
      
    #initial neutral H2 density - CONSTANT
    n_H2[2:-2] = (10**10)*10**6*np.exp(-0.3*(z/z[0])) + 1000000
    n_H2[0:2]=(10**10)*10**6*np.exp(-0.3*(gb_z/z[0])) + 1000000
    n_H2[-2:]=(10**10)*10**6*np.exp(-0.3*(ge_z/z[0])) + 1000000
    rho_H2[2:-2] = n_H2[2:-2] * m_H2
    rho_H2[0:2]=n_H2[0:2] * m_H2
    rho_H2[-2:]=n_H2[-2:] * m_H2
    
    #initial neutral H density - CONSTANT
    n_H[2:-2] = (10**8)*10**6*np.exp(-0.3*(z/z[0])) + 100000
    n_H[0:2]=(10**8)*10**6*np.exp(-0.3*(gb_z/z[0])) + 100000
    n_H[-2:]=(10**8)*10**6*np.exp(-0.3*(ge_z/z[0])) + 100000
    rho_H[2:-2] = n_H[2:-2] * m_H
    rho_H[0:2]=n_H[0:2] * m_H
    rho_H[-2:]=n_H[-2:] * m_H
    
    #initial neutral H2O density - CONSTANT
    n_H2O[2:-2] = (10**2)*10**6*np.exp(-0.3*(z/z[0])) + 100000
    n_H2O[0:2]=(10**2)*10**6*np.exp(-0.3*(gb_z/z[0])) + 100000
    n_H2O[-2:]=(10**2)*10**6*np.exp(-0.3*(ge_z/z[0])) + 100000
    rho_H2O[2:-2] = n_H2O[2:-2] * m_H2O
    rho_H2O[0:2]=n_H2O[0:2] * m_H2O
    rho_H2O[-2:]=n_H2O[-2:] * m_H2O
   
    #initial neutral He density - CONSTANT
    n_He[2:-2] = 10**5*10**6*np.exp(-0.5*(z/z[0])) + 1000000
    n_He[0:2]=10**5*10**6*np.exp(-0.5*(gb_z/z[0])) + 1000000
    n_He[-2:]=10**5*10**6*np.exp(-0.5*(ge_z/z[0])) + 1000000
    rho_He[2:-2] = n_He[2:-2] * m_He
    rho_He[0:2]=n_He[0:2] * m_He
    rho_He[-2:]=n_He[-2:] * m_He
     
    # H+ mass production rate    - CONSTANT
    S_H_plus[2:-2] =n_H_plus[2:-2,0] * 0.001 * m_H_plus *np.exp(-z/10000000)
    S_H_plus[0:2]=n_H_plus[0:2,0] * 0.001* m_H_plus *np.exp(-gb_z/10000000)
    S_H_plus[-2:]=n_H_plus[-2:,0] * 0.001* m_H_plus *np.exp(-ge_z/10000000)
   
    # H3+ mass production rate      - CONSTANT
    S_H3_plus[2:-2] = n_H3_plus[2:-2,0] * 0.001* m_H3_plus *np.exp(-z/10000000)
    S_H3_plus[0:2]=n_H3_plus[0:2,0] * 0.001* m_H3_plus *np.exp(-gb_z/10000000)
    S_H3_plus[-2:]=n_H3_plus[-2:,0] * 0.001* m_H3_plus *np.exp(-ge_z/10000000)
    S_e[2:-2] = (S_H3_plus[2:-2]/m_H3_plus + S_H_plus[2:-2]/m_H_plus) * m_e
    S_e[0:2] = (S_H3_plus[0:2]/m_H3_plus + S_H_plus[0:2]/m_H_plus) * m_e
    S_e[-2:] = (S_H3_plus[-2:]/m_H3_plus + S_H_plus[-2:]/m_H_plus) * m_e  
    
    P_H_plus[2:-2,0] = pw.plasma_pressure(n_H_plus[2:-2,0],k_b,T_H_plus[2:-2,0])
    P_H_plus[0:2,0]= pw.plasma_pressure(n_H_plus[0:2,0],k_b,T_H_plus[0:2,0])
    P_H_plus[-2:,0]= pw.plasma_pressure(n_H_plus[-2:,0],k_b,T_H_plus[-2:,0])
        
    P_H3_plus[2:-2,0] = pw.plasma_pressure(n_H3_plus[2:-2,0],k_b,T_H3_plus[2:-2,0])
    P_H3_plus[0:2,0]= pw.plasma_pressure(n_H3_plus[0:2,0],k_b,T_H3_plus[0:2,0])
    P_H3_plus[-2:,0]= pw.plasma_pressure(n_H3_plus[-2:,0],k_b,T_H3_plus[-2:,0])
       
    #Initial electron temperature profile and hence pressure 
    T_e[2:-2,0] = T_H_plus[2:-2,0] +T_H3_plus[2:-2,0]
    T_e[0:2,0]=T_H_plus[0:2,0] +T_H3_plus[0:2,0]
    T_e[-2:,0]=T_H_plus[-2:,0] +T_H3_plus[-2:,0]
    
    P_e[2:-2,0] = pw.plasma_pressure(n_e[2:-2,0],k_b,T_e[2:-2,0])
    P_e[0:2,0]=pw.plasma_pressure(n_e[0:2,0],k_b,T_e[0:2,0])
    P_e[-2:,0]=pw.plasma_pressure(n_e[-2:,0],k_b,T_e[-2:,0])
       
    #Initial neutral temperature profile and hence pressure
    T_neut[2:-2] = np.ones(np.size(x))* 400
    T_neut[0:2]=np.ones(np.size(gb_x))* 400
    T_neut[-2:]=np.ones(np.size(ge_x))* 400
    u_neut = np.zeros(np.size(T_neut))
    P_neut[2:-2] = (n_H2[2:-2]+n_He[2:-2]) * k_b * T_neut[2:-2]
    P_neut[0:2]=(n_H2[0:2]+n_He[0:2]) * k_b * T_neut[0:2]
    P_neut[-2:]=(n_H2[-2:]+n_He[-2:]) * k_b * T_neut[-2:]
          
    # heat conductivities
    kappa_H_plus[2:-2,0] = pw.heat_conductivity(T_H_plus[2:-2,0],e_charge,m_H_plus,m_p)
    kappa_H_plus[0:2,0]=pw.heat_conductivity(T_H_plus[0:2,0],e_charge,m_H_plus,m_p)
    kappa_H_plus[-2:,0]=pw.heat_conductivity(T_H_plus[-2:,0],e_charge,m_H_plus,m_p)
    
    kappa_H3_plus[2:-2,0] = pw.heat_conductivity(T_H3_plus[2:-2,0],e_charge,m_H3_plus,m_p)
    kappa_H3_plus[0:2,0]=pw.heat_conductivity(T_H3_plus[0:2,0],e_charge,m_H3_plus,m_p)
    kappa_H3_plus[-2:,0]=pw.heat_conductivity(T_H3_plus[-2:,0],e_charge,m_H3_plus,m_p)
    
    kappa_e[2:-2,0] = pw.heat_conductivity_electrons(T_e[2:-2,0],e_charge,gamma)
    kappa_e[0:2,0]=pw.heat_conductivity_electrons(T_e[0:2,0],e_charge,gamma)
    kappa_e[-2:,0]=pw.heat_conductivity_electrons(T_e[-2:,0],e_charge,gamma)
    
    #Neutral gas polarizabilities
    lambda_H2 = 0.82*10**-30 # m^3 from 0.82*10**-24 cm^3 Schunk&Nagy 2000
    lambda_He = 0.21*10**-30 # m^3 from 0.21*10**-24 cm^3 Schunk&Nagy 2000
    lambda_H = 0.67*10**-30 # m^3 from 0.21*10**-24 cm^3 Schunk&Nagy 2000
    lambda_H2O = 1.48*10**-30 # m^3 from 0.21*10**-24 cm^3 Schunk&Nagy 2000
    
    #==========================================================================
    #     FILL DICTIONARIES
    #==========================================================================

    # =============ions================
    #H_plus
    ions = {1:{"name":"H plus",
               "mass":m_H_plus,
            "T": T_H_plus,
            "u": u_H_plus,
            "n": n_H_plus,
            "rho": rho_H_plus,
            "P": P_H_plus,
            "S": S_H_plus,
            "kappa":kappa_H_plus
            },
    #H3_plus
            2:{"name":"H3 plus",
               "mass":m_H3_plus,
            "T": T_H3_plus,
            "u": u_H3_plus,
            "n": n_H3_plus,
            "rho": rho_H3_plus,
            "P": P_H3_plus,
            "S": S_H3_plus,
            "kappa":kappa_H3_plus
            }
    }
    
    #combine all ionic species into one big dictionary
 
    
    #============electrons==============
    electrons = {
            "name":"Electrons",
            "mass": m_e,
            "T": T_e,
            "u": u_e,
            "n": n_e,
            "rho": rho_e,
            "P": P_e,
            "S": S_e,
            "kappa":kappa_e
    }   
    #============neutrals==============  
    #Molecular Hydrogen
    neutrals = {1: {"name":"H2",
                    "mass":m_H2,
            "T": T_neut,
            "u": u_neut,
            "n": n_H2,
            "rho": rho_H2,
            "P":P_neut,
            "lambda":lambda_H2
            },
    #Helium
            2:{"name":"He",
               "mass":m_He,
            "T": T_neut,
            "u": u_neut,
            "n": n_He,
            "rho": rho_He,
            "P":P_neut,
            "lambda":lambda_He
            },
    #Atomic hydrogen
            3:{"name":"H",
               "mass":m_H,
            "T": T_neut,
            "u": u_neut,
            "n": n_H,
            "rho": rho_H,
            "P":P_neut,
            "lambda":lambda_H
            },
    #Water
            4:{"name":"H2O",
               "mass":m_H2O,
            "T": T_neut,
            "u": u_neut,
            "n": n_H2O,
            "rho": rho_H2O,
            "P":P_neut,
            "lambda":lambda_H2O
            }
    }
    
    
    return consts,A,FAC,ions,electrons,neutrals
#==============================================================================

# JUPITER SECTION edited by H.S. Joyce @ Lancaster University

#==============================================================================
# import variables set in 1Dsinglefieldline 
def jupiter(its,phys_consts,z,z_ext,x,ghosts,b_temp,j,n_den_H_plus,n_den_H3_plus):
    # import necessary modules
    import numpy as np
    import ISORRS_equations as pw
    
    dens = 'scales' # options are 'carley' and 'scales'
    
    # constants for Jupiter
    radius= 7.1492e7
    mass_planet = 1.898e27 
    b0=776.6e-6 # magnetic field at surface
    rot_period=(9+(56/60))*3600
    dipole_offset=10 # doesnt do anything yet
    g = 24.79 # graviational acceleration
    
    # combine to read out - need this to read into 1Dsinglefieldline
    consts=[radius,mass_planet,b0,rot_period,dipole_offset,g]
    
    # unpack physical constants from 1Dsinglefieldline
    m_e = phys_consts[0]
    m_p = phys_consts[1]
    k_b = phys_consts[2]
    e_charge = phys_consts[3]
    gamma = phys_consts[4]
    
    # masses of neutral and ionic species
    m_H2 = 2*m_p
    m_He = 4*m_p
    m_H = m_p
    m_H_plus = m_p
    m_H3_plus = 3*m_p

    # unpack ghost points from 1Dsinglefieldline
    gb_z = ghosts[0]
    ge_z = ghosts[1]
    gb_x = ghosts[2]
    ge_x = ghosts[3]
    
    # prefill fixed variable arrays
    A = np.empty([len(z)+4,])
    n_H2 = np.empty([len(z)+4,])
    rho_H2 = np.empty([len(z)+4,])
    n_He = np.empty([len(z)+4,])
    rho_He = np.empty([len(z)+4,])
    n_H = np.empty([len(z)+4,])
    rho_H = np.empty([len(z)+4,])
    S_H_plus = np.empty([len(z)+4,])
    S_H3_plus = np.empty([len(z)+4,])
    S_e = np.empty([len(z)+4,])
    T_neut = np.empty([len(z)+4,])
    P_neut = np.empty([len(z)+4,])
    # prefill iterative variable arrays
    rho_H_plus = np.empty([len(z)+4,its])
    rho_H3_plus = np.empty([len(z)+4,its])
    rho_e = np.empty([len(z)+4,its])
    n_H_plus = np.empty([len(z)+4,its])
    n_H3_plus = np.empty([len(z)+4,its])
    n_e = np.empty([len(z)+4,its])
    u_H_plus = np.empty([len(z)+4,its])
    u_H3_plus = np.empty([len(z)+4,its])
    u_e = np.empty([len(z)+4,its])
    P_H_plus = np.empty([len(z)+4,its])
    P_H3_plus = np.empty([len(z)+4,its])
    P_e = np.empty([len(z)+4,its])
    T_H_plus = np.empty([len(z)+4,its])
    T_H3_plus = np.empty([len(z)+4,its])
    T_e = np.empty([len(z)+4,its])
    kappa_H_plus = np.empty([len(z)+4,its])
    kappa_H3_plus = np.empty([len(z)+4,its])
    kappa_e = np.empty([len(z)+4,its])
    
    #==========================================================================
    #     INITIAL CONDITIONS
    #==========================================================================
    
    # DIPOLE
    # cross sectional area A
    alp= (1/(z[0]**3)) *0.00000001 #1e-8
    A[2:-2] = alp*(z**3) +0.0001
    A[0:2]=alp*(gb_z**3)+0.0001
    A[-2:]=alp*(ge_z**3)+0.0001
    
    # JRM09 Mag field A
    # p = np.array([  2.58766763e-22,  -6.62512843e-19,   7.93276617e-16, -3.33492496e-13,   9.01334060e-11,   7.10497874e-09])
    # A[2:-2] = np.polyval(p, x)
    # A[0:2]=np.polyval(p,gb_x)
    # A[-2:]=np.polyval(p,ge_x)
    #breakpoint()
    
    # field aligned currents, where j is an assigned value from 1Dsinglefileline
    FAC = j *(A[2]/A) # 1-7 microAm-1 from Ray+ 2009 # 7 for auroral, 3 for run 1
    
    # ----------ION TEMPERATURE-------------
    # initial H+ temperature profile and hence pressure 
    T_H_plus[2:-2,0] = 25* (x/20)**2*np.exp(-0.1*(x/20)) +b_temp +200*np.log(0.001*x+1) #20
    T_H_plus[0,0]=b_temp
    T_H_plus[1,0]=b_temp
    T_H_plus[-2:,0]= 25* (ge_x/20)**2*np.exp(-0.1*(ge_x/20)) +b_temp +200*np.log(0.001*ge_x+1)
    
    # initial H3+ temperature profile and hence pressure
    T_H3_plus[2:-2,0] = 30* (x/15)**2*np.exp(-0.1*(x/15)) +b_temp +200*np.log(0.001*x+1) #15
    T_H3_plus[0,0]=b_temp
    T_H3_plus[1,0]=b_temp
    T_H3_plus[-2:,0]= 30* (ge_x/15)**2*np.exp(-0.1*(ge_x/15)) +b_temp +200*np.log(0.001*ge_x+1)
    
    # ----------ION VELOCITY, DEPENDENT ON ION TEMPERATURE-------------
    # initial H+ ion velocity - 1eV proton eV2vel(1,m_H_plus)
    u_H_plus[2:-2,0] = pw.T2v(T_H_plus[2:-2,0],m_H_plus)
    u_H_plus[0,0]= pw.T2v(T_H_plus[0,0],m_H_plus)
    u_H_plus[1,0]= pw.T2v(T_H_plus[1,0],m_H_plus)
    u_H_plus[-1,0]= pw.T2v(T_H_plus[-1,0],m_H_plus)
    u_H_plus[-2,0]= pw.T2v(T_H_plus[-2,0],m_H_plus)
        
    # initial H3+ ion velocity - 1eV ion
    u_H3_plus[2:-2,0] = pw.T2v(T_H3_plus[2:-2,0],m_H3_plus)
    u_H3_plus[0,0]= pw.T2v(T_H3_plus[0,0],m_H3_plus)
    u_H3_plus[1,0]= pw.T2v(T_H3_plus[1,0],m_H3_plus)
    u_H3_plus[-1,0]= pw.T2v(T_H3_plus[-1,0],m_H3_plus)
    u_H3_plus[-2,0]= pw.T2v(T_H3_plus[-2,0],m_H3_plus)
    
    # ----------ION NUMBER DENSITY-------------
    
    # density uses formula for scale height n = n0exp(-z/H), where H = kT/mg
    # here we try to n0*exp(-z/(H_scale*(mass ratio*b_temp/background T)))
    if dens == 'scales': 
        # scale height for H+
        H_plus_scale = 5*((k_b*b_temp)/(m_H_plus*g))
        # floor validated by Constable Vlasov simulations
        n_H_plus_floor = 6*10**4 #7
        # initial H+ ion density defined by scale height
        n_H_plus[2:-2,0] = n_den_H_plus*np.exp(-(z-z[0])/H_plus_scale) + n_H_plus_floor
        n_H_plus[0:2,0] = n_den_H_plus*np.exp(-(gb_z-z[0])/H_plus_scale) +  n_H_plus_floor
        n_H_plus[-2:,0] =  n_H_plus_floor #n_den_H_plus*np.exp(ge_z/H_plus_scale) +
        # scale height for H3+
        H3_plus_scale = 5*((k_b*b_temp)/(m_H3_plus*g))
        # minimum (floor) for H3+
        n_H3_plus_floor = 10**5 #8
        # initial H3+ ion density defined by scale height
        n_H3_plus[2:-2,0] = n_den_H3_plus*np.exp(-(z-z[0])/H3_plus_scale) + n_H3_plus_floor
        n_H3_plus[0:2,0] = n_den_H3_plus*np.exp(-(gb_z-z[0])/H3_plus_scale) +  n_H3_plus_floor
        n_H3_plus[-2:,0] = n_H3_plus_floor #n_den_H3_plus*np.exp(ge_z/H3_plus_scale) +
        
    elif dens == 'carley':
        # carley number densities - H+ w/ minimum floor imposed
        n_H_plus[2:-2,0] = n_den_H_plus*np.exp(-0.5*(z/z[0])) + 6*10**4
        n_H_plus[0:2,0]=n_den_H_plus*np.exp(-0.5*(gb_z/z[0])) + 6*10**4
        n_H_plus[-2:,0]=n_den_H_plus*np.exp(-0.5*(ge_z/z[0])) + 6*10**4
        # carley number densities - H3+ w/ minimum floor imposed
        n_H3_plus[2:-2,0] = n_den_H3_plus*np.exp(-0.4*(z/z[0])) + 10**5
        n_H3_plus[0:2,0]=n_den_H3_plus*np.exp(-0.4*(gb_z/z[0])) + 10**5
        n_H3_plus[-2:,0]=n_den_H3_plus*np.exp(-0.4*(ge_z/z[0])) + 10**5
    else:
        print('Invalid selection for Number Density function')
        # old number densities
        # #H_scale = 2*z[0]
        # n_H_plus[2:-2,0] = n_den_H_plus*np.exp(-(z/(H_scale*((m_H_plus/m_H_plus)*(b_temp/3000))))) + n_H_plus_floor #2*10**9 (subauroral, run 1) #1*10**10 (auroral), #5*10**8 (nonauroral)
        # n_H_plus[0:2,0] =n_den_H_plus*np.exp(-(gb_z/(H_scale*((m_H_plus/m_H_plus)*(b_temp/3000))))) + n_H_plus_floor
        # n_H_plus[-2:,0]=n_den_H_plus*np.exp(-(ge_z/(H_scale*((m_H_plus/m_H_plus)*(b_temp/3000)))))+ n_H_plus_floor
        # old number densities
        # n_H3_plus[2:-2,0] = n_den_H3_plus*np.exp(-(z/(H_scale*((m_H_plus/m_H3_plus)*(b_temp/3000))))) + n_H3_plus_floor #1*10**10 (subauroral, run 1) #5*10**10 (auroral) #1*10**9 (nonauroral)
        # n_H3_plus[0:2,0]=n_den_H3_plus*np.exp(-(gb_z/(H_scale*((m_H_plus/m_H3_plus)*(b_temp/3000))))) + n_H3_plus_floor
        # n_H3_plus[-2:,0]=n_den_H3_plus*np.exp(-(ge_z/(H_scale*((m_H_plus/m_H3_plus)*(b_temp/3000))))) + n_H3_plus_floor
    
    # ----------ION MASS DENSITY-------------
    # mass density
    rho_H_plus[2:-2,0] = n_H_plus[2:-2,0] * m_H_plus 
    rho_H_plus[0:2,0]=n_H_plus[0:2,0] * m_H_plus
    rho_H_plus[-2:,0]=n_H_plus[-2:,0] * m_H_plus

    # mass density
    rho_H3_plus[2:-2,0] = n_H3_plus[2:-2,0] * m_H3_plus
    rho_H3_plus[0:2,0]=n_H3_plus[0:2,0] * m_H3_plus
    rho_H3_plus[-2:,0]=n_H3_plus[-2:,0] * m_H3_plus
     
    # ----------ELECTRON DENSITY, DEPENDENT ON ION NUMBER DENSITY-------------
    #initial electron density - quasi-neutrality
    n_e[2:-2,0] = n_H3_plus[2:-2,0]+n_H_plus[2:-2,0]
    n_e[0:2,0]=n_H3_plus[0:2,0]+n_H_plus[0:2,0]
    n_e[-2:,0]=n_H3_plus[-2:,0]+n_H_plus[-2:,0]
    # mass density
    rho_e[2:-2,0] = n_e[2:-2,0] * m_e 
    rho_e[0:2,0]=n_e[0:2,0]* m_e
    rho_e[-2:,0]=n_e[-2:,0]* m_e
        
    # ----------ELECTRON VELOCITY, DEPENDENT ON ION NUMBER DENSITY, VELOCITY AND ELECTRON NUMBER DENSITY-------------
    # initial electron velocity - calcualted from ion velocities and densities - quasinauetrality
    u_e[2:-2,0] = (1/n_e[2:-2,0]) * (n_H3_plus[2:-2,0] * u_H3_plus[2:-2,0] + n_H_plus[2:-2,0]*u_H_plus[2:-2,0] - FAC[2:-2]/e_charge) 
    u_e[0:2,0]=(1/n_e[0:2,0]) * (n_H3_plus[0:2,0] * u_H3_plus[0:2,0] + n_H_plus[0:2,0]*u_H_plus[0:2,0] - FAC[0:2]/e_charge)
    u_e[-2:,0]=(1/n_e[-2:,0]) * (n_H3_plus[-2:,0] * u_H3_plus[-2:,0] + n_H_plus[-2:,0]*u_H_plus[-2:,0] - FAC[-2:]/e_charge)
      
    # ----------NEUTRAL ION DENSITY-------------
    # initial neutral H2 density - CONSTANT
    if dens == 'scales':
         # initial neutral H2 density - CONSTANT
         H2_scale = 5*((k_b*b_temp)/(m_H2*g))
         n_H2[2:-2] = 10e16*np.exp(-(z-z[0])/H2_scale) + 1000000
         n_H2[0:2] =  10e16*np.exp(-(gb_z-z[0])/H2_scale) +  1000000
         n_H2[-2:] =  1000000 #10e16*np.exp(ge_z/H2_scale) +  
         # initial neutral H density - CONSTANT
         H_scale = (5*(k_b*b_temp)/(m_H*g))
         n_H[2:-2] = 10e15*np.exp(-(z-z[0])/H_scale) + 1000000
         n_H[0:2] =  10e15*np.exp(-(gb_z-z[0])/H_scale) +  1000000
         n_H[-2:] =  1000000  #10e15*np.exp(ge_z/H_scale) +
         #initial neutral He density - CONSTANT
         He_scale = 5*((k_b*b_temp)/(m_He*g))
         n_He[2:-2] = 10e14*np.exp(-(z-z[0])/He_scale) + 1000000
         n_He[0:2] =  10e14*np.exp(-(gb_z-z[0])/He_scale) +  1000000
         n_He[-2:] =  1000000 #10e14*np.exp(ge_z/He_scale) + 
         
    elif dens == 'carley':
        # carley number densities - initial neutral H2 density - CONSTANT
        n_H2[2:-2] = (10**10)*10**6*np.exp(-0.3*(z/z[0])) + 1000000
        n_H2[0:2]=(10**10)*10**6*np.exp(-0.3*(gb_z/z[0])) + 1000000
        n_H2[-2:]=(10**10)*10**6*np.exp(-0.3*(ge_z/z[0])) + 1000000
        
        # carley number densities - initial neutral H density - CONSTANT
        n_H[2:-2] = (10**9)*10**6*np.exp(-0.3*(z/z[0])) + 100000
        n_H[0:2]=(10**9)*10**6*np.exp(-0.3*(gb_z/z[0])) + 100000
        n_H[-2:]=(10**9)*10**6*np.exp(-0.3*(ge_z/z[0])) + 100000
    
        # carley number densities - initial neutral He density - CONSTANT
        n_He[2:-2] = (10**8)*10**6*np.exp(-0.5*(z/z[0])) + 1000000
        n_He[0:2]=(10**8)*10**6*np.exp(-0.5*(gb_z/z[0])) + 1000000
        n_He[-2:]=(10**8)*10**6*np.exp(-0.5*(ge_z/z[0])) + 1000000
        
    else:
        print('Invalid selection for Number Density function')

        # old number densities
        # n_H2[2:-2] = 10e16 *np.exp(-(z/(H_scale*((m_H_plus/m_H2)*(b_temp/3000))))) + 1000000
        # n_H2[0:2]=10e16*np.exp(-(gb_z/(H_scale*((m_H_plus/m_H2)*(b_temp/3000))))) + 1000000
        # n_H2[-2:]=10e16*np.exp(-(ge_z/(H_scale*((m_H_plus/m_H2)*(b_temp/3000))))) + 1000000 
        # old number densities
        # n_H[2:-2] = 10e15*np.exp(-(z/(H_scale*((m_H_plus/m_H)*(b_temp/3000))))) + 1000000
        # n_H[0:2]=10e15*np.exp(-(gb_z/(H_scale*((m_H_plus/m_H)*(b_temp/3000))))) + 1000000
        # n_H[-2:]=10e15*np.exp(-(ge_z/(H_scale*((m_H_plus/m_H)*(b_temp/3000))))) + 1000000
        #old number densities
        # n_He[2:-2] = 10e14*np.exp(-(z/(H_scale*((m_H_plus/m_He)*(b_temp/3000))))) + 1000000 
        # n_He[0:2]=10e14*np.exp(-(gb_z/(H_scale*((m_H_plus/m_He)*(b_temp/3000))))) + 1000000
        # n_He[-2:]=10e14*np.exp(-(ge_z/(H_scale*((m_H_plus/m_He)*(b_temp/3000))))) + 1000000
        
    # ----------NEUTRAL MASS DENSITY-------------
    # mass density
    rho_H2[2:-2] = n_H2[2:-2] * m_H2
    rho_H2[0:2]=n_H2[0:2] * m_H2
    rho_H2[-2:]=n_H2[-2:] * m_H2 
   
    # mass density
    rho_H[2:-2] = n_H[2:-2] * m_H
    rho_H[0:2]=n_H[0:2] * m_H
    rho_H[-2:]=n_H[-2:] * m_H
   
    # mass density
    rho_He[2:-2] = n_He[2:-2] * m_He
    rho_He[0:2]=n_He[0:2] * m_He
    rho_He[-2:]=n_He[-2:] * m_He
     
    # ----------ION MASS PRODUCTION RATE-------------
    # H+ mass production rate    - CONSTANT
    S_H_plus[2:-2] =n_H_plus[2:-2,0] * 0.001 * m_H_plus *np.exp(-z/10000000) 
    S_H_plus[0:2]=n_H_plus[0:2,0] * 0.001* m_H_plus *np.exp(-gb_z/10000000)
    S_H_plus[-2:]=n_H_plus[-2:,0] * 0.001* m_H_plus *np.exp(-ge_z/10000000)
   
    # H3+ mass production rate      - CONSTANT
    S_H3_plus[2:-2] = n_H3_plus[2:-2,0] * 0.001* m_H3_plus *np.exp(-z/10000000) 
    S_H3_plus[0:2]=n_H3_plus[0:2,0] * 0.001* m_H3_plus *np.exp(-gb_z/10000000)
    S_H3_plus[-2:]=n_H3_plus[-2:,0] * 0.001* m_H3_plus *np.exp(-ge_z/10000000)
    
    # electron mass production rate    - CONSTANT
    S_e[2:-2] = (S_H3_plus[2:-2]/m_H3_plus + S_H_plus[2:-2]/m_H_plus) * m_e
    S_e[0:2] = (S_H3_plus[0:2]/m_H3_plus + S_H_plus[0:2]/m_H_plus) * m_e
    S_e[-2:] = (S_H3_plus[-2:]/m_H3_plus + S_H_plus[-2:]/m_H_plus) * m_e  

    # ----------ION PRESSURE----------
    # H+ pressure - calculated using number density and temperature
    P_H_plus[2:-2,0] = pw.plasma_pressure(n_H_plus[2:-2,0],k_b,T_H_plus[2:-2,0]) 
    P_H_plus[0:2,0]= pw.plasma_pressure(n_H_plus[0:2,0],k_b,T_H_plus[0:2,0])
    P_H_plus[-2:,0]= pw.plasma_pressure(n_H_plus[-2:,0],k_b,T_H_plus[-2:,0])
        
    # H3+ pressure - calculated using number density and temperature 
    P_H3_plus[2:-2,0] = pw.plasma_pressure(n_H3_plus[2:-2,0],k_b,T_H3_plus[2:-2,0]) 
    P_H3_plus[0:2,0]= pw.plasma_pressure(n_H3_plus[0:2,0],k_b,T_H3_plus[0:2,0])
    P_H3_plus[-2:,0]= pw.plasma_pressure(n_H3_plus[-2:,0],k_b,T_H3_plus[-2:,0])
       
    # ----------ELECTRON TEMPERAURE, DEPENDENT ON ION TEMPERATURES----------
    # initial electron temperature profile and hence pressure 
    T_e[2:-2,0] = T_H_plus[2:-2,0] +T_H3_plus[2:-2,0] 
    T_e[0:2,0]=T_H_plus[0:2,0] +T_H3_plus[0:2,0]
    T_e[-2:,0]=T_H_plus[-2:,0] +T_H3_plus[-2:,0]
    
    # ----------ELECTRON PRESSURE,DEPENDENT ON ELETRON TEMPERATURE----------
    # electron pressure
    P_e[2:-2,0] = pw.plasma_pressure(n_e[2:-2,0],k_b,T_e[2:-2,0]) 
    P_e[0:2,0]=pw.plasma_pressure(n_e[0:2,0],k_b,T_e[0:2,0])
    P_e[-2:,0]=pw.plasma_pressure(n_e[-2:,0],k_b,T_e[-2:,0])
       
    # ----------NEUTRAL TEMPERATURE----------
    # initial neutral temperature profile
    T_neut[2:-2] = np.ones(np.size(x))* b_temp 
    T_neut[0:2]=np.ones(np.size(gb_x))* b_temp
    T_neut[-2:]=np.ones(np.size(ge_x))* b_temp
    # ----------NEUTRAL VELOCITY, DEPENDENT ON NEUTRAL TEMPERATURE----------
    # initial neutral velocity
    u_neut = np.zeros(np.size(T_neut))
    # ----------NEUTRAL PRESSURE, DEPENDENT ON NEUTRAL TEMPERATURE, NEUTRAL NUMBER DENSITY----------
    # initial neutral pressure - using H2, He and H
    P_neut[2:-2] = (n_H2[2:-2]+n_He[2:-2]+n_H[2:-2]) * k_b * T_neut[2:-2]
    P_neut[0:2]=(n_H2[0:2]+n_He[0:2]+n_H[0:2]) * k_b * T_neut[0:2]
    P_neut[-2:]=(n_H2[-2:]+n_He[-2:]+n_H[-2:]) * k_b * T_neut[-2:]
          
    
    # ----------KAPPA, DEPENDENT ON TEMPERATURE----------
    # heat conductivity for H+
    kappa_H_plus[2:-2,0] = pw.heat_conductivity(T_H_plus[2:-2,0],e_charge,m_H_plus,m_p)
    kappa_H_plus[0:2,0]=pw.heat_conductivity(T_H_plus[0:2,0],e_charge,m_H_plus,m_p)
    kappa_H_plus[-2:,0]=pw.heat_conductivity(T_H_plus[-2:,0],e_charge,m_H_plus,m_p)
    
    # heat conductivty for H3+
    kappa_H3_plus[2:-2,0] = pw.heat_conductivity(T_H3_plus[2:-2,0],e_charge,m_H3_plus,m_p)
    kappa_H3_plus[0:2,0]=pw.heat_conductivity(T_H3_plus[0:2,0],e_charge,m_H3_plus,m_p)
    kappa_H3_plus[-2:,0]=pw.heat_conductivity(T_H3_plus[-2:,0],e_charge,m_H3_plus,m_p)
    
    # heat conductivity for electrons
    kappa_e[2:-2,0] = pw.heat_conductivity_electrons(T_e[2:-2,0],e_charge,gamma)
    kappa_e[0:2,0]=pw.heat_conductivity_electrons(T_e[0:2,0],e_charge,gamma)
    kappa_e[-2:,0]=pw.heat_conductivity_electrons(T_e[-2:,0],e_charge,gamma)
    
    # ----------------LAMBDA FOR NEUTRALS-----------
    # neutral gas polarizabilities
    lambda_H2 = 0.82*10**-30 # m^3 from 0.82*10**-24 cm^3 Schunk&Nagy 2000
    lambda_He = 0.21*10**-30 # m^3 from 0.21*10**-24 cm^3 Schunk&Nagy 2000
    lambda_H = 0.67*10**-30 # m^3 from 0.21*10**-24 cm^3 Schunk&Nagy 2000
    
    #==========================================================================
    #     FILL DICTIONARIES
    #==========================================================================
    
    # ============= ions ================
    # H_plus
    ions = {1:{"name":"H plus",
               "mass":m_H_plus,
            "T": T_H_plus,
            "u": u_H_plus,
            "n": n_H_plus,
            "rho": rho_H_plus,
            "P": P_H_plus,
            "S": S_H_plus,
            "kappa":kappa_H_plus
            },
    # H3_plus
            2:{"name":"H3 plus",
               "mass":m_H3_plus,
            "T": T_H3_plus,
            "u": u_H3_plus,
            "n": n_H3_plus,
            "rho": rho_H3_plus,
            "P": P_H3_plus,
            "S": S_H3_plus,
            "kappa":kappa_H3_plus
            }
    }
    
    # combine all ionic species into one big dictionary
 
    
    #============ electrons ==============
    electrons = {
            "name":"Electrons",
            "mass": m_e,
            "T": T_e,
            "u": u_e,
            "n": n_e,
            "rho": rho_e,
            "P": P_e,
            "S": S_e,
            "kappa":kappa_e
    }   
    #============neutrals==============  
    # molecular Hydrogen
    neutrals = {1: {"name":"H2",
                    "mass":m_H2,
            "T": T_neut,
            "u": u_neut,
            "n": n_H2,
            "rho": rho_H2,
            "P":P_neut,
            "lambda":lambda_H2
            },
    # helium
            2:{"name":"He",
               "mass":m_He,
            "T": T_neut,
            "u": u_neut,
            "n": n_He,
            "rho": rho_He,
            "P":P_neut,
            "lambda":lambda_He
            },
    # atomic hydrogen
            3:{"name":"H",
               "mass":m_H,
            "T": T_neut,
            "u": u_neut,
            "n": n_H,
            "rho": rho_H,
            "P":P_neut,
            "lambda":lambda_H
            }
    }
    
    # return variables needed for main module
    return consts,A,ions,electrons,neutrals

# END OF MODULE
