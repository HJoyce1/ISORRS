"""
Created on Mon Mar  4 11:12:03 2019
This function will generate a dipole field given the planet radius, rP, mass, mP, ionosphere magnetic field, b0, l/m-shell, period (in seconds) and number of points.
Also generated are field aligned centrifugal acceleration, ac, and gravity, ag.
Assumes alignment of spin and magnetic dipoles (will change this in future version).
@author: dconstable
"""

# radius, radius
def dipolefield(rE,mE,b0,lshell,period,spin_dipole_offset,numpoints,saveflag):
    # import necessary modules
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.constants

    #Values for testing - Earth
#    rE = 6.371e6
#    b0 = 3.12e-5
#    mE = 5.9736e24
#    lshell = 7
#    numpoints = 500
#    spin_dipole_offset = 11.5
#    period = 24*3600
#    saveflag = 1

    #Values for testing - Jupiter
#    rE = 7.1492e7
#    mE = 1.898e27
#    b0 = 776.6e-6
#    lshell = 30
#    numpoints = 1001
#    spin_dipole_offset = 10
#    period = (9+(56/60))*3600

    omega = 2*np.pi/period                            #Angular frequency of planet rotation

    #Calculation of relevant positions and angles
    latE = np.rad2deg(np.arccos(np.sqrt(1/lshell)))   #Latitude of L-shell crossing in degrees
    lat_deg = np.linspace(0,latE,numpoints)           #Linear array of latitudes in degrees
    lat_rad = lat_deg*np.pi/180                       #Linear array of latitudes in radians
    r = lshell*rE*np.cos(lat_rad)**2                  #Radial position of L-shell points
    r_norm = r/rE                                     #Radial position norm. to planet radius
    x = r* np.cos(lat_rad)                            #x axis aligned with equator of dipole
    y = r* np.sin(lat_rad)                            #y axis aligned with dipole axis
    phi = np.arctan(np.gradient(y)/np.gradient(x))    #Phi is angle of the field line wrt equator

    #Calculation of dipole magnetic field
    bmag = (b0/lshell**3)*np.sqrt(1+3*np.sin(lat_rad)**2)/np.cos(lat_rad)**6
    dB = np.gradient(bmag)

    #Plot to check magnetic field
#    plt.plot(r_norm, bmag, '.'),
#    plt.xlabel('r [r$_E$]'),
#    plt.ylabel('B-Field [T]'),
#    plt.xlim([np.min(r_norm),np.max(r_norm)])

    #Calculate length of individual segments of field line, ds
    dlat_rad = np.gradient(lat_rad)                   #Gradient in points in latitude (radians)
    dr = np.gradient(r)                               #Gradient in radial position
    ds = np.zeros(np.size(r))                         #Zeros for ds
    for i in range(0,np.size(ds)):
        ds[i] = np.sqrt(dr[i]**2 + r[i]**2*dlat_rad[i]**2)

    #Calculate sum of all segment lengths to get 1-D z-axis for field line
    z = np.zeros(np.size(ds))
    for i in range(0,np.size(z)):
        z[i] = sum(ds[0:i])

    #Calculate field aligned acceleration due to gravity and centrifugal accelerations
    ag = (scipy.constants.G*mE/r**2)* np.cos(phi)
    ac = omega**2 * -r * np.cos(phi)

    if saveflag == 1:
        np.savetxt('grav_data.dat',ag)
        np.savetxt('centrifugal_data.dat',ac)
        data = np.array([bmag, dB]).T
        np.savetxt('bfield_data.dat',data)
        bbound = np.zeros([2,2])
        bbound = [data[0,:],data[-1,:]]
        np.savetxt('bbound_data.dat',bbound)
    return z,lat_deg,latE,r,bmag,dB,ag,ac

def dipolefieldlimits(rE,rmin,rmax,mE,b0,lshell,period,spin_dipole_offset,numpoints,saveflag):
    #Import some modules
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.constants

    #Values for testing - Earth
#    rE = 6.371e6
#    b0 = 3.12e-5
#    mE = 5.9736e24
#    lshell = 7
#    numpoints = 500
#    spin_dipole_offset = 11.5
#    period = 24*3600
#    saveflag = 1

    #Values for testing - Jupiter
#    rE = 7.1492e7
#    rmin = 20*rE
#    rmax = 30*rE
#    mE = 1.898e27
#    b0 = 776.6e-6 # magnetic field at 1RJ
#    lshell = 30
#    numpoints = 1001
#    spin_dipole_offset = 10
#    period = (9+(56/60))*3600

    #Values for testing - Saturn
#    rE = 5.8232e7
#    rmin = 1400000
#    rmax = 61325000
#    mE = 5.683e26
#    b0 = 21e-6 # magnetic field at 1RJ
#    lshell = 10
#    numpoints = 800
#    spin_dipole_offset = 0
#    period = (10+(42/60))*3600

    omega = 2*np.pi/period                            #Angular frequency of planet rotation

    #Calculation of relevant positions and angles
    latE = np.rad2deg(np.arccos(np.sqrt(1/lshell)))   #Latitude of L-shell crossing in degrees
    lat_deg = np.linspace(0,latE,numpoints)           #Linear array of latitudes in degrees
    lat_rad_old = lat_deg*np.pi/180                       #Linear array of latitudes in radians
    r_old = lshell*rE*np.cos(lat_rad_old)**2                  #Radial position of L-shell points
    indices = np.where((r_old > rmin) & (r_old < rmax))
    lat_rad = lat_rad_old[indices]
    lat_rad = np.linspace(np.min(lat_rad),np.max(lat_rad),numpoints)
    r = lshell*rE*np.cos(lat_rad)**2
    r_norm = r/rE                                     #Radial position norm. to planet radius
    x = r* np.cos(lat_rad)                            #x axis aligned with equator of dipole
    y = r* np.sin(lat_rad)                            #y axis aligned with dipole axis
    x_old = r_old* np.cos(lat_rad_old)                            #x axis aligned with equator of dipole
    y_old = r_old* np.sin(lat_rad_old)                            #y axis aligned with dipole axis

    xnew = np.linspace(np.min(x_old),np.max(x_old),numpoints)
    ynew = np.linspace(np.max(y_old),np.max(y_old),numpoints)
    phi = np.arctan(np.gradient(y)/np.gradient(x))    #Phi is angle of the field line wrt equator


    #Calculation of dipole magnetic field
    bmag = (b0/lshell**3)*np.sqrt(1+3*np.sin(lat_rad)**2)/np.cos(lat_rad)**6
    dB = np.gradient(bmag)

    #Plot to check magnetic field
#    plt.plot(r_norm, bmag, '.'),
#    plt.xlabel('r [r$_E$]'),
#    plt.ylabel('B-Field [T]'),
#    plt.xlim([np.min(r_norm),np.max(r_norm)])

    #Calculate length of individual segments of field line, ds
    dlat_rad = np.gradient(lat_rad)                   #Gradient in points in latitude (radians)
    dr = np.gradient(r)                               #Gradient in radial position
    ds = np.zeros(np.size(r))                         #Zeros for ds
    for i in range(0,np.size(ds)):
        ds[i] = np.sqrt(dr[i]**2 + r[i]**2*dlat_rad[i]**2)

    #Calculate sum of all segment lengths to get 1-D z-axis for field line
    z = np.zeros(np.size(ds))
    for i in range(0,np.size(z)):
        z[i] = sum(ds[0:i])

    #Calculate field aligned acceleration due to gravity and centrifugal accelerations
    ag = (scipy.constants.G*mE/r**2)* np.cos(phi)
    ac = omega**2 * -r * np.cos(phi)

    if saveflag == 1:
        np.savetxt('grav_data.dat',ag)
        np.savetxt('cent_data.dat',ac)
        data = np.array([bmag, dB]).T
        np.savetxt('bfield_data.dat',data)
        bbound = np.zeros([2,2])
        bbound = [data[0,:],data[-1,:]]
        np.savetxt('bbound_data.dat',bbound)
    return z,x,y,lat_deg,latE,r,bmag,dB,ag,ac

# radius, inner boundary, outer boundary, mass of planet, magnetic field at 1 RJ, l shell value, rotation period, dipole offset, number of grid points
def dipolefield_ISORRS(rE,rmin,rmax,mE,b0,lshell,period,spin_dipole_offset,numpoints):
    # import needed modules
    import numpy as np
    import scipy.constants

    #Values for testing - Earth
#    rE = 6.371e6
#    b0 = 3.12e-5
#    mE = 5.9736e24
#    lshell = 7
#    numpoints = 500
#    spin_dipole_offset = 11.5
#    period = 24*3600
#    saveflag = 1

    #Values for testing - Jupiter
#    rE = 7.1492e7
#    rmin = 20*rE
#    rmax = 30*rE
#    mE = 1.898e27
#    b0 = 776.6e-6 # magnetic field at 1RJ
#    lshell = 30
#    numpoints = 1001
#    spin_dipole_offset = 10
#    period = (9+(56/60))*3600

    #Values for testing - Saturn
#    rE = 5.8232e7
#    rmin = 1400000
#    rmax = 5*5.8232e7
#    mE = 5.683e26
#    b0 = 21e-6 # magnetic field at 1RJ
#    lshell = 10
#    numpoints = 800
#    spin_dipole_offset = 0
#    period = (10+(42/60))*3600

    # --------- CALCULATIONS ---------
    # generates field line based on input variables
    # generates a point for each point on the grid
    # then generates angle of field line to equator in order to calculate the
    # gravitational acceleration and centrifugal acceleration generated due to
    # this field line positioning

    omega = 2*np.pi/period                            # angular frequency of planet rotation

    # calculation of relevant positions and angles
    latE = np.rad2deg(np.arccos(np.sqrt(1/lshell)))   # latitude of L-shell crossing in degrees
    lat_deg = np.linspace(latE,0,numpoints)           # linear array of latitudes in degrees (counting downward from ~70)
    lat_rad_old = lat_deg*np.pi/180                   # linear array of latitudes in radians
    r_old = lshell*rE*np.cos(lat_rad_old)**2          # radial position of L-shell points - mapped out
    indices = np.where((r_old > rmin) & (r_old < rmax))
    lat_rad = lat_rad_old[indices]
    lat_rad = np.linspace(np.max(lat_rad),np.min(lat_rad),numpoints)
    r = lshell*rE*np.cos(lat_rad)**2
    x = r* np.cos(lat_rad)                            # x axis aligned with equator of dipole
    y = r* np.sin(lat_rad)                            # y axis aligned with dipole axis
    phi = np.arctan(np.gradient(y)/np.gradient(x))    # phi is angle of the field line wrt equator

    # calculate field aligned acceleration due to gravity and centrifugal accelerations
    ag = (scipy.constants.G*mE/r**2)* np.cos(phi)
    ac = omega**2 * r * np.cos(phi)

    # return variables needed for main module
    return phi,ag,ac
