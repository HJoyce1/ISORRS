#Code written by H. Joyce (2021)
#Module designed to replace the command line if desired
#=========================================================

from pw_1Dsinglefieldline_FAC import polar_wind_outflow as pwm

# planet, time step, iterations, distance (RJ), Lshell, FAC?, CF?, plots?, saves?, name of file
pwm('jupiter',0.01,10000,3,25,0,0,1,1,'FAC1_dusk_700k')
