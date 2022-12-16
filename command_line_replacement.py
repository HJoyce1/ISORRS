# replacement module for the command line

# main file with up to date commentary and working functions
from ISORRS_1Dsinglefieldline import bulk_outflow as ISORRS
#from ISORRS_1Dsinglefieldline_asymmetries import bulk_outflow as ISORRS

# # current work-in-progress file used to test new features
# from ISORRS_1Dsinglefieldline_lax_boundaries import bulk_outflow as ISORRS

# planet name, time step, iterations, distance (RJ), Lshell, FAC?, CF?, plots?, saves?, name of file
# FAC = 0 yes, 1 no, CF = 1 yes, 0 no, plots 1 yes, 0 no, saves yes 1, 0 no
ISORRS('jupiter',0.01,10000,2,10,0,0,1,1,'auroral_dawn')
