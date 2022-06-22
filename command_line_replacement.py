#Code written by H. Joyce (2021)
#Module designed to replace the command line if desired
#=========================================================
from pw_1Dsinglefieldline_asymmetries import bulk_outflow as asym_out
from pw_1Dsinglefieldline_FAC import bulk_outflow as FAC_out

# planet, time step, iterations, distance (RJ), Lshell, FAC?, CF?, plots?, saves?, name of file
# FAC = 0 yes, 1 no, CF = 1 yes, 0 no, plots 1 yes, 0 no, saves yes 1, 0 no
asym_out('jupiter',0.01,10000,2,25,0,0,1,1,'test')
