# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 11:57:07 2021

@author: joyceh1
"""

#import pw_20190502_hannahedit as pw_HJ
from pw_1Dsinglefieldline_march2 import polar_wind_outflow as pwm
from pw_20190502 import polar_wind_outflow as pwm2
from pw_1Dsinglefieldline_FAC import polar_wind_outflow as pwm3

# planet, time step, iterations, distance (RJ), Lshell, FAC?, CF?, plots?, saves?, name of file
#pw_HJ.polar_wind_outflow('jupiter',0.01,10000,3,33,0,0,1,1,'TestA')
#pwm('jupiter',0.01,10000,3,24,0,0,1,1,'breakpoint')
#pwm2('jupiter',0.01,10000,3,24,0,0,1,1,'pw_2019test')
pwm3('jupiter',0.01,10000,3,24,0,0,1,1,'results_test_8')

'''
NOTES:
    

'''
