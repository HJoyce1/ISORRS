"""
Created on Mon Nov 22 11:57:07 2021

@author: joyceh1
"""

# import pw_20190502_hannahedit as pw_HJ
import pw_1Dsinglefieldline_wipfix as pwm

# pw_20190502_hannahedit.polar_wind_outflow('jupiter',0.01,100,3,33,0,0,1,1,'TestA')

pwm.polar_wind_outflow('jupiter',0.01,50000,3,33,0,0,1,1,'1C')
# 0.01,10000,3,16,0,0,1,0

'''
notes:  
    
    
    current issues: inputs looks mostly correct now! not 100% sure if the changes made are going to fix issues though as it was
    playing with the numbers more than anything else
    results plot: electron and ion fluxes increase drastically at the point they should be flattening out? yes there is an 
    initial small increase after the minimum but then should decrease and flatten
    
    all electric fields have the wrong starting point, but the right shape? 
    
    unsure about temp and number density as have nothing to compare to
    
'''