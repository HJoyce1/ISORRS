# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 11:57:07 2021

@author: joyceh1
"""

#import pw_20190502_hannahedit as pw_HJ
from pw_1Dsinglefieldline_testing import polar_wind_outflow as pwm
#import pw_plotting_tools_editt as pwpt

#pw_HJ.polar_wind_outflow('jupiter',0.01,10000,3,33,0,0,1,1,'TestA')

pwm('jupiter',0.01,10000,3,12,1,0,1,1,'fixrhotest1')
# 0.01,10000,3,16,0,0,1,0

#pwm('planet_name', dt, its, rs, lshell, FAC_flag, CF_flag, plots, saves, run_name)

#pwpt.plot_me_quick(e_flux,zz,'xlabel','Distance','on')


# print(pwm.e_flux)

'''
notes:  
    
    
  currently n is initially a factor of 10 to low and drops off by another 10**2 by the end
  is only changed by rho, so either changes in rho are the problem or the initial values
    
'''