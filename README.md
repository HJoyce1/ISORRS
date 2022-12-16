# IonoSpheric Outflow in Rapidly Rotating Systems (ISORRS) Model

Repository for the development of the model <br>
More scientfic information available in doi.org/10.1029/2019JA027728 and doi.org/10.1029/2019JA027727 <br>
This README file is focused more on the construction of and how to run the model <br>

## Get Started!
You can run the `ISORRS_1D_singlefieldline.py` as a function or from command line with arguments:
`python3 ISORRS_1D_singlefieldline_FAC.py 'jupiter' 0.01 100 3 25 0 0 1 1 'save file name'`
Check the docstring to see the argument meanings.
Same with the other versions of 'MAIN' modules.

### NOTE: some files contain arguments that are not incorporated into the command line prompt yet, such as changing from 'dusk' to 'dawn' parameters within the 'ISORRS_1Dsinglefieldline_asymmetries.py' module. Please check files to ensure they are in the correct configuration before you run them. Additionally, some older modules ask for input modules that are not the current most up to date modules, so please check what modules you need to be able to run the MAIN module you are interested in.


# Contents

## Current Model Version Files:

### ISORRS_1Dsinglefieldline.py
MAIN Module set up to run a series of iterations along a single field line at given position on planet.

### ISORRS_dipolefield.py
Module developed by Dave Constable at Lancaster University (2019) to generate the dipole field for the planet, centrifugal acceleration and graviational acceleration.

### ISORRS_planet_scale_heights.py
Module originally developed by Carley Martin at Lancaster University (2019) to contain initial conditions and functions associated with Jupiter and Saturn. Has been updated by Hannah Joyce at Lancaster University (2022) to contain updated initial conditions better representing the scale-heights of ions and neutrals in Jupiter's ionosphere. 

### ISORRS_plotting_tools_cb.py
Module to aid with plotting originally designed by Carley Martin but modified by Hannah Joyce. Colourblind friendly version of pw_plotting_tools with no grid lines and inward tick markers.

### ISORRS_equations.py
Module containing the physical calculations used in the model.



## Current Work in Progress Modules

### ISORRS_1Dsinglefileline_lax.py
Current iteration of the MAIN module working on intergrating a more robust numerical method for resolving the transport equations.

### ISORRS_lax_method.py
Lax methodology intended to be merged into the equations module at a future date if this numerical method proves to improve the current state of the model.



## Files for Specific Functions

### ISORRS_1Dsinglefieldline_asymmetries.py
MAIN Module created to iterate over different conditons in the auroral region. Currently allows for temperature, width, number density and field-aligned current strength to be varied over 'dawn' or 'dusk' conditions.

### ISORRS_planet_den_asym.py
Version of the PLANET_SCALE_HEIGHTS module specifically designed to work with the ASYMMETRIES MAIN module above. Contains initial conditions for Jupiter and Saturn.

### pw_20190502.py
MAIN Module - generic file to run overall - need to look into this more. I believe this was written so Carley could print out all
the files needed for the repository for the acknowledgements kept in the new Lancaster DOI stuff. 

### pw_equatorialFlux.py
MAIN Module which calculates the amount of flux all the way out where the magnetic field crosses the equatorial region.
Obviously with no Dungey and Vasyliunas cycle consideration.



## Older Files for Reference

### pw_plotting_tools.py
Module to aid in plotting.
 
### planet.py
Original module written by Carley Martin containing the original ntial conditions/functions and general variables for the model associated with either Jupiter or Saturn.

### pw_plotting_tools_x4.py
Module to aid in plotting with E field and acceleration terms for ions multiplied by 4 in overview plots to represent results in the Jupiter paper. Not confident with why this has been done, but this file exists just for example for now.

### Note that pw.py and dipolefield.py have been renamed ISORRS_equations.py and ISORRS_dipolefield.py but have not been updated to an extend they will not work with older MAIN modules.


## Other Information

## Licence Info
Chosen the GNU Lesser license as it's the closest to the CC-BY that Lancaster required available through github. 
This license allows free distribution, use and integration with the understanding that the user will credit the authors of this model. 
At the moment, citing doi.org/10.1029/2019JA027727 is enough credit.

## References
See individual programs.
