# IonoSpheric Outflow in Rapidly Rotating Systems (ISORRS) Model

Repository for the development of the model <br>
More scientfic information available in doi.org/10.1029/2019JA027728 and doi.org/10.1029/2019JA027727 <br>
This README file is focused more on the construction of and how to run the model <br>

## Get Started!
You can run the `ISORRS_1D_singlefieldline.py` as a function or from command line with arguments:
`python3 ISORRS_1D_singlefieldline_FAC.py 'jupiter' 0.01 100 3 25 0 0 1 1 'save file name'`
Check the docstring to see the argument meanings.
Same with the other versions of 'MAIN' modules.

## Contents
### dipolefield.py
Module developed by Dave Constable at Lancaster University (2019)

### planet.py
Module contains intial conditions/functions and general variables associated with either Jupiter or Saturn.

### ISORRS_planet_scale_heights.py
Updated version of planet with an updated number density profile to better incorporate scale heights in (currently only) Jupiter's atmosphere.

### ISORRS_1Dsinglefieldline.py
MAIN Module set up to run a series of iterations along a single field line at given position on planet.

### ISORRS_1Dsinglefieldline_asymmetries.py
MAIN Module created to iterate over different conditons in the auroral region. Currently allows for temperature, width and field-aligned current strength to be varied over 'dawn' or 'dusk' conditions.

### pw_20190502.py
MAIN Generic file to run overall - need to look into this more. I believe this was written so Carley could print out all
the files needed for the repository for the acknowledgements kept in the new Lancaster DOI stuff. 

### pw_equatorialFlux.py
MAIN Module which calculates the amount of flux all the way out where the magnetic field crosses the equatorial region.
Obviously with no Dungey and Vasyliunas cycle consideration.

### pw_plotting_tools.py
Module to aid in plotting.

### ISORRS_plotting_tools_cb.py
Colourblind friendly version of pw_plotting_tools with no grid lines and inward tick markers.

### pw_plotting_tools_x4.py
Module to aid in plotting with E field and acceleration terms for ions multiplied by 4 in overview plots to represent results in the Jupiter paper. Not confident with why this has been done, but this file exists just for example for now.

### ISORRS_equations.py
Module containing the physical calculations used in the model.

## Licence Info
Chosen the GNU Lesser license as it's the closest to the CC-BY that Lancaster required available through github. 
This license allows free distribution, use and integration with the understanding that the user will credit the authors of this model. 
At the moment, citing doi.org/10.1029/2019JA027727 is enough credit.

## References
See individual programs.
