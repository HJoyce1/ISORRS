# Polar Outer Planet Ionospheric Outflow Model

Repository for the development of the Polar Outer Planet Ionospheric Outflow model <br>
More scientfic information available in doi.org/10.1029/2019JA027728 and doi.org/10.1029/2019JA027727 <br>
This README file is focused more on the construction of and how to run the model <br>
(yes I know polar is spelled wrong in the repo name, its a feature now, oops) <br>
(Also absolutely open to calling it the POPI Model now, it's cute?) <br>
(Also definitely requires some updating which I will do just so it's readable and PEP8 (ish))

## Get Started!
Instructions on how to run here.

## Contents
### dipolefield.py
Module developed by Dave Constable at Lancaster University (2019)

### planet.py
Module contains intial conditions and general variables associated with either Jupiter or Saturn.

### pw_1Dsinglefieldline.py
Module set up to run a series of iterations along a single field line at given position and planet.

### pw_20190502.py
Generic file to run overall - need to look into this more.

### pw_equatorialFlux.py
Module which calculates the amount of flux all the way out where the magnetic field crosses the equatorial region.

### pw_plotting_tools.py
Module to aid in plotting.

### pw.py
Module containing the physical calculations used in the model.

## Licence Info
Chosen the GNU Lesser license as it's the closest to the CC-BY that Lancaster required available through github. 
This license allows free distribution, use and integration with the understanding that the user will credit the authors of this model. 
At the moment, citing doi.org/10.1029/2019JA027727 is enough credit.

## References
