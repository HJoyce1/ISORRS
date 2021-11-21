# Polar Outer Planet Ionospheric Outflow Model

Repository for the development of the Polar Outer Planet Ionospheric Outflow model <br>
More scientfic information available in doi.org/10.1029/2019JA027728 and doi.org/10.1029/2019JA027727 <br>
This README file is focused more on the construction of and how to run the model <br>
(yes I know polar is spelled wrong in the repo name, its a feature now, oops) <br>
(Also absolutely open to calling it the POPI Model now, it's cute?) <br>
(Also definitely requires some updating which I will do just so it's readable and PEP8 (ish))

## Get Started!
You can run the `pw_20190502.py` as a function or from command line with arguments:
`python3 pw_20190502.py 'jupiter' 0.01 100 3 33 0 0 1 0 'anything'`
Check the docstring to see the argument meanings.
Same with the other two versions.

## Contents
### dipolefield.py
Module developed by Dave Constable at Lancaster University (2019)

### planet.py
Module contains intial conditions/functions and general variables associated with either Jupiter or Saturn.

### pw_1Dsinglefieldline.py
Module set up to run a series of iterations along a single field line at given position on planet.

### pw_20190502.py
Generic file to run overall - need to look into this more. I believe this was written so I could print out all
the files needed for the repository for the acknowledgements kept in the new Lancaster doi stuff. 

### pw_equatorialFlux.py
Module which calculates the amount of flux all the way out where the magnetic field crosses the equatorial region.
Obviously with no Dungey and Vasyliunas cycle consideration.

### pw_plotting_tools.py
Module to aid in plotting.

### pw.py
Module containing the physical calculations used in the model.

## Licence Info
Chosen the GNU Lesser license as it's the closest to the CC-BY that Lancaster required available through github. 
This license allows free distribution, use and integration with the understanding that the user will credit the authors of this model. 
At the moment, citing doi.org/10.1029/2019JA027727 is enough credit.

## References
See individual programs.
