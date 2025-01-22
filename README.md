# The Vaganov-Shashkin tree-ring growth model (VSM)

This repository provides the Vaganov-Shashkin tree-ring growth model (VSM) in MATLAB (Octave-compatible).  VSM, originally written in Fortran, mimics subdaily and daily resolution processes of cambial growth in trees as a function of soil moisture, air temperature, and insolation, with environmental forcing modeled as the principle of limiting factors.  Our re-implementation in a high level interpreted language, while sacrificing speed, provides opportunities to systematically evaluate model parameters, generate large ensembles of simulated tree-ring chronologies, and embed proxy system modeling within data assimilation approaches to climate reconstruction.  We also provide examples of model applications which permit process-level understanding of tree ring width variations in response to environmental variations and boundary conditions. 

## Versions, updates, and bugfixes

* Version 1.1 reflects a bugfix with how the growth rate effect on transpiration was set.  In the previous version the first day of growth (`fday`) specifically when used to determine when growth started to influence transpiration did not have the index to specify the year.  The result was that as the simulation proceeded it was more likely that growth rates in the very first days of the year (if they were non-zero growth rates) influenced transpiration.  Because these growth rates tend to be small (or zero) and this bug only influenced the early part of the year, in most cases this bug did not have any significant effect on the simulation; however, in some cases this bug could cause a 1 cell difference from a corrected version that properly uses the `fday` from the specific year.  The bug was identified by [Dr. Feng Wang](https://www.researchgate.net/profile/Feng-Wang-167) and the correction and update made January 22nd, 2025. 

* Version 1.0 is the version released with the publication of Anchukaitis et al. (2020). 

## Citation

Please cite this work when using this model:

Anchukaitis, K.J., M.N. Evans, M. K. Hughes, and E. Vaganov, An interpreted language implementation of the Vaganov-Shashkin tree-ring proxy system model, *Dendrochronologia*, 60, 125677, doi:[10.1016/j.dendro.2020.125677](https://doi.org/10.1016/j.dendro.2020.125677), 2020

A earlier preprint is still available at [EartharXiv](https://doi.org/10.31223/osf.io/ruyad). 

Please inform Kevin Anchukaitis if you identify any bugs.  

## Basic Applications

The core model function is `vsm`. This function runs both the Environmental and Growth blocks using daily temperature and precipitation input as well as latitude (for daylength calculations).  The function call is:

```matlab
[output] = vsm(T,P,phi,syear,eyear,parameters,varargin)
```

where `T` is the daily temperature (in ), `P` is the daily precipitation (in ), `phi` is the latitude of the tree-ring site or the latitude of the source of the daily meteorological data, `syear` is the starting year of the simulation (the first year of the daily meteorological data), `eyear` is the last year of the simulation, and `parameters` is the name of the structure containing the parameters. `varargin` allows addition variables to be passed to the model without affecting the standard, required set (the first 6 input variables). `output` is a structure containing the annual, daily, and cell-specific outputs from the simulation. 

`vsm` expects an input structure `parameters` that contains the tunable parameters that determine model behavior. Several example parameter sets are included in this distribution (`generic_parameters`,`wm1_parameters`,`e06_parameters`, and `a06_parameters`)

## Missing Input Data

Although the code has been built to attempt to catch errors related to missing input data, you should estimate or fill missing temperature and precipitation data values prior to running `vsm`.  You may use built-in MATLAB functions like `fillmissing`, but we also include a utility function `vsm_fillmiss` that mimics the internal missing value filling in the original FORTRAN code. 

## Speed
Initial testing of the code in this repository was primarily done on a Macbook Pro (Retina, 15-inch, Mid 2015, 2.8 GHz Intel Core i7, 16 GB 1600 MHz DDR3) running MATLAB R2017b on MacOS 10.13.6 (High Sierra).  A single 22 year simulation (as for the White Mountain bristlecone example) takes approximately 0.20 seconds to complete. A ~80 year simulation (as for the Mohonk Hemlock simulation) take approximately 1.1 seconds to complete. 

## Demonstration Files

This repository comes with three demonstration files that recreate the figures in the manuscript and demonstrate the functionality of `vsm`. The first (`dendrochronologia_demo_1.m`) reproduces the White Mountain (USA) bristlecone simulations from the manuscript.  This script contains a loop over a range of drainage rate parameter value that can take up 40 second or more to complete depending on your system, but the code is annotated if you wish to reduce the number of drainage rates used for the simulation (default is 191 different rates). The second script (`dendrochronologia_demo_2.m`) reproduces the analysis of the temperature parameter sensitivity at Mohonk (USA) using a Latin Hypercube design (`lhsdesignbnd`).  The manuscript itself uses a design with 1000 draws, but this can require anywhere from 15 to 30 minutes to complete, do the file is set to start with a design of 100 members initially.  The third demonstration script uses syntax for setting the characteristics of `bar` charts introduced in MATLAB R2014B. 

## Known Limitations

The original FORTRAN model was developed for Northern Hemisphere applications.  As a consequence the definition of a year in the model is from January 1st to December 31st, which is unlikely to be adequate for Southern Hemisphere tree-ring applications. 

The original FORTRAN also had additional modules for calculating cell size and which could also be used to derive density-type metrics.  This module has not be migrated to MATLAB. 

## Octave Functionality

The core model function `vsm` works in Octave without modification (avoiding use of `nansum` and `nanmean` that are not part of core Octave distributions). Some of the demo scripts and the helper functions therein, however, use commands not available or with a different syntax in Octave and in general Octave (testing on version 4.4.1 running on a Macbook Pro with Darwin Kernel Version 17.7.0) was found to be substantially slower than the same operations in MATLAB. 

## License

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
