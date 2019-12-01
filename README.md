# The Vaganov-Shashkin tree-ring growth model (VSM)

This repository provides the Vaganov-Shashkin tree-ring growth model (VSM) in MATLAB (Octave-compatible).  VSM, originally written in Fortran, mimics subdaily and daily resolution processes of cambial growth in trees as a function of soil moisture, air temperature, and insolation, with environmental forcing modeled as the principle of limiting factors.  Our re-implementation in a high level interpreted language, while sacrificing speed, provides opportunities to systematically evaluate model parameters, generate large ensembles of simulated tree-ring chronologies, and embed proxy system modeling within data assimilation approaches to climate reconstruction.  We also provide examples of model applications which permit process-level understanding of tree ring width variations in response to environmental variations and boundary conditions. 

## Citation

Please cite this work when using this model:

Anchukaitis, K.J., M.N. Evans, M. K. Hughes, and E. Vaganov, An interpreted language implementation of the Vaganov-Shashkin tree-ring proxy system model, submitted to *Dendrochronologia*, 2019

A preprint is available at EarthaArXiv. 

Please inform Kevin Anchukaitis if you identify any bugs.  

## Basic Applications

The code model function is `vsm`. This function runs both the Environmental and Growth blocks using daily temperature and precipitation input as well as latitude (for daylength calculations).  The function call is:

```matlab
[output] = vsm(T,P,phi,syear,eyear,parameters,varargin)
```

where `T` is the daily temperature (in ), `P` is the daily precipitation (in ), `phi` is the latitude of the tree-ring site or the latitude of the source of the daily meteorological data, `syear` is the starting year of the simulation (the first year of the daily meteorological data), `eyear` is the last year of the simulation, and `parameters` is the name of the structure containing the parameters. `varargin` allows addition variables to be passed to the model without affecting the standard, required set (the first 6 input variables). `output` is a structure containing the annual, daily, and cell-specific outputs from the simulation. 

`vsm` expects an input structure `parameters` that contains the tunable parameters that determine model behavior. Several example parameter sets are included in this distribution (`generic_parameters`,`e06_parameters`, and `a06_parameters)

## Missing Input Data

Although the code has been built to attempt to catch errors related to missing input data, you should estimate or fill missing temperature and precipitation data values prior to running `vsm`.  You may use built-in MATLAB functions like `fillmissing`, but we also include a utility function `vsm_fillmiss` that mimics the internal missing value filling in the original FORTRAN code. 

## Known Limitations

The original FORTRAN model was developed for Northern Hemisphere applications.  As a consequence the definition of a year in the model is from January 1st to December 31st, which is unlikely to be adequate for Southern Hemisphere tree-ring applications. 

## Octave Functionality

The core model function `vsm` works in Octave without modification, largely by avoiding use of `nansum` and `nanmean`. Some of the demo scripts and the helper functions therein, however, use commands not available or with a different syntax in Octave and in general Octave (testing on version 4.4.1 running on a Macbook Pro with Darwin Kernel Version 17.7.0) was found to be substantially slower than the same operations in MATLAB> 

## License

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
