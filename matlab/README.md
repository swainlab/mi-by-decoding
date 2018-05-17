# mi-by-decoding

Matlab scripts for estimating mutual information by decoding stochastic time series data.

## Installation

Add the `external`, `MutualInformationCode`, and `DataProcessingFunctions`
directories to your matlab path, e.g., change to this directory in Matlab and
run `addpath(genpath('.'));` at the Matlab command prompt.

In order to run the example script, you will also need to install the JSONlab
toolbox that can be found in the `external` directory. Simply open the
`jsonlab-1.5.mltbx` file in Matlab to install it.

## Getting started

The function that performs mutual information estimation is `MIdecoding`.
Type `help MIdecoding` at the Matlab command prompt for more information on
how to use it. 

Check the `ExampleScript` directory for a worked example using data from
Granados, et al. (2018) Distributed and dynamic intracellular organization of
extracellular information.

The licence for this software is described in `LICENCE.txt` located in the
parent directory.
