# Optimising Stability of Walking Robots by Justin Z. Tang
This repository contains MATLAB code for the published paper entitled "Invariant Funnels for Underactuated Dynamic Walking Robots: New Phase Variable and Experimental Validation" by Justin Z. Tang et al.  The paper was [published in the IEEE International Conference for Robotics and Automation](https://ieeexplore.ieee.org/document/7989400).

## Software requirements
Please install the following software packages prior to running the code
- Mosek numerical solver (www.mosek.com)
- Spotless (https://github.com/spot-toolbox/spotless)

## Suggested usage
After installing the above dependencies, please run the script “RUNME.m” from the SampleCode folder.  It has been written to guide the user through the multiple scripts involved in attaining the stability funnel.

This code has been tested on Linux (Ubuntu) and Mac OS X with MATLAB 2015b.

## Acknowledgements:
- “arrow.m” is an unmodified helper function obtained from Matlab file exchange here: http://au.mathworks.com/matlabcentral/fileexchange/278-arrow-m
- “write_fcn.m” is an unmodified helper function written by Benjamin Morris.
