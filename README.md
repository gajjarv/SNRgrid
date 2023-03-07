## Description

This code takes a number of filterbank files as input and measures the signal-to-noise ratio (SNR) of a spectral line in each file. It creates a grid of measured SNR as a function of Right Ascension (RA) and Declination (DEC).

## Usage


The code takes the following positional and optional arguments:

* `fil` - Path to the .fil file

Optional arguments:

* `--f_start_on` - Starting frequency for on-line
* `--f_stop_on` - Stopping frequency for on-line
* `--f_start_off` - Starting frequency for off-line
* `--f_stop_off` - Stopping frequency for off-line

## Output

The code generates a plot of SNR as a function of RA and DEC on the left and also generates three spectra lines with Max, Min, and Half Power SNR files. The plot can be saved as a .png file for later reference.

![alt text](https://github.com/gajjarv/SNRgrid/blob/main/Grid_pattern_RA_DEC_Freq_Power.png)
