# Scripts

## automatic-x.R/automatic-average.R and refit-x.R/refit-average.R

These scripts determine the noise amplitudes and exponents as a function of temperature.  There are two versions.  One version is for a noise power spectra along a single axis (automatic-x.R and refit-x.R).  The other version is for the average over three axes (automatic-average.R and refit-average.R).

The following steps assume the exponent as a function of temperature is needed.

First, the automatic-x.R script should be run.  Most of the tunable parameters are at the beginning of the file.  The low-frequency limit for fitting is set on line 88 (Pi needs to be manually set depending on how noisy the data is).  The fit parameters are output to new-fit-x.dat.

If automatic-x.R does not automatically find the correct fit range, the noise amplitude and exponents as a function of temperature will show a sudden jump.

refit-x.R will fit to the data using the range in modified-range-x.dat.  This frequency ranges in this file can be manually changed.  After running refit-x.R, the noise amplitudes and exponents as a function of temperature will be in refit-x.dat.
