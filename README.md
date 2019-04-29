# aperdrift
Tools for Apertif driftscans

Produce drift scan schedule. Still in the development phase.  

The program reads the beam pattern from a txt file in ancillary/ and creates a drift scan schedule for a user specified calibrator and time**, and user specified number of drifts per beam***.

## Usage:
```
python3 make_imaging_sched.py [options]
```

## Notes:
Run with -h to see the list of options.

** The program will make an observing schedule for any calibrator/time combination.  If a calibrator is not actually above the observing horizon (HA = +/-6 hours), it will throw a warning but still write the csv file.

*** Number of drifts per beam is an imprecise description.  It really means the number of drifts between beam centers + 1 (counting one of the beam centers). 
