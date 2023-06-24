This directory contains all scripts for the "APPLICATION TO DATA"

The main script for analysing the i^th species of flora is:
'elith_data_analysis_iterated.R'
this can be used within a wrapper function to run all species
by iterating through job = 1:29.
I sent these to a HPC using 'run_analysis_elith_iterated.pbs'

'collate elith results.R' - consolidates the individual species results and creates the plots from the manuscript