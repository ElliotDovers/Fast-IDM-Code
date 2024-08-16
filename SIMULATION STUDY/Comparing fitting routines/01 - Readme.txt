This directory contains all scripts for the simulation section: "Comparing scampr and INLA"

The main script for analysing the i^th simulation is:
'simulation_analysis_iterated.R'
this can be used within a wrapper function to run all 100 simulations x 10 models x 2 scenarios
by iterating through job = 1:1000.
I sent these to a HPC using 'run_simulation_analysis_iterated.pbs'

SIMULATION:
'sim_PO_PA_data.R' - simulates the presence-only and presence/absence datasets, as well as a regular grid over the domain representing truth.

MODEL FITTING AND EVALUATION FUNCTIONS:
'scampr_all.R' - fits the PA only, PO only and IDM via scampr with optimised basis functions
'scampr_fixed_all.R' - fits the PA only, PO only and IDM via scampr with a fixed set of basis functions
'inla_pa.R' - fits the PA only model via INLA
'inla_po.R' - fits the PO only model via INLA
'inla_idm.R' - fits the IDM via INLA

'collate raw results.R' - collates the raw results saved in /Results and produces the plot included in the manuscript.