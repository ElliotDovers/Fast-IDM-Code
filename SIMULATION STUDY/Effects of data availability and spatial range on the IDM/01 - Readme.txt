This directory contains all scripts for the simulation section: "Effects of data availability and spatial range on the IDM"

The main script for analysing the i^th simulation is:
'simulation_analysis_iterated.R'
this can be used within a wrapper function to run all:
100 simulations x 10 PO biasing predictor ranges x 2 env. predictor ranges x 2 shared latent field ranges x 3 PO data availability scenarios x 3 PA data availability scenarios
(see Table 1 in the manuscript)
by iterating through job = 1:18000.
I sent these to a HPC using 'run_simulation_analysis_iterated.pbs'

SIMULATION:
'sim_PO_PA_data.R' - simulates the presence-only and presence/absence datasets, as well as a regular grid over the domain representing truth.

PLOTTING AND RESULT CONSOLIDATION:
'basis_search_plots.R' - produces Figure 2 in the manuscript (demo basis function search)
'sim_setup_plots.R' - produces the plots required for Figure 1 in the manuscript

within the Results folder:

'collate raw results.R' - collates the raw results saved in /Results/raw and produces the Figures included in the manuscript.
'collated raw results.RDATA' - are the collated raw results, as I have excluded the 18000 individual results files.