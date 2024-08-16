# make job array script

# for single scenario

model_to_test <- c("MGCV", "MGCV FIXED", "SCAMPR", "SCAMPR FIXED", "INLA PA", "INLA PO", "INLA IDM")
sim <- 1:100

# expand out all combinations
tab <- data.frame(expand.grid(sim, model_to_test))
colnames(tab) <- c("sim", "fit_model")
tab$job <- 1:nrow(tab)

write.csv(tab, file = "job_array.csv", row.names = F)

# for all scenarios

model_to_test <- c("MGCV", "SCAMPR", "INLA")
# intercepts_pa <- c(-3.825, -2.007, 1.054) # expected 50, 200, 800 presences
# intercepts_po <- c(-6.685, -5.325, -3.9525) # expected 50, 200, 800 points
# intercepts_pa <- c(-2.007, 1.054) # expected 200, 800 presences
# intercepts_po <- c(-5.325, -3.9525) # expected 200, 800 points
bias_range <- c(15, 25, 30)
env_range <- c(15, 25, 30)
lat_range <- c(15, 25, 30)
sim <- 1:100

# expand out all combinations
# tab <- data.frame(expand.grid(intercepts_pa, intercepts_po, bias_range, env_range, lat_range, sim, model_to_test))
# colnames(tab) <- c("intercept_pa", "intercept_po", "bias_range", "env_range", "lat_range", "sim", "fit_model")
tab <- data.frame(expand.grid(bias_range, env_range, lat_range, sim, model_to_test))
colnames(tab) <- c("bias_range", "env_range", "lat_range", "sim", "fit_model")
tab$job <- 1:nrow(tab)

write.csv(tab, file = "job_array_all.csv", row.names = F)
