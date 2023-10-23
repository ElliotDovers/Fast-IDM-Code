## Analyse Scenario ##

home.wd <- getwd()

# initialise the result storage
dat <- NULL

# change into appropriate result folder
setwd(paste(home.wd, "raw", sep = "/"))
# get a list of all the individual result files
res.list <- list.files()[grepl("res_", list.files(), fixed = T)]
# inner loop through individual sim-by-model files
for (job in res.list) {
  # get the job number
  job.no <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(job, ".", fixed = T), function(x){x[1]})), "_", fixed = T), function(x){x[length(x)]})))
  # load in the individual simulation results
  load(job)
  # add in the job number
  res_tab$job <- job.no
  # add in the scenario
  res_tab$scen <- paste0("PO", res_tab$scen_po, "_PA", res_tab$scen_pa)
  # add to the data
  dat <- rbind(dat, res_tab)
  # ranges <- rbind(ranges, attr(res_tab, "correlation ranges"))
  rm(res_tab)
}

setwd(home.wd)

save(list = "dat", file = "collated raw results.RDATA")