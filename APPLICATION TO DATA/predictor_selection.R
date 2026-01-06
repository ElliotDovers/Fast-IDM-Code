home.wd <- getwd()

## Install required packages ###################################################
if(!require(scampr, quietly = T)){
  # scampr package can be installed from source provided in the code zip
  setwd("..")
  install.packages(paste0(getwd(), "/scampr_0.0.0.9000.tar.gz"), repos = NULL, type="source")
  library(scampr)
  setwd(home.wd)
}
if(!require(disdat, quietly = T)){
  install.packages("disdat")
  library(disdat)
}
if(!require(MASS, quietly = T)){
  install.packages("MASS")
  library(MASS)
}
if(!require(pROC, quietly = T)){
  install.packages("pROC")
  library(pROC)
}
################################################################################

# Perform the spatial k-fold cross-validation #

# TOGGLE TO DETERMINE SIMULATION/JOB NUMBER (THESE CORRESPOND TO SPECIES 1-29)
# job = 1 # run the first job for example
# determine job number from pbs script
job = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

source("get_and_format_flora_data.R")

# set up the spatial CV folds #
K <- 4

# over the full region
domain$fold <- make.spatial.folds(domain, k = K)
# over the po data
presences$fold <- make.spatial.folds(presences, rangeX = range(domain$x), rangeY = range(domain$y), k = K)
# over the background points
background$fold <- make.spatial.folds(background, rangeX = range(domain$x), rangeY = range(domain$y), k = K)
# over the presence/absence data
dat_pa$fold <- make.spatial.folds(dat_pa, rangeX = range(domain$x), rangeY = range(domain$y), k = K)

# perform some checks on the CV folds
po_folds <- NULL
for (i in species) {
  po_folds <- rbind(po_folds, data.frame(sp = i, t(as.vector(table(presences[presences$spid == i, "fold"])))))
}
pa_folds <- NULL
for (i in species) {
  tmp.pa <- dat_pa[dat_pa$spid == i, ]
  fold_sums <- NULL
  for (j in levels(tmp.pa$fold)) {
    fold_sums <- c(fold_sums, sum(tmp.pa[tmp.pa$fold == j , "occ"]))
  }
  pa_folds <- rbind(pa_folds, data.frame(sp = i, t(fold_sums)))
}
# Relative numbers
# cbind(po_folds,  po_total = apply(po_folds[ , -1], 1, sum), pa_folds, pa_total = apply(pa_folds[ , -1], 1, sum))
missing.in.fold <- NULL
for (k in 1:K) {
  missing.in.fold <- cbind(missing.in.fold, po_folds[ , k + 1] == 0 & pa_folds[ , k + 1] == 0)
}
missing.sp <- data.frame(species, missing.in.fold)
colnames(missing.sp) <- c("species", paste0("f", 1:K))
missing.sp$n_po <- as.numeric(table(presences$spid))
missing.sp$n_pa <- aggregate(dat_pa$occ, by = list(dat_pa$spid), FUN = sum)$x
missing.sp$is.missing <- apply(missing.in.fold, 1, any)
# missing.sp$n_po_ge_20 <- missing.sp$n_po < 20

# plot up the CV
# nsw.destinations <- cbind.data.frame(x = c(151.2093, 150.8931, 151.7817, 153.1139, 153.6105, 150.9293, 151.6523, 152.8975, 152.0185),
#                                      y = c(-33.8688, -34.4278, -32.9283, -30.2962, -28.6419, -31.0900, -30.5036, -31.4580, -29.0574),
#                                      name = c("Sydney", "Wollongong", "Newcastle", "Coffs Harbour", "Byron Bay", "Tamworth", "Armidale", "Port Macquarie", "Tenterfield")
# )
# plot.res <- 500
# png(filename = paste0(getwd(), "/app_cv_folds.png"), width = 5.3 * plot.res, height = 6.3 * plot.res, res = plot.res)
# par(mar = c(0, 0, 1.8, 1))
# plot(vec2im(domain$fold, domain$x, domain$y), box = F, main = "Northern NSW\nSpatially Blocked four-fold CV")
# text(nsw.destinations$x, nsw.destinations$y, labels = nsw.destinations$name, col = "black")
# dev.off()

## Model the particular species ##

# obtain the species
s = species[job]

# subset presence records of species for this species
sp_presence <- presences[presences$spid == s, ]
# add background data
dat_po <- rbind(sp_presence, background)
# add in the quadrature weights
dat_po$quad.size <- region_size / nrow(background)
dat_po$quad.size[dat_po$occ == 1] <- 0 # set quadrature weight to zero at the presence records

# subsets the PA dataset to this particular species
dat_pa <- dat_pa[dat_pa$spid == s, ]

## perform CV ##

# split the data according to folds
fold.list_pa <- split(dat_pa, dat_pa$fold)
fold.list_po <- split(dat_po, dat_po$fold)
# initialise some storage lists
pa.forward <- list()
res_pa <- list()
res_po <- list()
res_idm <- list()
pa0_preds <- list()
pa_preds <- list()
po_preds <- list()
idm0_preds <- list()
idm_preds <- list()
po_preds_lambda <- list()
pa_preds_lambda <- list()
idm_preds_lambda <- list()
for (k in 1:K) {
  # combine the training data
  tmp.dat_pa <- do.call("rbind", fold.list_pa[-k])
  tmp.dat_po <- do.call("rbind", fold.list_po[-k])
  # check that there are presences and absences in this fold
  if (length(unique(tmp.dat_pa$occ)) < 2) {
    warning(paste0("training data for fold ", k, " within job ", job, " (species = ", s, ") has only one response in the PA data."))
  }
  if (length(unique(tmp.dat_po$occ)) < 2) {
    warning(paste0("training data for fold ", k, " within job ", job, " (species = ", s, ") has only quadrature points in the PO data."))
  }
  # perform forward selection on the fixed effects
  pa.glm0 <- glm(occ ~ 1,
                 family = binomial(link = "cloglog"), data = tmp.dat_pa
  )
  pa.forward[[k]] <- stepAIC(pa.glm0, scope = as.formula(paste0("occ ~ ", paste(preds, collapse = " + "), " + ", paste(paste0("I(", preds[!preds %in% c("disturb", "soilfert")], "^2)"), collapse = " + "))),
                             trace = F, direction = "forward"
  )
}
save(pa.forward, file = paste0("predictor_selected_models/", s, "_job_", job, ".RDATA"))