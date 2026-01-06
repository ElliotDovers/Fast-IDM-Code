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
job = 2 # run the first job for example
# determine job number from pbs script
# job = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

# get the data in the appropriate format
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

if (file.exists(paste0(getwd(), "/predictor_selected_models/", s, "_job_", job, ".RDATA"))) {
  load(paste0(getwd(), "/predictor_selected_models/", s, "_job_", job, ".RDATA"))
} else {
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
  save(pa.forward, file = paste0(getwd(), "/predictor_selected_models/", s, "_job_", job, ".RDATA"))
}

# loop through folds to fit the models
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
  
  if (fit_type == "scampr") {
    # fit the base scampr models
    base_po <- scampr(pa.forward[[k]]$formula, data = tmp.dat_po, include.sre = F, model.type = "PO", sre.approx = "laplace")
    base_pa <- scampr(pa.forward[[k]]$formula, data = tmp.dat_pa, include.sre = F, model.type = "PA", sre.approx = "laplace")
    base_idm <- scampr(pa.forward[[k]]$formula, data = tmp.dat_po, bias.formula = ~ 1, pa.data = tmp.dat_pa, include.sre = F, model.type = "IDM", sre.approx = "laplace", latent.po.biasing = F)
    
    # fit the PA model
    res_pa[[k]] <- basis.search.pa(base_pa, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10)
    # fit the PO model
    res_po[[k]] <- basis.search.po(base_po, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10)
    
    # fit the IDM
    if (is.null(res_pa[[k]]$basis.functions) & is.null(res_po[[k]]$basis.functions)) {
      res_idm[[k]] <- base_idm
    } else if (!is.null(res_pa[[k]]$basis.functions) & is.null(res_po[[k]]$basis.functions)) {
      res_idm[[k]] <- do.call("update", list(object = base_idm, include.sre = T, basis.functions = res_pa[[k]]$basis.functions, latent.po.biasing = F))
    } else if (is.null(res_pa[[k]]$basis.functions) & !is.null(res_po[[k]]$basis.functions)) {
      res_idm[[k]] <- do.call("update", list(object = base_idm, include.sre = T, basis.functions = res_po[[k]]$basis.functions, po.biasing.basis.functions = res_po[[k]]$basis.functions))
    } else {
      res_idm[[k]] <- do.call("update", list(object = base_idm, include.sre = T, basis.functions = res_pa[[k]]$basis.functions, po.biasing.basis.functions = res_po[[k]]$basis.functions))
    }
    
    # make predictions for each model on the test data
    pa0_preds[[k]] <- predict(pa.forward[[k]], newdata = fold.list_pa[[k]])
    pa_preds[[k]] <- predict(res_pa[[k]], newdata = fold.list_pa[[k]])
    po_preds[[k]] <- predict(res_po[[k]], newdata = fold.list_pa[[k]])
    idm_preds[[k]] <- predict(res_idm[[k]], newdata = fold.list_pa[[k]])
    idm0_preds[[k]] <- predict(base_idm, newdata = fold.list_pa[[k]])
    po_preds_lambda[[k]] <- predict(res_po[[k]], newdata = fold.list_po[[k]], include.bias.accounting = T)
    pa_preds_lambda[[k]] <- predict(res_pa[[k]], newdata = fold.list_po[[k]], include.bias.accounting = T)
    idm_preds_lambda[[k]] <- predict(res_idm[[k]], newdata = fold.list_po[[k]], include.bias.accounting = T)
    
  } else if (fit_type == "mgcv") {
    
  } else if (fit_type == "inla") {
    
  } else {
    stop("No recognised method for fitting")
  }

  # fit the base scampr models
  base_po <- scampr(pa.forward[[k]]$formula, data = tmp.dat_po, include.sre = F, model.type = "PO", sre.approx = "laplace")
  base_pa <- scampr(pa.forward[[k]]$formula, data = tmp.dat_pa, include.sre = F, model.type = "PA", sre.approx = "laplace")
  base_idm <- scampr(pa.forward[[k]]$formula, data = tmp.dat_po, bias.formula = ~ 1, pa.data = tmp.dat_pa, include.sre = F, model.type = "IDM", sre.approx = "laplace", latent.po.biasing = F)
  
   # fit the PA model
  res_pa[[k]] <- basis.search.pa(base_pa, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10)
  # fit the PO model
  res_po[[k]] <- basis.search.po(base_po, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10)
  # if (length(unique(tmp.dat_po$occ)) < 2) {
  #   res_po[[k]] <- base_po
  # } else {
  #   
  # }
  # fit the IDM
  if (is.null(res_pa[[k]]$basis.functions) & is.null(res_po[[k]]$basis.functions)) {
    res_idm[[k]] <- base_idm
  } else if (!is.null(res_pa[[k]]$basis.functions) & is.null(res_po[[k]]$basis.functions)) {
    # res_idm[[k]] <- do.call("update", list(object = base_idm, include.sre = T, basis.functions = res_pa[[k]]$basis.functions, po.biasing.basis.functions = res_pa[[k]]$basis.functions))
    res_idm[[k]] <- do.call("update", list(object = base_idm, include.sre = T, basis.functions = res_pa[[k]]$basis.functions, latent.po.biasing = F))
  } else if (is.null(res_pa[[k]]$basis.functions) & !is.null(res_po[[k]]$basis.functions)) {
    res_idm[[k]] <- do.call("update", list(object = base_idm, include.sre = T, basis.functions = res_po[[k]]$basis.functions, po.biasing.basis.functions = res_po[[k]]$basis.functions))
  } else {
    res_idm[[k]] <- do.call("update", list(object = base_idm, include.sre = T, basis.functions = res_pa[[k]]$basis.functions, po.biasing.basis.functions = res_po[[k]]$basis.functions))
  }
  # res_idm[[k]] <- do.call("basis.search.idm", list(object = base_idm, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10))
  
  # make predictions for each model on the test data
  pa0_preds[[k]] <- predict(pa.forward[[k]], newdata = fold.list_pa[[k]])
  pa_preds[[k]] <- predict(res_pa[[k]], newdata = fold.list_pa[[k]])
  po_preds[[k]] <- predict(res_po[[k]], newdata = fold.list_pa[[k]])
  idm_preds[[k]] <- predict(res_idm[[k]], newdata = fold.list_pa[[k]])
  idm0_preds[[k]] <- predict(base_idm, newdata = fold.list_pa[[k]])
  po_preds_lambda[[k]] <- predict(res_po[[k]], newdata = fold.list_po[[k]], include.bias.accounting = T)
  pa_preds_lambda[[k]] <- predict(res_pa[[k]], newdata = fold.list_po[[k]], include.bias.accounting = T)
  idm_preds_lambda[[k]] <- predict(res_idm[[k]], newdata = fold.list_po[[k]], include.bias.accounting = T)
}
# create a temporary df to store the test data responses in same order as predictions
tmp <- data.frame(ys = do.call("c", lapply(fold.list_pa, function(x){x$occ})))
# combine the predictions for each model into a vector
tmp$pres.prob_pa <- do.call("c", lapply(pa_preds, function(x){1-exp(-exp(x))}))
tmp$pres.prob_po <- do.call("c", lapply(po_preds, function(x){1-exp(-exp(x))}))
tmp$pres.prob_idm <- do.call("c", lapply(idm_preds, function(x){1-exp(-exp(x))}))

# create a temporary df to store the test data responses in same order as predictions
tmp_po <- data.frame(ys = do.call("c", lapply(fold.list_po, function(x){x$occ})))
# combine the predictions for each model into a vector
tmp_po$pres.prob_pa <- do.call("c", lapply(pa_preds_lambda, function(x){1-exp(-exp(x))}))
tmp_po$pres.prob_po <- do.call("c", lapply(po_preds_lambda, function(x){1-exp(-exp(x))}))
tmp_po$pres.prob_idm <- do.call("c", lapply(idm_preds_lambda, function(x){1-exp(-exp(x))}))

# calculate the OOS AUC
res <- data.frame(model = c("pa", "po", "idm"), job = job, spid = s,
                  auc = c(as.numeric(auc(with(tmp, roc(ys, as.vector(pres.prob_pa), quiet = T)))),
                          as.numeric(auc(with(tmp, roc(ys, as.vector(pres.prob_po), quiet = T)))),
                          as.numeric(auc(with(tmp, roc(ys, as.vector(pres.prob_idm), quiet = T))))
                  ),
                  auc_lambda = c(as.numeric(auc(with(tmp_po, roc(ys, as.vector(pres.prob_pa), quiet = T)))),
                          as.numeric(auc(with(tmp_po, roc(ys, as.vector(pres.prob_po), quiet = T)))),
                          as.numeric(auc(with(tmp_po, roc(ys, as.vector(pres.prob_idm), quiet = T))))
                  )
)

## example fits ################################################################

source("idm_scampr.R")
source("idm_inla.R")

# set the model formula
form <- occ ~ cti + disturb + mi + rainann + raindq + rugged + soildepth +
                 soilfert + solrad + tempann + tempmin + topo + I(cti^2) +
                 I(mi^2) + I(rainann^2) + I(raindq^2) + I(rugged^2) + I(soildepth^2) +
                 I(solrad^2) + I(tempann^2) + I(tempmin^2) + I(topo^2)

scampr.time <- system.time(assign("m_scampr", idm_scampr(form = form, dat_pa = dat_pa, dat_po = dat_po, domain = dat)))
# m_scampr$timing <- scampr.time
attr(res, "scampr.time") <- scampr.time

inla.time <- system.time(assign("m_inla", idm_inla(form = form, dat_pa = dat_pa, dat_po = dat_po, domain = dat)))
# m_inla$timing <- inla.time
attr(res, "inla.time") <- inla.time

# save(list = c("res", "m_scampr", "m_inla"), file = paste0("Results/", s, "_job_", job, ".RDATA"))
save(res, file = paste0("Results/", s, "_job_", job, ".RDATA"))
