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

# some hard-coded info
r = "NSW"
r_size <- 76.18 * 1000 # in meters^2
# will additionally include ordinal factors as categoricals in NSW
categoricalvars <- c("vegsys", "disturb", "soilfert")
flora_groups <- c("ot", "ou", "rt", "ru") # open-forest trees, open-forest understorey plants, rainforest trees, rainforest understorey plants

# reading presence-only and background species data for this region, one file per region:
presences <- disPo(r)
background <- disBg(r)
# reduce to just species of flora
presences <- presences[presences$group %in% flora_groups, ]

# get the relevant predictors
preds <- disPredictors(r)
preds <- preds[-13] # remove vegsys as this factor is likely correlated with the response (as it broadly categorises the main plant types in the area)

# names for all species - from supplementary material associated with Elith el al 2020
species <- unique(presences$spid)
# species.info <- unique(presences[,c("spid", "group")])
# row.names(species.info) <- 1:nrow(species.info)
# species.info$full_name <- c("Angophora costata",
#                             "Corymbia gummifera",
#                             "Corymbia intermedia",
#                             "Eucalyptus blakelyi",
#                             "Eucalyptus carnea",
#                             "Eucalyptus fastigata",
#                             "Eucalyptus campanulata",
#                             "Eucalyptus nova-anglica",
#                             "Cassinia quinquefarina",
#                             "Lepidosperma laterale",
#                             "Glycine clandestina",
#                             "Marsdenia liisae",
#                             "Imperata cylindrica",
#                             "Poa sieberiana",
#                             "Eustrephus latifolius",
#                             "Acrotriche aggregata",
#                             "Alectryon subdentantus",
#                             "Cupaniopsis anacardioides",
#                             "Diploglottis australis",
#                             "Heritiera actinophylla",
#                             "Schizomeria ovata",
#                             "Syzygium luehmanii",
#                             "Syzygium luehmanii",
#                             "Corokia whiteana",
#                             "Cyathea leichhardtiana",
#                             "Desmodium acanthocladum",
#                             "Dicksonia antarctia",
#                             "Elatostema reticulatum",
#                             "Tasmannia purpurascens"
#                             )
# save(list = "species.info", file = "sp_info.RDATA")
load("sp_info.RDATA")

# obtain the presence/absence data
pa_list <- list()
for (grp in flora_groups) {
  pa_list[[grp]] <- cbind(disEnv(r, grp), disPa(r, grp)[-(1:4)])
}
# Checks that the above is doing what we think
# dat_pa1 <- disEnv(r, flora_groups[1])
# pa1 <- disPa(r, flora_groups[1])
# dat_pa2 <- disEnv(r, flora_groups[2])
# pa2 <- disPa(r, flora_groups[2])
# dat_pa3 <- disEnv(r, flora_groups[3])
# pa3 <- disPa(r, flora_groups[3])
# dat_pa4 <- disEnv(r, flora_groups[4])
# pa4 <- disPa(r, flora_groups[4])
# all(match(pa1$siteid, dat_pa1$siteid) == 1:nrow(dat_pa1))
# all(match(pa2$siteid, dat_pa2$siteid) == 1:nrow(dat_pa2))
# all(match(pa3$siteid, dat_pa3$siteid) == 1:nrow(dat_pa3))
# all(match(pa4$siteid, dat_pa4$siteid) == 1:nrow(dat_pa4))

# get the full domain data grid
load("nsw_grid.RDATA") # load the pre-prepared full domain data

# set up the spatial CV folds #
K <- 4

# over the full region
dat$fold <- make.spatial.folds(dat, k = K)
# over the po data
presences$fold <- make.spatial.folds(presences, rangeX = range(dat$x), rangeY = range(dat$y), k = K)
# over the background points
background$fold <- make.spatial.folds(background, rangeX = range(dat$x), rangeY = range(dat$y), k = K)
# over the presence/absence data
for (i in 1:length(pa_list)) {
  pa_list[[i]]$fold <- make.spatial.folds(pa_list[[i]], rangeX = range(dat$x), rangeY = range(dat$y), k = K)
}

# perform some checks on the CV folds
po_folds <- NULL
for (i in species) {
  po_folds <- rbind(po_folds, data.frame(sp = i, t(as.vector(table(presences[presences$spid == i, "fold"])))))
}
pa_folds <- NULL
for (i in species) {
  tmp.pa <- pa_list[[unique(presences[presences$spid == i, "group"])]]
  fold_sums <- NULL
  for (j in levels(tmp.pa$fold)) {
    fold_sums <- c(fold_sums, sum(tmp.pa[tmp.pa$fold == j , i]))
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

# plot up the CV
# nsw.destinations <- cbind.data.frame(x = c(151.2093, 150.8931, 151.7817, 153.1139, 153.6105, 150.9293, 151.6523, 152.8975, 152.0185),
#                                      y = c(-33.8688, -34.4278, -32.9283, -30.2962, -28.6419, -31.0900, -30.5036, -31.4580, -29.0574),
#                                      name = c("Sydney", "Wollongong", "Newcastle", "Coffs Harbour", "Byron Bay", "Tamworth", "Armidale", "Port Macquarie", "Tenterfield")
# )
# plot.res <- 500
# png(filename = paste0(getwd(), "/app_cv_folds.png"), width = 5.3 * plot.res, height = 6.3 * plot.res, res = plot.res)
# par(mar = c(0, 0, 1.8, 1))
# plot(vec2im(dat$fold, dat$x, dat$y), box = F, main = "Northern NSW\nSpatially Blocked four-fold CV")
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
dat_po$quad.size <- r_size / nrow(background)
dat_po$quad.size[dat_po$occ == 1] <- 0 # set quadrature weight to zero at the presence records

# identify the group of the species
grp <- sp_presence[, "group"][1]

# select out the corresponding PA dataset
dat_pa <- pa_list[[grp]]
# rename the species response variable to match the PO data (for scampr models)
dat_pa$occ <- dat_pa[ , s]

# scale the covariates (according to entire range) and
# convert categorical vars to factor in both training and testing data. We use the package forcats to ensure that the levels of the factor in the evaluation data match those in the training data, regardless of whether all levels are present in the evaluation data. 

preds.dat <- rbind(dat_po[,preds], dat_pa[,preds]) # this was previously the full NSW data but I have reduced for anonymous code submission
for(i in preds){
  if(i %in% categoricalvars){
    fac_col <- i
    if (fac_col == "soilfert") { # combining the soil fertility ratings beyond 3
      dat_po[ ,fac_col][dat_po[ ,fac_col] %in% c(4,5)] <- 3
      dat_pa[ ,fac_col][dat_pa[ ,fac_col] %in% c(4,5)] <- 3
    }
    dat_po[ ,fac_col] <- as.factor(dat_po[ ,fac_col])
    dat_pa[ ,fac_col] <- as.factor(dat_pa[ ,fac_col])
    dat_pa[ ,fac_col] <- forcats::fct_expand(dat_pa[,fac_col], levels(dat_po[,fac_col]))
    dat_pa[ ,fac_col] <- forcats::fct_relevel(dat_pa[,fac_col], levels(dat_po[,fac_col]))
  } else {
    dat_po[ , i] <- (dat_po[ , i] - mean(preds.dat[ , i], na.rm = T)) / sd(preds.dat[ , i], na.rm = T)
    dat_pa[ , i] <- (dat_pa[ , i] - mean(preds.dat[ , i], na.rm = T)) / sd(preds.dat[ , i], na.rm = T)
  }
}

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
