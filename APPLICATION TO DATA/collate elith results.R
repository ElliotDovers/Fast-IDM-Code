### Analysis of Elith Data Results #############################################
################################################################################
home.wd <- getwd()

library(disdat)
library(raster)
library(dplyr)
library(tidyr)

res.dir <- "Results"

# some hard-coded info
r = "NSW"
r_size <- 76.18
categoricalvars <- c("vegsys", "disturb", "soilfert")
flora_groups <- c("ot", "ou", "rt", "ru") # open-forest trees, open-forest understorey plants, rainforest trees, rainforest understorey plants

# reading presence-only and background species data for this region, one file per region:
presences <- disPo(r)
background <- disBg(r)
# reduce to just species of flora
presences <- presences[presences$group %in% flora_groups, ]

# extract names and info for all species
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

pa_list <- list()
for(s in species){
  # subset presence records of species for this species
  sp_presence <- presences[presences$spid == s, ]

  grp <- sp_presence[, "group"][1]
  dat_pa <- disEnv(r, grp)
  pa <- disPa(r, grp)
  dat_pa$occ <- pa[, s] # ED: attach the presence/absence response
  
  pa_list[[which(s == species)]] <- dat_pa
  
}

# result data frame template
nsw_res <- data.frame(species.info, m = unlist(lapply(pa_list, nrow)), n_pa = unlist(lapply(pa_list, function(x){sum(x$occ)})), n_po = as.vector(table(presences$spid)))

# get the result files
sp_res_files <- list.files(paste(home.wd, res.dir, sep = "/"))[grepl("nsw", list.files(paste(home.wd, res.dir, sep = "/")), fixed = T) & !grepl("OLD", list.files(paste(home.wd, res.dir, sep = "/")), fixed = T)]

# initialise the other result columns
nsw_res$job <- NA
nsw_res$pa <- NA
nsw_res$po <- NA
nsw_res$idm <- NA

for (sp in species) {
  # check results exist
  if (any(grepl(sp, sp_res_files, fixed = T))) {
    # load the data
    load(paste0(res.dir, "/", sp_res_files[grepl(sp, sp_res_files, fixed = T)]))
    # attach to the appropriate row of results
    nsw_res[nsw_res$spid == unique(res$spid), c("job", "pa", "po", "idm")] <- pivot_wider(res, id_cols = job, values_from = auc, names_from = model)
  }
}

# retain the raw AUCs
raw_mets <- nsw_res[ , c("pa", "po", "idm")]
nsw_res$auc_idm <- nsw_res$idm
nsw_res$auc_pa <- nsw_res$pa

# remove the PO model results as this complicates the explanation
nsw_res$po <- NULL

# go through the results to bold the best of competing models
nsw_res$pa <- round(nsw_res$pa, 3)
nsw_res$idm <- round(nsw_res$idm, 3)
for (i in 1:nrow(nsw_res)) {
  m.id <- which.max(nsw_res[i , c("pa", "idm")])
  nsw_res[i , c("pa", "idm")][ , m.id] <- paste0("textbf{", nsw_res[i , c("pa", "idm")][ , m.id], "}")
}

# calculate the differences
nsw_res$auc_diff <- nsw_res$auc_idm - nsw_res$auc_pa

# create LaTex table for appendix
library(xtable)
print(xtable(nsw_res[order(nsw_res$auc_diff, decreasing = T) , !colnames(nsw_res) %in% c("job", "auc_idm", "auc_pa", "auc_diff")], label = "append:tab:app", digits = 3,
             caption = c("Results from our spatial-fold cross-validation, including predicted AUC on the PA data survey sites. Models compared are: ``PA only'' --- binomial GLMM on the PA dataset (with GRF); IDM --- the proposed integrated distribution model. We found that AUC on the survey sites was improved by the integrating the datasets in just over half of the species analysed ($16$ out of $29$).")),
      include.rownames=F, caption.placement = "top")

# set the plot resolution
plot.res <- 500

library(ggplot2)
nsw_res$pres_rate <- (nsw_res$n_po / nsw_res$n_pa)
nsw_res$pres_rate[nsw_res$pres_rate < 0.05] <- 0.04
png(filename = paste0(getwd(), "/app_results.png"), width = 4.8 * plot.res, height = 5.7 * plot.res, res = plot.res)
ggplot(data = nsw_res, aes(x = auc_pa, y = auc_idm, color = pres_rate)) +
  geom_point(aes(size = n_pa)) + theme_classic() +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "black") +
  scale_color_gradient(low = "blue", high = "red", name = expression(frac(n[PO], n[PA])), trans = "log", breaks = c(0.05, 0.25, 1, 4), labels = c("\u2264 0.05", 0.25, 1, 4)) +
  scale_x_continuous(name = "PA only") + scale_y_continuous(name = "IDM") + 
  scale_size_continuous(name = expression(n[PA]), breaks = c(5,50,500)) +
  coord_equal() +
  ggtitle("AUC (out-of-sample)")
dev.off()

# for text
max.mets <- apply(raw_mets[ , -2], 1, which.max)
table(max.mets)
16/29
table(nsw_res$auc_idm - nsw_res$auc_pa >= 0.02)

png(filename = paste0(getwd(), "/append_app_results.png"), width = 6.2 * plot.res, height = 3.7 * plot.res, res = plot.res)
par(mfrow = c(1,4))
par(mar = c(5.1,4,2.1,0))
with(nsw_res, plot(pres_rate, auc_diff, ylab = expression(paste(Delta," AUC (IDM - PA only)")), xlab = "", log = "x", xaxt = "n", main = "A"))
axis(1, at = c(0.05, 0.25, 1, 4), labels = c("0.05", "0.25", "1", "4"))
title(xlab = expression(frac(n[PO], n[PA])), line = 3.7)
abline(h = 0, col = "red")
par(mar = c(5.1,4,2.1,0))
with(nsw_res, plot(n_pa, auc_diff, yaxt = "n", xlab = expression(n[PA]), log = "x", ylab = "", xaxt = "n", main = "B"))
axis(1, at = c(5, 20, 100, 500))
abline(h = 0, col = "red")
par(mar = c(5.1,4,2.1,0))
with(nsw_res, plot(n_po, auc_diff, yaxt = "n", xlab = expression(n[PO]), log = "x", ylab = "", xaxt = "n", main = "C"))
axis(1, at = c(2, 5, 20, 50))
abline(h = 0, col = "red")
par(mar = c(5.1,4,2.1,0))
with(nsw_res, plot(m, auc_diff, yaxt = "n", xlab = expression(m), log = "x", ylab = "", xaxt = "n", main = "D"))
axis(1, at = c(1000, 1400, 1900))
abline(h = 0, col = "red")
dev.off()

png(filename = paste0(getwd(), "/append_app_results.png"), width = 4.8 * plot.res, height = 5.7 * plot.res, res = plot.res)
par(mfrow = c(2,2))
par(mar = c(0,4.1,5.1,0))
with(nsw_res, plot(log(pres_rate), auc_pa, ylab = "PA Only: AUC", xaxt = "n", ylim = c(0.5,1)))
par(mar = c(0,0,5.1,4.1))
with(nsw_res, plot(n_pa, auc_pa, yaxt = "n", xaxt = "n", ylim = c(0.5,1)))
par(mar = c(5.1,4.1,0,0))
with(nsw_res, plot(log(pres_rate), auc_idm, ylab = "IDM: AUC", xlab = "", ylim = c(0.5,1), yaxt = "n", xaxt = "n"))
title(xlab = expression(frac(n[PO], n[PA])), line = 3.7)
axis(2, at = seq(0.5, 0.9, by = 0.1))
axis(1, at = log(c(0.01,0.1,0.5,4)), labels = c(0.01,0.1,0.5,4))
par(mar = c(5.1,0,0,4.1))
with(nsw_res, plot(n_pa, auc_idm, yaxt = "n", xlab = expression(n[PA]), ylim = c(0.5,1)))
dev.off()

# create appendix plot of the spatial CV
load("NSW.RDATA")
# create spatial folds over the full region
dat$fold <- scampr::make.spatial.folds(dat, k = 4)
# get some locations
nsw.destinations <- cbind.data.frame(x = c(151.2093, 150.8931, 151.7817, 153.1139, 153.6105, 150.9293, 151.6523, 152.8975, 152.0185),
                                     y = c(-33.8688, -34.4278, -32.9283, -30.2962, -28.6419, -31.0900, -30.5036, -31.4580, -29.0574),
                                     name = c("Sydney", "Wollongong", "Newcastle", "Coffs Harbour", "Byron Bay", "Tamworth", "Armidale", "Port Macquarie", "Tenterfield")
)
plot.res <- 500
png(filename = paste0(getwd(), "/app_cv_folds.png"), width = 5.3 * plot.res, height = 6.3 * plot.res, res = plot.res)
par(mar = c(0, 0, 1.8, 1))
plot(scampr::vec2im(dat$fold, dat$x, dat$y), box = F, main = "Northern NSW\nSpatially Blocked four-fold CV")
text(nsw.destinations$x, nsw.destinations$y, labels = nsw.destinations$name, col = "black")
dev.off()