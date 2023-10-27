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

# add in space for results for either prediction type
nsw_res <- rbind(cbind(nsw_res, pred_type = "mu(s)"), cbind(nsw_res, pred_type = "lambda(s)")) 

timing.scampr <- NULL
timing.inla <- NULL
for (sp in species) {
  # check results exist
  if (any(grepl(sp, sp_res_files, fixed = T))) {
    # load the data
    load(paste0(res.dir, "/", sp_res_files[grepl(sp, sp_res_files, fixed = T)]))
    # attach to the appropriate row of results
    nsw_res[nsw_res$spid == unique(res$spid) & nsw_res$pred_type == "mu(s)", c("job", "pa", "po", "idm")] <- pivot_wider(res, id_cols = job, values_from = auc, names_from = model)
    nsw_res[nsw_res$spid == unique(res$spid) & nsw_res$pred_type == "lambda(s)", c("job", "pa", "po", "idm")] <- pivot_wider(res, id_cols = job, values_from = auc_lambda, names_from = model)
    timing.scampr[which(sp == species)] <- attr(res, "scampr.time")[3]
    timing.inla[which(sp == species)] <- attr(res, "inla.time")[3]
  } else {
    timing.scampr[which(sp == species)] <- NA
    timing.inla[which(sp == species)] <- NA
  }
}

# retain the raw AUCs
nsw_res$auc_idm <- nsw_res$idm
nsw_res$auc_pa <- nsw_res$pa
nsw_res$auc_po <- nsw_res$po

# go through the results to bold the best of competing models
nsw_res$pa <- round(nsw_res$pa, 3)
nsw_res$po <- round(nsw_res$po, 3)
nsw_res$idm <- round(nsw_res$idm, 3)
nsw_res_pa <- nsw_res[nsw_res$pred_type == "mu(s)", ]
nsw_res_po <- nsw_res[nsw_res$pred_type == "lambda(s)", ]
for (i in 1:nrow(nsw_res_pa)) {
  m.id <- which.max(nsw_res_pa[i , c("auc_pa", "auc_idm")])
  nsw_res_pa[i , c("pa", "idm")][ , m.id] <- paste0("textbf{", nsw_res_pa[i , c("pa", "idm")][ , m.id], "}")
  m.id <- which.max(nsw_res_po[i , c("auc_po", "auc_idm")])
  nsw_res_po[i , c("po", "idm")][ , m.id] <- paste0("textbf{", nsw_res_po[i , c("po", "idm")][ , m.id], "}")
}

# calculate the differences
# nsw_res$auc_diff <- nsw_res$auc_idm - nsw_res$auc_pa
nsw_res_pa$auc_diff <- nsw_res_pa$auc_idm - nsw_res_pa$auc_pa
nsw_res_po$auc_diff <- nsw_res_po$auc_idm - nsw_res_po$auc_po

# create LaTex table for appendix
library(xtable)
print(xtable(cbind(nsw_res_pa[order(nsw_res_pa$auc_diff, decreasing = T) , !colnames(nsw_res_pa) %in% c("job", "pred_type", "po", "auc_idm", "auc_pa", "auc_po", "auc_diff")], nsw_res_po[order(nsw_res_pa$auc_diff, decreasing = T) , !colnames(nsw_res_po) %in% c("spid", "group", "full_name", "m", "n_pa", "n_po", "job", "pred_type", "pa", "auc_idm", "auc_pa", "auc_po", "auc_diff")]), label = "append:tab:app", digits = 3,
             caption = c("Results from our spatial-fold cross-validation, including predicted AUC on the PA data survey sites. Models compared are: ``PA only'' --- binomial GLMM on the PA dataset (with GRF); IDM --- the proposed integrated distribution model. We found that AUC on the survey sites was improved by the integrating the datasets in just over half of the species analysed ($16$ out of $29$).")),
      include.rownames=F, caption.placement = "top")

# re-order the data (according to relative presence rate)
nsw_res_pa$pres_rate <- (nsw_res_pa$n_po / nsw_res_pa$n_pa)
nsw_res_pa$pres_rate[nsw_res_pa$pres_rate < 0.05] <- 0.04
nsw_res_po$pres_rate <- (nsw_res_po$n_po / nsw_res_po$n_pa)
nsw_res_po$pres_rate[nsw_res_po$pres_rate < 0.05] <- 0.04
nsw_res_pa <- nsw_res_pa[order(nsw_res_pa$pres_rate), ]
nsw_res_po <- nsw_res_po[order(nsw_res_po$pres_rate), ]

# create point sizes (according to n_po/n_pa)
pt.sizes_pa <- log(nsw_res_pa$n_pa) / max(log(nsw_res_pa$n_pa))
pt.sizes_po <- nsw_res_po$n_po / max(nsw_res_po$n_po)

plot.res <- 500
colfn <- colorRampPalette(c("blue1", "purple", "red1"))
png(filename = paste0(getwd(), "/app_results.png"), width = 6.2 * plot.res, height = 4 * plot.res, res = plot.res)
# par(mfrow = c(1, 2))
layout(matrix(1:3,ncol=3), width = c(2,2,0.5),height = c(1))
par(mar = c(5.1,4.1,4.1,0))
plot(nsw_res_pa$auc_pa, nsw_res_pa$auc_idm, pch = 21,
     bg = colfn(length(nsw_res_pa$pres_rate))[order(nsw_res_pa$pres_rate)],
     col = colfn(length(nsw_res_pa$pres_rate))[order(nsw_res_pa$pres_rate)],
     cex = pt.sizes_pa *3,
     # xlim = c(min(nsw_res_pa$auc_pa), 1),
     xlim = c(min(c(nsw_res_pa$auc_idm, nsw_res_pa$auc_pa)), 1),
     ylim = c(min(c(nsw_res_pa$auc_idm, nsw_res_po$auc_idm)), 1),
     # ylim = c(min(c(nsw_res_pa$auc_idm, nsw_res_po$auc_idm, nsw_res_po$auc_po, nsw_res_pa$auc_pa)), 1),
     xlab = "PA only", ylab = "IDM"#, main = "AUC (out-of-sample)", sub = parse(text = expression("atop(hat(p)[i]~on~the,Presence/Absence~survey~sites)"))
)
mtext(side = 3, line = 2, adj=0, cex=1.25, expression(AUC~(out-of-sample~~hat(p)[i])))
mtext(side = 3, line = 0.15, adj=0, cex=0.7, text = expression(IDM~italic("vs.")~PA~only))
# mtext(side = 3, line = 3, adj=0, cex=1, "AUC (out-of-sample)")
# mtext(side = 3, line = 1.25, adj=0, cex=0.7, text = expression(hat(p)[i]~on~the~presence/absence))
# mtext(side = 3, line = 0.25, adj=0, cex=0.7, text = "survey sites")
abline(a = 0, b = 1, lty = "dashed", lwd = 1.5)
par(mar = c(5.1,2.05,4.1,2.05))
plot(nsw_res_pa$auc_po, nsw_res_pa$auc_idm, pch = 22,
     bg = colfn(length(nsw_res_po$pres_rate))[order(nsw_res_po$pres_rate)],
     col = colfn(length(nsw_res_po$pres_rate))[order(nsw_res_po$pres_rate)],
     cex = (pt.sizes_po *3) + 0.5,
     # cex = pt.sizes_pa *3,
     # xlim = c(min(nsw_res_po$auc_po), 1),
     xlim = c(min(c(nsw_res_po$auc_idm, nsw_res_po$auc_po)), 1),
     # ylim = c(min(c(nsw_res_pa$auc_idm, nsw_res_po$auc_idm)), 1),
     ylim = c(min(c(nsw_res_pa$auc_idm, nsw_res_po$auc_idm)), 1),
     yaxt = "n", ylab = "",
     xlab = "PO only"#, main = parse(text = expression("atop(hat(p)[j]^(PO)*~on~the,Presence-Only~Data~and~Quadrature)"))#, ylab = "IDM"
)
mtext(side = 3, line = 0.15, adj=0, cex=0.7, text = expression(IDM~italic("vs.")~PO~only))
# mtext(side = 3, line = 1, adj=0, cex=0.7, text = expression(hat(p)[j]^(PO)*~on~the~presence-only))
# mtext(side = 3, line = 0.15, adj=0, cex=0.7, text = "data and quadrature")
abline(a = 0, b = 1, lty = "dashed", lwd = 1.75)
par(mar = c(0,0,0,0))
legend_image <- as.raster(matrix(rev(colfn(20)), ncol=1))
plot(c(0,5),c(0,12),type = 'n', axes = F,xlab = '', ylab = '', main = '') # Dummy
text(x=1.25, y = seq(0,4,l=4), labels = c("\u2264 0.05", "0.25", "1", "4"), pos = 4, adj = 0)
rasterImage(legend_image, 0,0,1,4)
for (i in seq(0,4,l=4)) {lines(x=c(1,1.25), y = c(i, i))}
lines(x=c(0,1,1,0,0), y = c(0,0,4,4,0))
text(x=1, y = 4.25, labels = expression(frac(n[PO], n[PA])), pos = 3, adj = c(0.5, 1))
text(x=1, y = 8.75, labels = expression(n[PO]), pos = 3, adj = c(0.5, 1))
text(x=1, y = 11.75, labels = expression(n[PA]), pos = 3, adj = c(0.5, 1))
# legend(x = 0, y = 12, legend = c("5", "50", "500"), pch = 16, pt.cex = (sqrt(c(5,50,500)) / max(sqrt(nsw_res_pa$n_pa))) * 5,
legend(x = 0, y = 12, legend = c("5", "50", "500"), pch = 16, pt.cex = (log(c(5,50,500)) / max(log(nsw_res_pa$n_pa))) *3,
       xpd = T, bty = "n", y.intersp=1.4, x.intersp = 1.2)#, title = expression(n[PA]))
legend(x = 0, y = 9, legend = c("5", "30", "60"), pch = 15, pt.cex = ((c(5,30,60)/max(nsw_res_po$n_po)) *3) + 0.5,
       xpd = T, bty = "n", y.intersp=1.4, x.intersp = 1.2)#, title = expression(n[PO]))
par(mar = c(5.1,4.1,4.1,2.1))
# par(mfrow = c(1, 1))
layout(matrix(1), width = 1, height = 1)
dev.off()

################################################################################
# for text #####################################################################
round(sum(timing.inla)/60)
round(sum(timing.scampr)/60)
sum(abs(nsw_res_pa$auc_diff) <= 0.02)
sum(abs(nsw_res_pa$auc_idm - nsw_res_pa$po) < 0.02)
sum(nsw_res_po$auc_diff > 0.02)
table(nsw_res_pa$auc_idm - nsw_res_pa$auc_po > 0.02)
table(nsw_res_pa$auc_idm - nsw_res_pa$auc_pa > 0.02)

png(filename = paste0(getwd(), "/append_app_results_pa.png"), width = 6.2 * plot.res, height = 3.7 * plot.res, res = plot.res)
par(mfrow = c(1,4))
par(mar = c(5.1,4,2.1,0))
with(nsw_res_pa, plot(pres_rate, auc_diff, ylab = expression(paste(Delta," AUC (IDM - PA only)")), xlab = "", log = "x", xaxt = "n", main = "A"))
axis(1, at = c(0.05, 0.25, 1, 4), labels = c("0.05", "0.25", "1", "4"))
title(xlab = expression(frac(n[PO], n[PA])), line = 3.7)
abline(h = 0, col = "red")
par(mar = c(5.1,4,2.1,0))
with(nsw_res_pa, plot(n_pa, auc_diff, yaxt = "n", xlab = expression(n[PA]), log = "x", ylab = "", xaxt = "n", main = "B"))
axis(1, at = c(5, 20, 100, 500))
abline(h = 0, col = "red")
par(mar = c(5.1,4,2.1,0))
with(nsw_res_pa, plot(n_po, auc_diff, yaxt = "n", xlab = expression(n[PO]), log = "x", ylab = "", xaxt = "n", main = "C"))
axis(1, at = c(2, 5, 20, 50))
abline(h = 0, col = "red")
par(mar = c(5.1,4,2.1,0))
with(nsw_res_pa, plot(m, auc_diff, yaxt = "n", xlab = expression(m), log = "x", ylab = "", xaxt = "n", main = "D"))
axis(1, at = c(1000, 1400, 1900))
abline(h = 0, col = "red")
dev.off()

png(filename = paste0(getwd(), "/append_app_results_po.png"), width = 6.2 * plot.res, height = 3.7 * plot.res, res = plot.res)
par(mfrow = c(1,4))
par(mar = c(5.1,4,2.1,0))
with(nsw_res_po, plot(pres_rate, auc_diff, ylab = expression(paste(Delta," AUC (IDM - PA only)")), xlab = "", log = "x", xaxt = "n", main = "A"))
axis(1, at = c(0.05, 0.25, 1, 4), labels = c("0.05", "0.25", "1", "4"))
title(xlab = expression(frac(n[PO], n[PA])), line = 3.7)
abline(h = 0, col = "red")
par(mar = c(5.1,4,2.1,0))
with(nsw_res_po, plot(n_pa, auc_diff, yaxt = "n", xlab = expression(n[PA]), log = "x", ylab = "", xaxt = "n", main = "B"))
axis(1, at = c(5, 20, 100, 500))
abline(h = 0, col = "red")
par(mar = c(5.1,4,2.1,0))
with(nsw_res_po, plot(n_po, auc_diff, yaxt = "n", xlab = expression(n[PO]), log = "x", ylab = "", xaxt = "n", main = "C"))
axis(1, at = c(2, 5, 20, 50))
abline(h = 0, col = "red")
par(mar = c(5.1,4,2.1,0))
with(nsw_res_po, plot(m, auc_diff, yaxt = "n", xlab = expression(m), log = "x", ylab = "", xaxt = "n", main = "D"))
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