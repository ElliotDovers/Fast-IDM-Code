home.wd <- getwd()
# Get the job array
tab <- read.csv("job_array.csv")

# determine job number from pbs script
job = 701

################################################################################
# Set parameters that define the scenarios:

# random seed / simulation number
seed = tab$sim_number[tab$job == job]
# environment range of effect
env.range = tab$env_corrleation_range[tab$job == job]
# latent range of effect
lat.range = tab$lat_corrleation_range[tab$job == job]
# bias field range of effect
bias.range = tab$bias_corrleation_range[tab$job == job]

################################################################################

source("sim_PO_PA_data.R")
dat_pa <- sim_occurrence_data(Intercept_po = -4,
                                       Intercept_pa = -1.25,
                                       rseed = seed,
                                       env.covariate.type = "random_field",
                                       presence.only.observer.bias.covariate.type = "random_field",
                                       presence.only.observer.bias.covariate.range = bias.range,
                                       env.covariate.range = env.range,
                                       latent.range = lat.range,
                                       latent.field = T,
                                       plotting = T
                                       
)
dat_po <- attr(dat_pa, "presence-only")
quad <- attr(dat_pa, "truth.grid")

# basis search
library(scampr)
quad$present <- 0
dat.scampr <- rbind(dat_po, quad)
pa0 <- scampr(formula = present ~ env, data = dat_pa, include.sre = F, sre.approx = "laplace", model.type = "PA")
pa <- basis.search.pa(pa0, domain.data = quad, return.model = F, max.basis.functions = 600)

po0 <- scampr(formula = present ~ env, data = dat.scampr, include.sre = F, sre.approx = "laplace", model.type = "PO")
po <- basis.search.po(po0, domain.data = quad, return.model = F, which.approx = "laplace")

plot.res <- 500

png(filename = paste0(home.wd, "/Results/Figures/bf_search_pa.png"), res = plot.res, width = 6.2 * plot.res, height = 4.5 * plot.res)
layout(mat = matrix(c(1,1,1,1,2,3,4,5), nrow = 2, ncol = 4, byrow = T), heights = c(0.7, 0.3), widths = rep(1/4, 4))
par(mar = c(4.1, 4.1, 4.1, 4.1))
plot(pa$k,
     pa$ll,
     xlab = "# Basis Functions", ylab = "PA log-Likelihood",#expression(l[PA](Y)),
     type = "l", main = "", col = "black", xaxt = "n" , bty = "n",
     ylim = range(pa$ll))
mtext("A", side = 3, line = 2, adj = 0)
points(pa$k, pa$ll)
abline(v = pa$k[pa$opt], col = "black", lty = "dashed")
axis(1, at = pa$k, labels = c("No GRF", paste(sqrt(pa$k[-1]), sqrt(pa$k[-1]), sep = "x")), gap.axis = 0.25)
legend(x = 0, y = max(pa$ll) + 5, legend = c("log-Likelihood", "Optimised Config.", "Domain Boundary", "Basis Functions"), lty = c("solid", "dashed", "solid", "solid"),
       col = c("black", "black", "red", "darkblue"), bty = "n", xpd = T, horiz = T)
par(new = T)
plot(pa$k, #pa$cpu,
     cumsum(pa$cpu),
     yaxt = "n", xaxt = "n", col = "darkgrey", xlab = "", ylab = "", bty = "n", type = "l",
     ylim = range(cumsum(pa$cpu)))
axis(4, ylim = range(cumsum(pa$cpu)), col = "darkgrey", col.axis = "darkgrey")
text(250, 0.6, "Cum. Comp. Time (sec)", srt = 90, xpd = TRUE, pos = 1, col = "darkgrey")
mtext("B", side = 1, line = 4, adj = 0)
# plot the basis function configurations
# par(mar = rep(0, 4))
par(mar = c(0, 0, 2.1, 0))
plims <- c(-50, 150)
for (i in 1:4) {
  if (i == 1) {
    plot(1, type = "n", axes = F, xlim = plims, ylim = plims, asp = 1, ylab = "", main = "", xlab = "")
    title(main = c("No GRF", "4x4", "5x5", "6x6")[i], line = -4.5)
    # text(x = 50, y = 50, labels = c("No GRF", "4x4", "5x5", "6x6")[i], adj = 0.5, cex = 1.5)
    lines(c(0,100,100,0,0), c(0,0,100,100,0), col = "red", lwd = 2)
    arrows(150, 50, x1 = 200, y1 = 50, length = 0.1, angle = 30)
  } else {
    plot(1, type = "n", axes = F, xlim = plims, ylim = plims, asp = 1, ylab = "", main = "", xlab = "")
    title(main = c("No GRF", "4x4", "5x5", "6x6")[i], line = 0.2)
    symbols(attr(pa, "bfs")[[i]]$x, attr(pa, "bfs")[[i]]$y, circles = attr(pa, "bfs")[[i]]$scale, add = T, inches = F, fg = "darkblue", lty = "dashed")
    lines(c(0,100,100,0,0), c(0,0,100,100,0), col = "red", lwd = 2)
    if (i != 4) {
      arrows(150, 50, x1 = 200, y1 = 50, length = 0.1, angle = 30)
    }
  }
}
dev.off()