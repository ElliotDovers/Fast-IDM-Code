# FUNCTION FOR VALIDATION

# result = model output
# resolution = c(x,y)
# join.stack = joint stack of data and predictions
# model_type = indication which model result came from
# dat1 = original spatial field (truth)
# unstructured_data 
# structured_data
# choose table and/or plot
# qsize - quadrat size
# absolute - absolute or relative to mean
# dim - dimensions of domain

validation_function <- function(result, resolution, join.stack=NULL, mesh = NULL, model_type = c("unstructured", "unstructuredcov", "unstructuredsf", 
                                                                               "structured", "joint", "jointcov", "jointtwo"), 
                                unstructured_data=unstructured_data, structured_data=structured_data, dat1,
                                plotting = FALSE, summary_results = FALSE, qsize = qsize, absolute = TRUE, dim = dim,
                                is.scampr=FALSE, full_quad = quad){
  
#All comparisons are on the same scale as the truth is logged too!

model_type = match.arg(model_type)

# get the true field for the current resolution  
source('make_truth_grid.R')
tmp.grid <- make_truth_grid(resolution, dat1, c(dim[1],dim[2]), type='grid', absolute=TRUE)
pred.data <- full_quad[match(paste(tmp.grid$x, tmp.grid$y, sep = ":"), paste(full_quad$x, full_quad$y, sep = ":")), ]
structured.truth <- cbind(full_quad[match(paste(structured_data$x, structured_data$y, sep = ":"), paste(full_quad$x, full_quad$y, sep = ":")), ], pa = structured_data$presence)
unstructured.truth <- cbind(full_quad[match(paste(unstructured_data$x, unstructured_data$y, sep = ":"), paste(full_quad$x, full_quad$y, sep = ":")), ])

if (!is.scampr) { # INLA MODELS:
  
  # create index to extract predictions
  index.pred.response <- inla.stack.index(join.stack, tag="pred.response")$data
  index.struct.response <- inla.stack.index(join.stack, tag="structured_data")$data
  index.unstruct.response <- inla.stack.index(join.stack, tag="unstructured_data")$data
  
  # find the mean of the result and the standard deviation of predictions
  m.prd <- result$summary.fitted.values$mean[index.pred.response]
  sd.prd <- result$summary.fitted.values$sd[index.pred.response] # not used (keep in for legacy of original code)
  m.prd_struct <- result$summary.fitted.values$mean[index.struct.response]
  m.prd_unstruct <- result$summary.fitted.values$mean[index.unstruct.response]
  
  ## Extract the true quadrature ##
  # find the nearest-neighbouring points FROM the full quadrature TO the mesh vertices
  nn_mesh <- data.frame(Reduce('cbind', nearest.pixel(mesh$x, mesh$y, scampr::vec2im(full_quad$log_lambda_po, full_quad$x, full_quad$y))))
  colnames(nn_mesh) <- c("x", "y")
  # create a mesh data frame for KL divergence metric calculations
  quad.truth <- full_quad[match(paste(nn_mesh$x, nn_mesh$y, sep = ":"), paste(full_quad$x, full_quad$y, sep = ":")), ]
  pres.quad.id <- c(rep(0, nrow(mesh)), rep(1, nrow(unstructured_data)))
  
  # create a fit flag (indicating potentially poor convergence) # ALWAYS 0 FOR INLA IF IT GET'S TO THIS STAGE
  flag.fit <- 0
  flag.se <- 0
  
} else { # SCAMPR MODELS:
  library(scampr)
  if (model_type == "structured") {
    m.prd <- predict(result, newdata = pred.data)
    m.prd_struct <- 1 - exp(-exp(result$fitted.values))
    m.prd_unstruct <- NULL
  } else if (model_type %in% c("unstructured", "unstructuredcov")) {
    m.prd <- predict(result, newdata = pred.data)
    m.prd_struct <- NULL
    m.prd_unstruct <- exp(result$fitted.values)
    # extract truth at the quadrature
    pres.quad.id <- result$pt.quad.id
    quad.truth <- result$data[pres.quad.id == 0, ]
    
  } else if (model_type %in% c("joint", "jointcov", "jointtwo")) {
    m.prd <- predict(result, newdata = pred.data)
    m.prd_struct <- 1 - exp(-exp(attr(result$fitted.values, "PA")))
    m.prd_unstruct <- exp(result$fitted.values)
    # extract truth at the quadrature
    pres.quad.id <- result$pt.quad.id
    quad.truth <- result$data[pres.quad.id == 0, ]
  } else {
    stop("model type unsupported")
  }

  # included for legacy code
  sd.prd <- rep(NA, length(m.prd))
  
  # create a fit flag (indicating potentially poor convergence)
  flag.fit <- as.numeric(result$convergence != 0)
  flag.se <- as.numeric(result$se.flag != 0)
}

## alternative metrics #########################################################

# divergence in fields:
if (length(m.prd) == 0) {
  KLdiv <- NA
} else {
  KLdiv <- as.numeric(pred.data$quad.size %*% (exp(pred.data$log_mu_A) * log(exp(pred.data$log_mu_A) / exp(m.prd)))) - as.numeric(pred.data$quad.size %*% (exp(pred.data$log_mu_A) - exp(m.prd)))
}

if (length(m.prd_unstruct) == 0) {
  KLdiv_quad <- NA
  # KLdiv_grid <- NA
  # KL.pres <- NA
  # KL.quad <- NA
} else {
  # KLdiv_po <- sum(unstructured.truth$log_lambda_po) - sum(log(m.prd_unstruct[pres.quad.id == 1])) - as.numeric(full_quad$quad.size %*% exp(full_quad$log_lambda_po)) + as.numeric(quad.truth$quad.size %*% m.prd_unstruct[pres.quad.id == 0])
  KLdiv_quad <- as.numeric(quad.truth$quad.size %*% (exp(quad.truth$log_lambda_po) * log(exp(quad.truth$log_lambda_po) / m.prd_unstruct[pres.quad.id == 0]))) - as.numeric(quad.truth$quad.size %*% (exp(quad.truth$log_lambda_po) - m.prd_unstruct[pres.quad.id == 0]))
  # KL.pres <- sum(exp(unstructured.truth$log_lambda_po) * log(exp(unstructured.truth$log_lambda_po) / m.prd_unstruct[pres.quad.id == 1]))
  # KL.quad <- as.numeric(quad.truth$quad.size %*% (exp(quad.truth$log_lambda_po) - m.prd_unstruct[pres.quad.id == 0]))
  # new.quad.size <- prod(dim) / nrow(tmp.grid)
  # KLdiv_grid <- as.numeric(rep(new.quad.size, nrow(tmp.grid)) %*% (exp(tmp.grid$abundance) * log(exp(tmp.grid$abundance) / exp(m.prd)))) - as.numeric(rep(new.quad.size, nrow(tmp.grid)) %*% (exp(tmp.grid$abundance) - exp(m.prd)))
  # KLdiv_po <- sum(log(exp(unstructured.truth$log_lambda_po) / m.prd_unstruct[pres.quad.id == 1])) - as.numeric(quad.truth$quad.size %*% (exp(quad.truth$log_lambda_po) - m.prd_unstruct[pres.quad.id == 0]))
}

# logLoss:
if (length(m.prd_struct) == 0) {
  KLdiv_presprob <- NA
} else {
  presprob.prd <- m.prd_struct # already transformed via cloglog
  presprob.truth <- 1 - exp(-exp(structured_data$log_mu_A))
  trunc.presprob <- 1e-6
  presprob.prd[presprob.prd <= trunc.presprob] <- trunc.presprob
  presprob.prd[1 - presprob.prd <= trunc.presprob] <- 1 - trunc.presprob
  presprob.truth[presprob.truth <= trunc.presprob] <- trunc.presprob
  presprob.truth[1 - presprob.truth <= trunc.presprob] <- 1 - trunc.presprob
  KLdiv_presprob <- sum(presprob.truth * log(presprob.truth / presprob.prd) + (1 - presprob.truth) * log((1 - presprob.truth) / (1 - presprob.prd)))
}

################################################################################
  
# calculate differences
source('make_truth_grid.R')
if(absolute == TRUE){truth_grid <- make_truth_grid(resolution, dat1, c(dim[1],dim[2]), type='truth', absolute=TRUE)} else 
{truth_grid <- make_truth_grid(resolution, dat1, c(dim[1],dim[2]), type='truth', absolute=FALSE)}

if(absolute == TRUE){
  differences <- m.prd-truth_grid # calculate differences
  method = "Absolute"
}
if(absolute == FALSE){
  differences <- (m.prd-mean(m.prd))-truth_grid
  m.prd <- m.prd - mean(m.prd)
  sd.prd <- sd.prd - mean(sd.prd)
  method = "Relative"
}
  
if(plotting == TRUE){
png(paste0(model_type, " ", method, " validation.png"), height = 1000, width = 1000, pointsize = 25)
par(mfrow=c(2,2))
par(mar = c(5.1, 4.1, 4.1, 3.5))
# Plot truth on grid scale
image.plot(seq(resolution[1]/2,dim[1],resolution[1]),seq(resolution[2]/2,dim[2],resolution[2]), 
             matrix(truth_grid, ncol=dim[2]/resolution[2], nrow=dim[1]/resolution[1]), col=tim.colors(), xlab='', ylab='',main="Averaged truth",asp=1)
#predicted mean
image.plot(seq(resolution[1]/2,dim[1],resolution[1]),seq(resolution[2]/2,dim[2],resolution[2]), 
           matrix(m.prd, ncol=dim[2]/resolution[2], nrow=dim[1]/resolution[1]), col=tim.colors(),xlab='', ylab='',main="Predicted mean intensity",asp=1)
image.plot(seq(resolution[1]/2,dim[1],resolution[1]),seq(resolution[2]/2,dim[2],resolution[2]), 
           matrix(sd.prd, ncol=dim[2]/resolution[2], nrow=dim[1]/resolution[1]), col=tim.colors(),xlab='', ylab='',main="Predicted sd intensity",asp=1)
# relative differences
image.plot(seq(resolution[1]/2,dim[1],resolution[1]),seq(resolution[2]/2,dim[2],resolution[2]), 
           matrix(differences, ncol=dim[2]/resolution[2], nrow=dim[1]/resolution[2]), col=tim.colors(),xlab='', ylab='',main=paste0(model_type, " ", method, "\ndifferences"),asp=1)
dev.off()
}

if(plotting == FALSE){
  output <- list(truth_grid, m.prd)
}

if(summary_results == TRUE){
  MAE_differences <- abs(differences)
  correlation <- cor(m.prd, truth_grid)
  grid <- make_truth_grid(c(10,10), dat1, c(dim[1],dim[2]), type='grid')
  if (!is.scampr) {
    coefficients <- result$summary.fixed[,c(1,3,5,6)]
    timing <- result$cpu.used[4]
  } else {
    coefficients <- data.frame(mean = result$fixed.effects[ , 1],
                               lower = result$fixed.effects[ , 1] - (qnorm(1 - (0.05 / 2)) * result$fixed.effects[ , 2]),
                               upper = result$fixed.effects[ , 1] + (qnorm(1 - (0.05 / 2)) * result$fixed.effects[ , 2]),
                               mode = rep(NA, nrow(result$fixed.effects))
    )
    timing <- result$cpu[1]
  }
  #ONLY want to transform predictions NOT coefficients
  tab <- data.frame(Model = model_type,
                    MAE = mean(MAE_differences))
  attr(tab, "alt") <- data.frame(Model = model_type, KL_div = KLdiv, KL_presprob = KLdiv_presprob, KL_quad_po = KLdiv_quad, timing = as.numeric(timing), fit_flag = flag.fit, se_flag = flag.se)
  summary_results = list(tab,
                     correlation = correlation,
               coefficients = coefficients
               )
  names(summary_results) <- c("Proto-table", "correlation", "coefficients")
  if(plotting == TRUE){return(summary_results)}else{return(c(summary_results, output))}
}

}
