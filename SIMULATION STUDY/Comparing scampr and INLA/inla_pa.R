## Function to run structured (PA only) models with simulated data

inla_pa <- function(structured_data, mesh, pred){
  
  # set the spde representation to be the mesh with defaults for parameters
  spde <- inla.spde2.matern(mesh)
  
  # make A matrix for structured data
  structured_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(structured_data[ , c("x","y")]))
  
  # make A matrix for the prediction points
  pred_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(pred[ , c("x","y")]))
  
  # PA Only model
  # Binomial GLM with a cloglog
  # One spatial field to account for missing predictors
  
  # set the number of prediction points
  np <- nrow(pred)

  # create stack including presence data from structured, have Ntrials instead of expected/weights
  stk_structured_data <- inla.stack(data=list(y=structured_data$present, Ntrials = rep(1, nrow(structured_data))),
                               effects=list(data.frame(interceptA=rep(1,nrow(structured_data)), env = structured_data$env), Bnodes=1:spde$n.spde),
                               A=list(1,structured_data_A),
                               tag="pa_data")
  
  # create the prediction stack
  stk_pred_response <- inla.stack(data=list(y=NA, Ntrials = rep(1,np)),
                                    effects = list(list(data.frame(interceptA=rep(1,np))), env = pred$env, list(Bnodes=1:spde$n.spde)),
                                    A=list(1,1, pred_data_A),
                                    tag='pred_response')

  # combine the stacks
  stk <- inla.stack(stk_structured_data, stk_pred_response)
  
  # fit the model
  result <- inla(y ~  interceptA + env + f(Bnodes, model = spde) -1,
                 family="binomial",
                 data=inla.stack.data(stk),
                 control.predictor=list(A=inla.stack.A(stk), compute=TRUE),
                 control.family = list(link = "cloglog"),
                 Ntrials = inla.stack.data(stk)$Ntrials,
                 control.compute = list(cpo=TRUE, dic = TRUE, waic = TRUE)
  )
  
  # create index to extract predictions
  index.pred.response <- inla.stack.index(stk, tag="pred_response")$data
  # get the predictions
  m.prd <- exp(result$summary.fitted.values$mean[index.pred.response])
  # calculate the KL divergence metric
  KLdiv <- as.numeric(pred$quad.size %*% (pred$mu * log(pred$mu / m.prd))) - as.numeric(pred$quad.size %*% (pred$mu - m.prd))
  # calculate the relative MAE metric from Simmonds et al 2020
  MAE <- mean(abs((m.prd-mean(m.prd))-(pred$mu - mean(pred$mu))))
  
  return(data.frame(MODEL = "PA", FIT = "INLA", KL = KLdiv, MAE = MAE, TIME = result$cpu.used[4], ALL_TIME = result$cpu.used[4]))
}
