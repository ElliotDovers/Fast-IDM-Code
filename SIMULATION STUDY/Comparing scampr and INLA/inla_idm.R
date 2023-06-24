## Function to run IDM with simulated data

inla_idm <- function(structured_data, unstructured_data, quad, mesh, pred){
  
  # set the spde representation to be the mesh with defaults for parameters
  spde <- inla.spde2.matern(mesh)
  
  # make A matrix for structured data
  structured_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(structured_data[ , c("x","y")]))
  
  # make A matrix for unstructured data
  unstructured_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(unstructured_data[ , c("x","y")]))
  
  # make A matrix for quadrature points (NOTE THIS WILL BE DIAGONAL IF MESH == QUAD)
  quad_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(quad[ , c("x","y")]))
  
  # make A matrix for the prediction points
  pred_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(pred[ , c("x","y")]))
  
  # Joint model
  
  # One spatial field for shared effects and one for PO biasing
  # Uses Simpson approach for PP data
  # Binomial model for PA data
  # Using cloglog
  
  # set the number of integration points, presence points in PO data and prediction points
  nq <- nrow(quad)
  n <- nrow(unstructured_data)
  np <- nrow(pred)
  
  
  # change data to include 0s for nodes and 1s for presences
  y.pp <- rep(0:1, c(nq, n))
  
  # add expectation vector (area for integration points/nodes and 0 for presences)
  e.pp <- c(quad$quad.size, rep(0, n))
  
  # make the full A matrix for the PPM component
  A.pp <- rbind(quad_data_A, unstructured_data_A)

  # unstructured data stack with integration points
  
  stk_unstructured_data <- inla.stack(data=list(y=cbind(y.pp, NA), e = e.pp),
                                      effects=list(list(data.frame(interceptB=rep(1,nq+n)), env = c(quad$env, unstructured_data$env)), list(uns_field=1:spde$n.spde, bias_field = 1:spde$n.spde)),
                                      A=list(1,A.pp),
                                      tag="po_data")	
  
  # stack for structured data
  # note intercept with different name
  
  stk_structured_data <- inla.stack(data=list(y=cbind(NA, structured_data$present), Ntrials = rep(1, nrow(structured_data))),
                                    effects=list(list(data.frame(interceptA=rep(1,length(structured_data$x)), env = structured_data$env)), list(str_field=1:spde$n.spde)),
                                    A=list(1,structured_data_A),
                                    tag="pa_data")


  # create the prediction stack
  stk_pred_response <- inla.stack(data=list(y=cbind(rep(NA, np), rep(NA, np))),
                                    effects = list(list(data.frame(interceptA=rep(1,np))), env = pred$env, list(uns_field=1:spde$n.spde)),
                                    A=list(1,1, pred_data_A),
                                    tag='pred_response')
  
  # combine the stacks
  stk <- inla.stack(stk_unstructured_data, stk_structured_data, stk_pred_response)
  
  # fit the model
  result <- inla(y ~  interceptA + interceptB + env + f(uns_field, model = spde) + f(str_field, copy = "uns_field", fixed = TRUE) + f(bias_field, model = spde) -1,
                 family=c("poisson", "binomial"),
                 data=inla.stack.data(stk),
                 control.predictor=list(A=inla.stack.A(stk), compute=TRUE),
                 control.family = list(list(link = "log"),
                                       list(link = "cloglog")),
                 E = inla.stack.data(stk)$e,
                 Ntrials = inla.stack.data(stk)$Ntrials,
                 control.compute = list(cpo=TRUE, waic = TRUE, dic = TRUE)
  )
  
  # create index to extract predictions
  index.pred.response <- inla.stack.index(stk, tag="pred_response")$data
  # get the predictions
  m.prd <- exp(result$summary.fitted.values$mean[index.pred.response])
  # calculate the KL divergence metric
  KLdiv <- as.numeric(pred$quad.size %*% (pred$mu * log(pred$mu / m.prd))) - as.numeric(pred$quad.size %*% (pred$mu - m.prd))
  # calculate the relative MAE metric from Simmonds et al 2020
  MAE <- mean(abs((m.prd-mean(m.prd))-(pred$mu - mean(pred$mu))))
  
  return(data.frame(MODEL = "IDM", FIT = "INLA", KL = KLdiv, MAE = MAE, TIME = result$cpu.used[4], ALL_TIME = result$cpu.used[4]))
}