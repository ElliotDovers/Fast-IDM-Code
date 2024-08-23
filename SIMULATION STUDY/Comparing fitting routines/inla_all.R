## Function to run IDM with simulated data TODO: ADAPT TO FIT ALL INLA MODELS IN ONE

inla_all <- function(structured_data, unstructured_data, quad, mesh, pred){
  
  # set the spde representation to be the mesh with defaults for parameters
  time.spde <- system.time(assign("spde", inla.spde2.matern(mesh)))
  
  # make A matrix for structured data
  time.projmesh_pa <- system.time(assign("structured_data_A", inla.spde.make.A(mesh = mesh, loc = as.matrix(structured_data[ , c("x","y")]))))
  
  # make A matrix for unstructured data
  time.projmesh_po <- system.time(assign("unstructured_data_A", inla.spde.make.A(mesh = mesh, loc = as.matrix(unstructured_data[ , c("x","y")]))))
  
  # make A matrix for quadrature points (NOTE THIS WILL BE DIAGONAL IF MESH == QUAD)
  time.projmesh_quad <- system.time(assign("quad_data_A", inla.spde.make.A(mesh = mesh, loc = as.matrix(quad[ , c("x","y")]))))
  
  # make A matrix for the prediction points
  time.projmesh_pred <- system.time(assign("pred_data_A", inla.spde.make.A(mesh = mesh, loc = as.matrix(pred[ , c("x","y")]))))
  
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
  
  # for model with PO only
  time.stk_po <- system.time(assign("stk_unstructured_data_po", inla.stack(data=list(y=y.pp, e = e.pp),
                                                                     effects=list(list(data.frame(interceptB=rep(1,nq+n)), env = c(quad$env, unstructured_data$env)), list(Bnodes=1:spde$n.spde)),
                                                                     A=list(1,A.pp),
                                                                     tag="po_data")))
  
  # for IDM
  time.stk_idm_po <- system.time(assign("stk_unstructured_data", inla.stack(data=list(y=cbind(y.pp, NA), e = e.pp),
                                      effects=list(list(data.frame(interceptB=rep(1,nq+n)), env = c(quad$env, unstructured_data$env)), list(uns_field=1:spde$n.spde, bias_field = 1:spde$n.spde)),
                                      A=list(1,A.pp),
                                      tag="po_data")))
  
  # stack for structured data
  # note intercept with different name
  
  # for model with PA only
  time.stk_pa <- system.time(assign("stk_structured_data_pa", inla.stack(data=list(y=structured_data$present, Ntrials = rep(1, nrow(structured_data))),
                                                                   effects=list(data.frame(interceptA=rep(1,nrow(structured_data)), env = structured_data$env), Bnodes=1:spde$n.spde),
                                                                   A=list(1,structured_data_A),
                                                                   tag="pa_data")))
  
  # for IDM
  time.stk_idm_pa <- system.time(assign("stk_structured_data", inla.stack(data=list(y=cbind(NA, structured_data$present), Ntrials = rep(1, nrow(structured_data))),
                                    effects=list(list(data.frame(interceptA=rep(1,length(structured_data$x)), env = structured_data$env)), list(str_field=1:spde$n.spde)),
                                    A=list(1,structured_data_A),
                                    tag="pa_data")))


  # create the prediction stacks
  
  # for model with PO only
  time.stk_pred1_po <- system.time(assign("stk_pred_response_po", inla.stack(data=list(y=NA),
                                                                       effects = list(list(data.frame(interceptB=rep(1,np))), env = pred$env, list(Bnodes=1:spde$n.spde)),
                                                                       A=list(1,1, pred_data_A),
                                                                       tag='pred_response')))
  
  # for model with PA only
  time.stk_pred1_pa <- system.time(assign("stk_pred_response_pa", inla.stack(data=list(y=NA, Ntrials = rep(1,np)),
                                                                       effects = list(list(data.frame(interceptA=rep(1,np))), env = pred$env, list(Bnodes=1:spde$n.spde)),
                                                                       A=list(1,1, pred_data_A),
                                                                       tag='pred_response')))
  
  # for IDM
  time.stk_pred1 <- system.time(assign("stk_pred_response", inla.stack(data=list(y=cbind(rep(NA, np), rep(NA, np))),
                                    effects = list(list(data.frame(interceptA=rep(1,np))), env = pred$env, list(uns_field=1:spde$n.spde)),
                                    A=list(1,1, pred_data_A),
                                    tag='pred_response')))
  
  # combine the stacks
  
  # for model with PO only
  time.stk_pred2_po <- system.time(assign("stk_po", inla.stack(stk_unstructured_data_po, stk_pred_response_po)))
  
  # for model with PA only
  time.stk_pred2_pa <- system.time(assign("stk_pa", inla.stack(stk_structured_data_pa, stk_pred_response_pa)))
  
  # for IDM
  time.stk0 <- system.time(assign("stk0", inla.stack(stk_unstructured_data, stk_structured_data)))
  time.stk_pred2 <- system.time(assign("stk", inla.stack(stk_unstructured_data, stk_structured_data, stk_pred_response)))

  # fit the model without predictions
  
  # for model with PO only
  result0_po <- inla(y ~ interceptB + env + f(Bnodes, model = spde) - 1,
                  family="poisson",
                  data=inla.stack.data(stk_unstructured_data_po),
                  control.predictor=list(A=inla.stack.A(stk_unstructured_data_po), compute=TRUE),
                  control.family = list(link = "log"),
                  E = inla.stack.data(stk_unstructured_data_po)$e,
                  control.compute = list(cpo=TRUE, waic = TRUE, dic = TRUE)
  )
  
  # for model with PA only
  result0_pa <- inla(y ~  interceptA + env + f(Bnodes, model = spde) -1,
                  family="binomial",
                  data=inla.stack.data(stk_structured_data_pa),
                  control.predictor=list(A=inla.stack.A(stk_structured_data_pa), compute=TRUE),
                  control.family = list(link = "cloglog"),
                  Ntrials = inla.stack.data(stk_structured_data_pa)$Ntrials,
                  control.compute = list(cpo=TRUE, dic = TRUE, waic = TRUE)
  )
  
  # for IDM
  result0 <- inla(y ~  interceptA + interceptB + env + f(uns_field, model = spde) + f(str_field, copy = "uns_field", fixed = TRUE) + f(bias_field, model = spde) -1,
                 family=c("poisson", "binomial"),
                 data=inla.stack.data(stk0),
                 control.predictor=list(A=inla.stack.A(stk0), compute=TRUE),
                 control.family = list(list(link = "log"),
                                       list(link = "cloglog")),
                 E = inla.stack.data(stk0)$e,
                 Ntrials = inla.stack.data(stk0)$Ntrials,
                 control.compute = list(cpo=TRUE, waic = TRUE, dic = TRUE)
  )
    
  # fit the model with predictions
  
  # for model with PO only
  result_po <- inla(y ~ interceptB + env + f(Bnodes, model = spde) - 1,
                        family="poisson",
                        data=inla.stack.data(stk_po),
                        control.predictor=list(A=inla.stack.A(stk_po), compute=TRUE),
                        control.family = list(link = "log"),
                        E = inla.stack.data(stk_po)$e,
                        control.compute = list(cpo=TRUE, waic = TRUE, dic = TRUE)
  )
  
  # for model with PA only
  result_pa <- inla(y ~  interceptA + env + f(Bnodes, model = spde) -1,
                 family="binomial",
                 data=inla.stack.data(stk_pa),
                 control.predictor=list(A=inla.stack.A(stk_pa), compute=TRUE),
                 control.family = list(link = "cloglog"),
                 Ntrials = inla.stack.data(stk_pa)$Ntrials,
                 control.compute = list(cpo=TRUE, dic = TRUE, waic = TRUE)
  )
  
  # for IDM
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
  index.pred.response_po <- inla.stack.index(stk_po, tag="pred_response")$data
  index.pred.response_pa <- inla.stack.index(stk_pa, tag="pred_response")$data
  index.pred.response <- inla.stack.index(stk, tag="pred_response")$data
  # get the predictions
  m.prd_po <- exp(result_po$summary.fitted.values$mean[index.pred.response_po])
  m.prd_pa <- exp(result_pa$summary.fitted.values$mean[index.pred.response_pa])
  m.prd <- exp(result$summary.fitted.values$mean[index.pred.response])
  # calculate the KL divergence metric
  KLdiv_po <- as.numeric(pred$quad.size %*% (pred$mu * log(pred$mu / m.prd_po))) - as.numeric(pred$quad.size %*% (pred$mu - m.prd_po))
  KLdiv_pa <- as.numeric(pred$quad.size %*% (pred$mu * log(pred$mu / m.prd_pa))) - as.numeric(pred$quad.size %*% (pred$mu - m.prd_pa))
  KLdiv <- as.numeric(pred$quad.size %*% (pred$mu * log(pred$mu / m.prd))) - as.numeric(pred$quad.size %*% (pred$mu - m.prd))
  # calculate the relative MAE metric from Simmonds et al 2020
  MAE_po <- mean(abs((m.prd_po-mean(m.prd_po))-(pred$mu - mean(pred$mu))))
  MAE_pa <- mean(abs((m.prd_pa-mean(m.prd_pa))-(pred$mu - mean(pred$mu))))
  MAE <- mean(abs((m.prd-mean(m.prd))-(pred$mu - mean(pred$mu))))
  
  # collate the results
  ret.obj <- data.frame(MODEL = c("PA", "PO", "IDM"), FIT = rep("INLA", 3), KL = c(KLdiv_pa, KLdiv_po, KLdiv), MAE = c(MAE_pa, MAE_po, MAE), #TIME = result$cpu.used[4], ALL_TIME = result$cpu.used[4]))
             TIME = c(mesh$timing.init + time.spde[3] + time.projmesh_pa[3] + time.stk_pa[3] + result0_pa$cpu.used[4],
                      mesh$timing.init + time.spde[3] + time.projmesh_po[3] + time.projmesh_quad[3] + time.stk_po[3] + result0_po$cpu.used[4],
                      mesh$timing.init + time.spde[3] + time.projmesh_pa[3] + time.projmesh_po[3] + time.projmesh_quad[3] + time.stk_idm_po[3] + time.stk_idm_pa[3] + time.stk0[3] + result0$cpu.used[4]),
             TIME_FINAL_FIT = c(result0_pa$cpu.used[4], result0_po$cpu.used[4], result0$cpu.used[4]),
             TIME_PRED = c(time.projmesh_pred[3] + time.stk_pred1_pa[3] + time.stk_pred2_pa[3] + result_pa$cpu.used[4],
                           time.projmesh_pred[3] + time.stk_pred1_po[3] + time.stk_pred2_po[3] + result_pa$cpu.used[4],
                           time.projmesh_pred[3] + time.stk_pred1[3] + time.stk_pred2[3] + result$cpu.used[4]),
             RHO_PA = c(sqrt(8)/exp(result_pa$summary.hyperpar["Theta2 for Bnodes", "mean"]),
                        NA,
                        sqrt(8)/exp(result$summary.hyperpar["Theta2 for uns_field", "mean"])),
             RHO_PO = c(NA,
                        sqrt(8)/exp(result_po$summary.hyperpar["Theta2 for Bnodes", "mean"]),
                        sqrt(8)/exp(result$summary.hyperpar["Theta2 for bias_field", "mean"])),
             BETA_ENV = c(result0_pa$summary.fixed["env", "mean"],
                          result0_po$summary.fixed["env", "mean"],
                          result0$summary.fixed["env", "mean"]),
             ACUTAL_DIM = c(rep(mesh$n, 2), mesh$n * 2)
  )
  row.names(ret.obj) <- 1:nrow(ret.obj)
  
  return(ret.obj)
}