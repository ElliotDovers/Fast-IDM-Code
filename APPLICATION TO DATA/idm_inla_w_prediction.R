## Function to run IDM with Elith data

idm_inla_w_prediction <- function(dat_pa, dat_po, domain){
  
  # set the mesh according to the domain points
  # mesh.time <- system.time(assign("mesh", inla.mesh.2d(loc.domain = domain[ , c("x","y")], max.edge=c(0.35,0.5), cutoff=0.05)))
  mesh <- inla.mesh.2d(loc.domain = domain[ , c("x","y")], max.edge=c(0.35,0.5), cutoff=0.05)
  
  # set the spde representation to be the mesh with defaults for parameters
  spde <- inla.spde2.matern(mesh)
  
  # make A matrix for structured data
  structured_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(dat_pa[ , c("x","y")]))
  
  # make A matrix for unstructured data
  unstructured_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(dat_po[dat_po$occ == 1 , c("x","y")]))
  
  # make A matrix for quadrature points
  quad_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(dat_po[dat_po$occ == 0 , c("x","y")]))
  
  # make A matrix for the prediction points
  pred_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(domain[ , c("x","y")]))
  
  # Joint model
  
  # One spatial field for shared effects and one for PO biasing
  # Uses Simpson approach for PP data
  # Binomial model for PA data
  # Using cloglog
  
  # set the number of integration points, presence points in PO data and prediction points
  nq <- sum(dat_po$occ == 0)
  n <- sum(dat_po$occ == 1)
  np <- nrow(domain)
  
  
  # change data to include 0s for nodes and 1s for presences
  y.pp <- dat_po$occ
  
  # add expectation vector (area for integration points/nodes and 0 for presences)
  e.pp <- dat_po$quad.size
  
  # make the full A matrix for the PPM component
  A.pp <- rbind(unstructured_data_A, quad_data_A)
  
  # need to make our own factor contrasts for INLA
  fixed.po <- data.frame(get.design.matrix(res_idm$formula, dat_po))
  fixed.pa <- data.frame(get.design.matrix(res_idm$formula, dat_pa))
  fixed.pred <- data.frame(get.design.matrix(res_idm$formula, domain))
  
  # unstructured data stack with integration points
  stk_unstructured_data <- inla.stack(data=list(y=cbind(y.pp, NA), e = e.pp),
                                      effects=list(list(data.frame(interceptB=rep(1,nq+n)),
                                                        cti = dat_po$cti,
                                                        disturb2 = fixed.po$disturb2,
                                                        disturb3 = fixed.po$disturb3,
                                                        disturb4 = fixed.po$disturb4,
                                                        mi = dat_po$mi,
                                                        rainann = dat_po$rainann,
                                                        raindq = dat_po$raindq,
                                                        rugged = dat_po$rugged,
                                                        soildepth = dat_po$soildepth,
                                                        soilfert2 = fixed.po$soilfert2,
                                                        soilfert3 = fixed.po$soilfert3,
                                                        solrad = dat_po$solrad,
                                                        tempann = dat_po$tempann,
                                                        tempmin = dat_po$tempmin,
                                                        topo = dat_po$topo
                                      ),
                                      list(uns_field=1:spde$n.spde, bias_field = 1:spde$n.spde)),
                                      A=list(1,A.pp),
                                      tag="po_data")
  
  # stack for structured data
  # note intercept with different name
  stk_structured_data <- inla.stack(data=list(y=cbind(NA, dat_pa$occ), Ntrials = rep(1, nrow(dat_pa))),
                                    effects=list(list(data.frame(interceptA=rep(1,nrow(dat_pa))),
                                                      cti = dat_pa$cti,
                                                      disturb2 = fixed.pa$disturb2,
                                                      disturb3 = fixed.pa$disturb3,
                                                      disturb4 = fixed.pa$disturb4,
                                                      mi = dat_pa$mi,
                                                      rainann = dat_pa$rainann,
                                                      raindq = dat_pa$raindq,
                                                      rugged = dat_pa$rugged,
                                                      soildepth = dat_pa$soildepth,
                                                      soilfert2 = fixed.pa$soilfert2,
                                                      soilfert3 = fixed.pa$soilfert3,
                                                      solrad = dat_pa$solrad,
                                                      tempann = dat_pa$tempann,
                                                      tempmin = dat_pa$tempmin,
                                                      topo = dat_pa$topo
                                    ),
                                    list(str_field=1:spde$n.spde)),
                                    A=list(1,structured_data_A),
                                    tag="pa_data")
  

  # create the prediction stack
  stk_pred_response <- inla.stack(data=list(y=cbind(rep(NA, np), rep(NA, np))),
                                    effects = list(list(data.frame(interceptA=rep(1,np)),
                                                        cti = domain$cti,
                                                        disturb2 = fixed.pred$disturb2,
                                                        disturb3 = fixed.pred$disturb3,
                                                        disturb4 = fixed.pred$disturb4,
                                                        mi = domain$mi,
                                                        rainann = domain$rainann,
                                                        raindq = domain$raindq,
                                                        rugged = domain$rugged,
                                                        soildepth = domain$soildepth,
                                                        soilfert2 = fixed.pred$soilfert2,
                                                        soilfert3 = fixed.pred$soilfert3,
                                                        solrad = domain$solrad,
                                                        tempann = domain$tempann,
                                                        tempmin = domain$tempmin,
                                                        topo = domain$topo
                                    ),
                                    list(uns_field=1:spde$n.spde)),
                                    A=list(1,1, pred_data_A),
                                    tag='pred_response')
  
  # combine the stacks
  stk <- inla.stack(stk_unstructured_data, stk_structured_data, stk_pred_response)
  
  # fit the model
  result <- inla(y ~ interceptA + interceptB + cti + disturb2 + disturb3 + disturb4 + mi + rainann + raindq + rugged + soildepth + 
                   soilfert2 + soilfert3 + solrad + tempann + tempmin + topo + I(cti^2) + 
                   I(mi^2) + I(rainann^2) + I(raindq^2) + I(rugged^2) + I(soildepth^2) + 
                   I(solrad^2) + I(tempann^2) + I(tempmin^2) + I(topo^2) + 
                   f(uns_field, model = spde) + f(str_field, copy = "uns_field", fixed = TRUE) + f(bias_field, model = spde) -1,
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
  # add in the predictions
  result$pred <- result$summary.fitted.values$mean[index.pred.response]
  
  return(result)
}