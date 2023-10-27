## Function to run IDM with Elith data

idm_scampr <- function(form, dat_pa, dat_po, domain){
  
  # # set the model formula
  # form <- occ ~ cti + disturb + mi + rainann + raindq + rugged + soildepth + 
  #                  soilfert + solrad + tempann + tempmin + topo + I(cti^2) + 
  #                  I(mi^2) + I(rainann^2) + I(raindq^2) + I(rugged^2) + I(soildepth^2) + 
  #                  I(solrad^2) + I(tempann^2) + I(tempmin^2) + I(topo^2)
  
  # fit the base model to initialise the basis search (this includes no spatial random effects)
  base_idm <- scampr(form, data = dat_po, bias.formula = ~ 1, pa.data = dat_pa, include.sre = F, model.type = "IDM", sre.approx = "laplace", latent.po.biasing = F)
  
  # fit the optimised IDM (using the basis search function)
  res_idm <- do.call("basis.search", list(object = base_idm, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10, max.basis.functions = 200))
  # res_idm <- basis.search(base_idm, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10, max.basis.functions = 200)
  
  # # predict the abundance across the domain
  # res_idm$pred <- predict(res_idm, newdata = domain)

  return(res_idm)
}