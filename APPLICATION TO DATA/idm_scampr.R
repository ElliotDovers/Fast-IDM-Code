## Function to run IDM with Elith data

idm_scampr <- function(form, dat_pa, pres, quad, domain){
  
  # setup the required data structure
  pres$quad.size <- 1e-6
  dat_po <- rbind(pres, quad)
  
  # fit the base model to initialise the basis search (this includes no spatial random effects)
  base_idm <- scampr(form, data = dat_po, bias.formula = ~ 1, pa.data = dat_pa, include.sre = F, model.type = "IDM", sre.approx = "laplace", latent.po.biasing = F)
  
  # fit the optimised IDM (using the basis search function)
  res_idm <- do.call("basis.search", list(object = base_idm, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10, max.basis.functions = 500))
  
  # predict the abundance across the domain
  res_idm$pred <- predict(res_idm, newdata = domain, include.bias.accounting = T)

  return(res_idm)
}