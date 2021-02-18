fit_mod <- function(model, dat, chains = 4, 
                    burnin = 5000, iterations = 20000, thinning = 2,
                    varnames = c("b", "prec_u", "u")) {
  
  mod <- jags.model(file=textConnection(model), data = dat, n.chains = chains)
  
  update(mod, burnin)
  
  samples <- coda.samples(mod, variable.names = varnames, 
                          n.iter = iterations, thin = thinning)
  
  return(samples)
}


