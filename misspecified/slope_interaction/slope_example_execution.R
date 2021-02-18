source("slope_component_functions.R")
source("slope_full_run.R")

start_time <- Sys.time()

load("doses.Rda")
load("seeds_inter.Rda")
load("starts_inter.Rda")
load("ends_inter.Rda")

tau2 <- 1
fit_type <- "lbc"
i <- 1
burnin <- 10
thin <- 2
iters <- 100

fit <- paste(fit_type, tau2, "inter", sep = "_")


start <- starts_inter[fit]
end <- ends_inter[fit]

seeds <- seeds_inter[start:end]


run_name <- paste(fit, "_", i, sep = "")

file_name <- paste(run_name, ".Rda", sep = "")

res <- fit_one_type_inter(tau2 = as.numeric(tau2), i = as.numeric(i), doses = doses, 
                          nA = 75, nB = 75,
                          nsubjA = 100, nsubjB = 100,
                          minsubj = 50, maxsubj = 200,
                          fit_type = fit_type, seed = seeds[i],
                          burnin = burnin, thin = thin, iters = iters)

save(res, file = file_name)

end_time <- Sys.time()  


