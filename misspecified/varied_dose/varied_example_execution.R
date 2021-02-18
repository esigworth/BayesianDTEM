source("varied_component_functions.R")
source("varied_full_run.R")

start_time <- Sys.time()

load("doses.Rda")
load("seeds.Rda")
load("starts.Rda")
load("ends.Rda")

tau2 <- 1
fit_type <- "lbn"
i <- 1
burnin <- 10
thin <- 2
iters <- 100

fit <- paste(fit_type, tau2, "varied", sep = "_")


start <- starts[fit]
end <- ends[fit]

seeds <- seeds[start:end]

ind <- as.numeric(i)

run_name <- paste(fit, "_", i, sep = "")

file_name <- paste(run_name, ".Rda", sep = "")

res <- fit_one_varied(tau2 = as.numeric(tau2), i = as.numeric(i), doses = doses, 
                        nA = 75, nB = 75,
                        nsubjA = 100, nsubjB = 100,
                        minsubj = 50, maxsubj = 200,
                        fit_type = fit_type, seed = seeds[ind],
                        burnin = burnin, thin = thin, iters = iters)

save(res, file = file_name)

end_time <- Sys.time()  


