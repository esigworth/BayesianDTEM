source("component_functions.R")
source("full_run.R")

start_time <- Sys.time()

load("doses.Rda")
load("seeds.Rda")
load("starts.Rda")
load("ends.Rda")

tau2 <- 1
fit_type <- "lbc"
i <- 1
burnin <- 10
thin <- 1
iters <- 150

fit <- paste(fit_type, tau2, sep = "_")

start <- starts[fit]
end <- ends[fit]

seeds <- seeds[start:end]

run_name <- paste(fit, "_", i, sep = "")

file_name <- paste(run_name, ".Rda", sep = "")

cat(paste("Started ", run_name, " at", start_time, "\n"))

res <- fit_one_type_tau(tau2 = as.numeric(tau2), i = i, doses = doses, 
                        nA = 75, nB = 75,
                        nsubjA = 100, nsubjB = 100,
                        minsubj = 50, maxsubj = 200,
                        fit_type = fit_type, seed = seeds[i],
                        burnin = burnin, thin = thin, iters = iters)

save(res, file = file_name)

end_time <- Sys.time()  

