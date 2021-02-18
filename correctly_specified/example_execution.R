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

res1 <- fit_one_type_tau(tau2 = as.numeric(tau2), i = 2, doses = doses, 
                         nA = 25, nB = 25,
                         nsubjA = 30, nsubjB = 30,
                         minsubj = 20, maxsubj = 50,
                        fit_type = fit_type, seed = seeds[2],
                        burnin = burnin, thin = thin, iters = iters)

res2 <- fit_one_type_tau(tau2 = as.numeric(tau2), i = 3, doses = doses, 
                         nA = 25, nB = 25,
                         nsubjA = 30, nsubjB = 30,
                         minsubj = 20, maxsubj = 50,
                         fit_type = "luc", seed = seeds[3],
                         burnin = burnin, thin = thin, iters = iters)

res3 <- fit_one_type_tau(tau2 = as.numeric(tau2), i = 4, doses = doses, 
                         nA = 25, nB = 25,
                         nsubjA = 30, nsubjB = 30,
                         minsubj = 20, maxsubj = 50,
                         fit_type = "luc", seed = seeds[4],
                         burnin = burnin, thin = thin, iters = iters)
