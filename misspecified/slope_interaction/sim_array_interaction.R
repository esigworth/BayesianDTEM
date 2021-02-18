source("generate_sim_data_interaction.R")
source("fit_sim_models_accre.R")
source("equiv_functions_interaction.R")
source("run_sim_functions_interaction.R")

start_time <- Sys.time()

load("doses.Rda")
load("seeds_inter.Rda")
load("starts_inter.Rda")
load("ends_inter.Rda")

tau2 <- Sys.getenv('tau2')
fit_type <- Sys.getenv('fit_type')
i <- Sys.getenv('SLURM_ARRAY_TASK_ID')
burnin <- as.numeric(Sys.getenv('burnin'))
thin <- as.numeric(Sys.getenv('thin'))
iters <- as.numeric(Sys.getenv('iters'))

fit <- paste(fit_type, tau2, "inter", sep = "_")


start <- starts_inter[fit]
end <- ends_inter[fit]

seeds <- seeds_inter[start:end]

ind <- as.numeric(i)

run_name <- paste(fit, "_", i, sep = "")

file_name <- paste(run_name, ".Rda", sep = "")

cat(paste("Started ", run_name, " at", start_time, "\n"))

res <- fit_one_type_inter(tau2 = as.numeric(tau2), i = as.numeric(i), doses = doses, 
                        nA = 75, nB = 75,
                        nsubjA = 100, nsubjB = 100,
                        minsubj = 50, maxsubj = 200,
                        fit_type = fit_type, seed = seeds[ind],
                        burnin = burnin, thin = thin, iters = iters)

save(res, file = file_name)

end_time <- Sys.time()  

cat(paste("Finished ", run_name, " at", end_time, "took", end_time-start_time, "\n"))

