source("util.R")

version_check("ZamCovid", "0.1.3")

## Restart dates (Mondays) for counterfactuals are as follows:
# Vaccination 2020-04-20: 100th day after first full genome sequence
# Healthcare 2020-06-15: dexamethasone approved for covid (2020-06-18)
# Healthcare 2020-10-05: before second wave
# Vaccination 2020-12-14: vaccination roll-out in HIC settings
# Healthcare 2021-04-05: before third wave
# Vaccination 2021-04-19: vaccination roll-out in Zambia
restart_date <- c("2020-04-20", "2020-06-15", "2020-10-05",
                  "2020-12-14", "2021-04-05", "2021-04-19")

## Load all parameters from the last run; creates priors, and updates
## new entries into the proposal matrix as needed.
pars <- load_mcmc_parameters(assumptions, deterministic)

baseline <- list(create_baseline(restart_date, pars, assumptions))
names(baseline) <- "kabwe"

message("Writing parameters_info.csv")
write_csv(pars$info, "parameters_info.csv")
message("Writing parameters_proposal.csv")
write_csv(pars$proposal, "parameters_proposal.csv")
message("Writing parameters_prior.csv")
write_csv(pars$prior, "parameters_prior.csv")

message("Writing parameters_base.rds")
saveRDS(baseline, "parameters_base.rds")

message("Writing parameters_transform.R")
fs::file_copy("R/transform.R",
              "parameters_transform.R", overwrite = TRUE)
