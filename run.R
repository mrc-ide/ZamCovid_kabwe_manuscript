## ---------------------------
## Develop and test locally
## ---------------------------
root_dir <- paste0(orderly::orderly_config()$root, "/src/")
## ---------------------------
short_run <- TRUE
date <- "2021-09-30"
assumptions <- "central"
## can be: central or
# c("sero_sens_high", "sero_sens_low",
#   "imm_waning_fast", "imm_waning_slow",
#   "serorev_slow", "serorev_fast",
#   "observed_low", "observed_high")
deterministic <- FALSE
env_keep <- c("root_dir", "short_run", "date", "assumptions",
              "deterministic", "env_keep")
# ---------------------------------


## 1. ZamCovid_kabwe_data ----

# Develop
orderly::orderly_develop_start("ZamCovid_kabwe_data", use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_kabwe_data"))
file.edit("script.R")
# tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
orderly::orderly_run("ZamCovid_kabwe_data", use_draft = "newer")
rm(list = setdiff(ls(), env_keep))

#----


## 2. ZamCovid_kabwe_parameters ----

# Develop
orderly::orderly_develop_start(
  "ZamCovid_kabwe_parameters",
  parameters = list(assumptions = assumptions, date = date,
                    deterministic = deterministic),
  use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_kabwe_parameters"))
file.edit("script.R")
# Tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
orderly::orderly_run("ZamCovid_kabwe_parameters", parameters = list(
  assumptions = assumptions, date = date, deterministic = deterministic),
  use_draft = "newer")
rm(list = setdiff(ls(), env_keep))

#----


## 3. ZamCovid_kabwe_fits ----

# Develop
orderly::orderly_develop_start("ZamCovid_kabwe_fits",
                               parameters = list(region = "kabwe",
                                                 date = date,
                                                 short_run = short_run,
                                                 assumptions = assumptions,
                                                 deterministic = deterministic),
                               use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_kabwe_fits"))
file.edit("script.R")
# tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
orderly::orderly_run("ZamCovid_kabwe_fits",
                     parameters = list(region = "kabwe",
                                       date = date,
                                       short_run = short_run,
                                       assumptions = assumptions,
                                       deterministic = deterministic),
                     use_draft = "newer")
rm(list = setdiff(ls(), env_keep))

#----


## 4. ZamCovid_kabwe_sens_analysis ----

# Develop
orderly::orderly_develop_start("ZamCovid_kabwe_sens_analysis",
                               parameters = list(date = date,
                                                 short_run = short_run,
                                                 deterministic = deterministic),
                               use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_kabwe_sens_analysis"))
file.edit("script.R")
# tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
orderly::orderly_run("ZamCovid_kabwe_sens_analysis",
                     parameters = list(date = date,
                                       short_run = short_run,
                                       deterministic = deterministic),
                     use_draft = "newer")
rm(list = setdiff(ls(), env_keep))

#----


## 5. ZamCovid_kabwe_simulation ----

# Develop
orderly::orderly_develop_start("ZamCovid_kabwe_simulation",
                               parameters = list(n_par = 20,
                                                 fit_short_run = short_run),
                               use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_kabwe_simulation"))
file.edit("script.R")
# tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
orderly::orderly_run("ZamCovid_kabwe_simulation",
                     parameters = list(n_par = 20,
                                       fit_short_run = short_run),
                     use_draft = "newer")
rm(list = setdiff(ls(), env_keep))

#----

# To remove failed tasks from draft folder
orderly::orderly_cleanup(failed_only = TRUE)
