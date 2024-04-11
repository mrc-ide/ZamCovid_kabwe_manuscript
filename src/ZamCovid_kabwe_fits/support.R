read_parameters <- function(path) {
  ret <- list(info = read_csv(file.path(path, "info.csv")),
              prior = read_csv(file.path(path, "prior.csv")),
              proposal = read_csv(file.path(path, "proposal.csv")))
  ret
}


pars_info_single <- function(region, info) {
  
  cols <- c("region", "include", "name", "initial", "max", "integer")
  stopifnot(all(cols %in% colnames(info)))
  
  if (!(region %in% info$region)) {
    stop(sprintf("Did not find region '%s' in parameter info", region))
  }
  info <- info[info$region == region & info$include, ]
  info <- info[setdiff(names(info), c("region", "include"))]
  rownames(info) <- NULL
  info
}


pars_prior_single <- function(region, info, prior) {
  prior_cols <- c("region", "type", "name", "gamma_scale", "gamma_shape",
                  "beta_shape1", "beta_shape2")
  stopifnot(all(prior_cols %in% colnames(prior)))
  prior <- prior[prior$region == region & prior$name %in% info$name, ]
  prior <- prior[setdiff(names(prior), "region")]
  rownames(prior) <- NULL
  prior
}


pars_proposal_single <- function(region, info, proposal, kernel_scaling) {
  proposal <- proposal[proposal$region == region, ]
  stopifnot(all(unique(proposal$name) %in% unique(info$name)))
  proposal <- as.matrix(proposal[match(info$name, proposal$name), info$name])
  rownames(proposal) <- info$name
  proposal * kernel_scaling
}


load_transform <- function(path, region, assumptions) {
  e <- new.env()
  sys.source(file.path(path, "transform.R"), e)
  stopifnot(is.function(e$make_transform))
  make_transform <- e$make_transform
  base <- readRDS(file.path(path, "base.rds"))[[region]]
  transform <- make_transform(base)
  list(base = base, transform = transform)
}


pars_mcmc_single <- function(info, prior, proposal, transform) {
  pars_mcmc <- Map(
    mcstate::pmcmc_parameter,
    name = info$name,
    initial = info$initial,
    min = info$min,
    max = info$max,
    integer = info$integer,
    prior = lapply(split(prior, prior$name), make_prior))
  
  ret <- mcstate::pmcmc_parameters$new(pars_mcmc, proposal, transform)
  
  ## Try and transform a single case and see if it works:
  ret$model(ret$initial())
  
  ret
}


make_prior <- function(d) {
  if (d$type == "gamma") {
    ## TODO: as_duration was droppd from here as never used, but if it
    ## is, then we'd transform p to 1/p
    shape <- d$gamma_shape
    scale <- d$gamma_scale
    function(p) {
      dgamma(p, shape = shape, scale = scale, log = TRUE)
    }
  } else if (d$type == "beta") {
    shape1 <- d$beta_shape1
    shape2 <- d$beta_shape2
    function(p) {
      dbeta(p, shape1 = shape1, shape2 = shape2, log = TRUE)
    }
  } else if (d$type == "null") {
    NULL
  } else {
    stop("Unknown prior type")
  }
}


fit_pars_load <- function(path, region, assumptions, short_run, deterministic) {
  
  if (short_run) {
    kernel_scaling <- 0.5
  } else {
    kernel_scaling <- ifelse(deterministic, 0.1, 0.2)
  }
  
  parameters <- read_parameters(path)
  info <- pars_info_single(region, parameters$info)
  prior <- pars_prior_single(region, info, parameters$prior)
  proposal <- pars_proposal_single(region, info, parameters$proposal,
                                   kernel_scaling)
  
  dat <- load_transform(path, region, assumptions)
  base <- dat$base
  transform <- dat$transform
  
  mcmc <- pars_mcmc_single(info, prior, proposal, transform)
  
  list(region = region,
       assumptions = assumptions,
       info = info,
       prior = prior,
       proposal = proposal,
       transform = transform,
       raw = parameters,
       base = base,
       mcmc = mcmc)
}


set_control <- function(short_run, deterministic, n_particles = 192,
                        workers = TRUE, n_threads = NULL, mcmc_path = NULL,
                        restart_date = NULL) {
  
  ## MCMC control (only applies if short_run = FALSE)
  if (deterministic) {
    burnin <- 1000
    n_mcmc <- 10000
    n_sample <- 1000
    n_chains <- 8
  } else {
    burnin <- 1000
    n_mcmc <- 3000
    n_sample <- 500
    n_chains <- 4
  }
  
  if (short_run) {
    burnin <- min(20, burnin)
    n_particles <- min(96, n_particles)
    n_mcmc <- min(50, n_mcmc)
    n_sample <- min(50, n_mcmc)
    n_chains <- min(4, n_chains)
  }
  
  
  n_steps_retain <- ceiling(n_sample / n_chains)
  parallel <- control_parallel(n_chains, workers, n_threads)
  if (!is.null(restart_date)) {
    restart_date <- ZamCovid:::numeric_date(restart_date)
  }
  
  pmcmc <- mcstate::pmcmc_control(n_steps = n_mcmc,
                                  n_chains = n_chains,
                                  n_burnin = burnin,
                                  save_state = TRUE,
                                  save_trajectories = TRUE,
                                  save_restart = restart_date,
                                  restart_match = TRUE,
                                  progress = interactive(),
                                  rerun_every = 100,
                                  n_threads_total = parallel$n_threads_total,
                                  n_workers = parallel$n_workers,
                                  use_parallel_seed = TRUE,
                                  nested_step_ratio = 1,
                                  rerun_random = TRUE,
                                  n_steps_retain = n_steps_retain,
                                  path = mcmc_path)

  n_threads <- parallel$n_threads_total / parallel$n_workers
  
  particle_filter <- list(n_particles = n_particles,
                          n_threads = n_threads,
                          seed = NULL)
  
  list(pmcmc = pmcmc,
       particle_filter = particle_filter)
}


control_parallel <- function(n_chains, workers, n_threads,
                             verbose = TRUE) {
  n_threads <- n_threads %||% control_cores()
  max_workers <- 4
  
  if (workers) {
    pos <- seq_len(max_workers)
    n_workers <- max(pos[n_threads %% pos == 0 & pos <= n_chains])
  } else {
    n_workers <- 1L
  }
  
  if (verbose) {
    message(sprintf("Running on %d workers with %d threads",
                    n_workers, n_threads))
  }
  list(n_threads_total = n_threads, n_workers = n_workers)
}


control_cores <- function() {
  as.integer(Sys.getenv("CONTEXT_CORES",
                        Sys.getenv("MC_CORES",
                                   getOption("mc.cores", 1))))
}


parse_data <- function(dat, fit_sero = FALSE, fit_deaths = FALSE,
                       fit_pcr = FALSE, sero_by_age = FALSE,
                       deaths_by_age = FALSE, pcr_by_age = FALSE) {
  
  # A bit of pre-processing to make sure minimum of expected columns are there
  sero_age_bands <- as.vector(
    outer(c("sero_pos_", "sero_tot_"),
          c("15_19", "20_29", "30_39", "40_49", "50_plus"),
          paste0))
  pcr_age_bands <- as.vector(
    outer(c("pcr_pos_", "pcr_tot_"),
          c("15_19", "20_29", "30_39", "40_49", "50_plus"),
          paste0))
  deaths_age_bands <- paste0(
    "deaths_", c("0_14", "15_39", "40_59", "60_plus"))
  
  expected <- c("date", "date_string", "region",
                "deaths_all", "sero_tot_over15", "sero_tot_all",
                "sero_pos_over15", "sero_pos_all", "pcr_pos_all",
                "pcr_tot_all",
                # Age-disaggregated columns
                sero_age_bands, pcr_age_bands, deaths_age_bands,
                # Some dummy columns - we might experiment with this later
                "hosp_admissions", "deaths_hosp")
  
  dat$region <- "kabwe"
  dat$date_string <- dat$date
  dat$date <- ZamCovid:::numeric_date(dat$date)
  dat$sero_tot_over15 <- dat$sero_tot_all
  dat$sero_pos_over15 <- dat$sero_pos_all
  
  missing <- setdiff(expected, colnames(dat))
  dat[, missing] <- NA_integer_
  dat <- dat[, expected]
  dat <- dat[as.Date(dat$date_string) <= as.Date("2021-09-30"), ]
  
  # Now check which data we are fitting to and avoid double-fitting to
  # age-disaggregated and aggregated data
  if (fit_sero) {
    
    if (sero_by_age) {
      stopifnot(!all(is.na(dat[, sero_age_bands])))
      dat[, c("sero_pos_over15", "sero_tot_over15")] <- NA_integer_
      
      dat <- dat %>%
        mutate(sero_pos_all = case_when(
          !is.na(sero_pos_50_plus) ~ NA_integer_,
          TRUE ~ sero_pos_all)) %>%
        mutate(sero_tot_all = case_when(
          !is.na(sero_tot_50_plus) ~ NA_integer_,
          TRUE ~ sero_tot_all))
      
    } else {
      stopifnot(!all(is.na(dat$sero_pos_over15)) && !all(is.na(dat$sero_tot_over15)))
      dat[, sero_age_bands] <- NA_integer_
    }
    
  } else {
    dat[, c("sero_pos_over15", "sero_tot_over15", sero_age_bands)] <- NA_integer_
  }
  
  
  if (fit_deaths) {
    
    if (deaths_by_age) {
      stopifnot(!all(is.na(dat[, deaths_age_bands])))
      dat$deaths_all <- NA_integer_
    } else {
      stopifnot(!all(is.na(dat$deaths_all)))
      dat[, deaths_age_bands] <- NA_integer_
    }
    
    ## NA deaths data in early 2020 to avoid inferring an earlier than
    ## expected start to the epidemic
    dat[as.Date(dat$date_string) < as.Date("2020-03-15"),
        c("deaths_all", deaths_age_bands)] <- NA_integer_
    
  } else {
    dat[, c("deaths_all", deaths_age_bands)] <- NA_integer_
  }
  
  
  if (fit_pcr) {
    
    if (pcr_by_age) {
      stopifnot(!all(is.na(dat[, pcr_age_bands])))
      dat <- dat %>%
        mutate(pcr_pos_all = case_when(
          !is.na(pcr_pos_50_plus) ~ NA_integer_,
          TRUE ~ pcr_pos_all)) %>%
        mutate(pcr_tot_all = case_when(
          !is.na(pcr_tot_50_plus) ~ NA_integer_,
          TRUE ~ pcr_tot_all))
    } else {
      stopifnot(!all(is.na(dat$pcr_pos_all)) && !all(is.na(dat$pcr_tot_all)))
      dat[, pcr_age_bands] <- NA_integer_
    }
    
  } else {
    dat[, c("pcr_pos_all", "pcr_tot_all", pcr_age_bands)] <- NA_integer_
  }
  
  # Final checks so data conforms to compare function
  ZamCovid::ZamCovid_check_data(dat)
  
  dat
}


ZamCovid_particle_filter <- function(data, pars, control, initial = NULL, 
                                     initial_date = 0, deterministic = FALSE) {
  
  # Transformed model parameters
  p <- pars$model(pars$initial())
  if (inherits(p, "multistage_parameters")) {
    p <- p[[1]]$pars
  }
  steps_per_day <- p$steps_per_day
  
  # Initial model states
  if (is.null(initial)) {
    initial <- ZamCovid::ZamCovid_initial
  } else if (is.matrix(initial)) {
    initial <- mcstate::particle_filter_initial(initial)
  } else {
    stop("initial must be a matrix")
  }
  
  # Data in particle filter format
  data <- mcstate::particle_filter_data(data, "date",
                                        steps_per_day, initial_date)
  
  # Particle filter
  if (deterministic) {
    ## Not checked yet! Default is FALSE
    mcstate::particle_deterministic$new(
      data = data,
      model = ZamCovid::ZamCovid,
      compare = ZamCovid::ZamCovid_compare,
      index = ZamCovid::ZamCovid_index,
      initial = initial,
      n_threads = control$n_threads)
  } else {
    mcstate::particle_filter$new(
      data = data,
      model = ZamCovid::ZamCovid,
      n_particles = control$n_particles,
      compare = ZamCovid::ZamCovid_compare,
      index = ZamCovid::ZamCovid_index,
      initial = initial,
      n_threads = control$n_threads,
      seed = control$seed)
  }
}


fit_run <- function(pars, filter, control) {
  message("Running chains - this will take a while!")
  initial <- replicate(control$n_chains,
                       pars$mcmc$propose(pars$mcmc$initial(), 1))
  ret <- mcstate::pmcmc(pars$mcmc, filter, initial = initial, control = control)
  ## Add some additional version information
  data <- pars$mcmc$model(ret$pars[1, ])
  base <- pars$base
  region <- base$region
  pars_names <- pars$mcmc$names()
  
  info <- ret$predict$filter$model$new(data, 0, 1)$info()
  
  ret$info <- list(version = packageVersion("ZamCovid"),
                   info = info,
                   data = data,
                   date = base$date,
                   restart_date = base$restart_date,
                   region = region,
                   beta_date = base$beta_date,
                   epoch_dates = base$epoch_dates,
                   pars = pars_names)
  
  ret
}
