ZamCovid_fit_process <- function(samples, parameters, data_full, data_fit) {
  
  region <- samples$info$region
  
  samples$trajectories$date <-
    samples$trajectories$time / samples$trajectories$rate
  date <- ZamCovid:::numeric_date_as_date(tail(samples$trajectories$date, 1))
  
  message("Computing Rt")
  rt <- calculate_ZamCovid_Rt(samples)
  
  message("Computing severity")
  severity <- calculate_severity(samples)
  
  message("Computing parameter MLE and covariance matrix")
  parameters_new <- fit_parameters(samples, parameters)
  
  if (!is.null(samples$restart)) {
    
    message("Computing restart information")
    restart <- process_restart(samples, parameters)
    restart_date <- max(restart$state$time)
    i <- samples$trajectories$date <= restart_date
    parent_rt <- rt_filter_time(rt, i)
    
    restart$parent <- list(
      trajectories = trajectories_filter_time(samples$trajectories, i),
      rt = parent_rt,
      data = data,
      prior = parameters$raw$prior)
    
  } else {
    restart <- NULL
  }
  
  list(fit = list(
    samples = samples,
    restart = restart,
    rt = rt,
    severity = severity,
    parameters = parameters_new,
    data = list(fitted = data_fit, full = data_full))
  )
}


calculate_severity <- function(samples) {
  
  time <- samples$trajectories$time
  dates <- samples$trajectories$date
  state <- samples$trajectories$state
  
  ret <- calculate_severity_region(state, time, dates)
  ret
}


calculate_severity_region <- function(state, time, date) {
  
  # String vectors to formulate severity trajectory names needed
  sev_traj <- grep("^ifr", rownames(state), value = TRUE)
  
  severity <- list()
  
  for (s in sev_traj) {
    tmp <- t(state[s, , ])
    severity[[s]] <- tmp
  }
  
  severity$time <- time
  severity$date <- date
  class(severity) <- "IFR_t_trajectories"
  
  severity
}


calculate_ZamCovid_Rt <- function(samples, weight_Rt = FALSE) {
  time <- samples$trajectories$time
  info <- samples$info$info
  epoch_dates <- samples$info$epoch_dates
  
  state <- samples$trajectories$state
  pars <- samples$pars
  transform <- samples$predict$transform
  

  ret <- calculate_ZamCovid_Rt_region(pars, state, transform,
                                      time, info, epoch_dates, weight_Rt)
  ret
}


calculate_ZamCovid_Rt_region <- function(pars, state, transform,
                                         time, info, epoch_dates, weight_Rt,
                                         keep_strains_Rt) {
  index_S <- grep("^S_", rownames(state))
  index_R <- grep("^R_", rownames(state))
  
  S <- state[index_S, , , drop = FALSE]
  R <- state[index_R, , , drop = FALSE]
  
  dates <- time / 4
  
  n_pars <- nrow(pars)
  
  rt <- list(time = numeric(0),
             date = numeric(0),
             beta = numeric(0),
             eff_Rt_general = numeric(0),
             Rt_general = numeric(0))
  
  pars_model <- lapply(seq_rows(pars), function(j)
    transform(pars[j, ]))
  
  n_vacc_classes <- pars_model[[1]]$n_vacc_classes
  
  suffix <- paste0("_", c(ZamCovid:::model_age_bins()$start))
  S_nms <- get_names("S", list(n_vacc_classes), suffix)
  S1 <- S[S_nms, , , drop = FALSE]
  
  rt <- ZamCovid::ZamCovid_Rt_trajectories(
    time, S1, pars_model, R = NULL)

  class(rt) <- c("Rt_trajectories", "Rt")
  rt
}


fit_parameters <- function(samples, parameters) {

  keep <- function(x, region) {
    is.na(x) | x %in% region
  }
  
  region <- samples$info$region
  info <- parameters$raw$info[keep(parameters$raw$info$region, region), ]
  rownames(info) <- NULL
  
  i <- which.max(samples$probabilities[, "log_posterior"])
  initial <- samples$pars[i, ]
  info$initial[match(names(initial), info$name)] <- unname(initial)
  
  covariance <- cov(samples$pars)
  rownames(covariance) <- NULL
  proposal <- data_frame(region = samples$info$region,
                         name = colnames(covariance),
                         covariance)
  
  prior <- parameters$prior
  prior$region <- samples$info$region
  prior <- prior[, c("region", setdiff(names(prior), "region"))]
  rownames(prior) <- NULL
  
  
  parameters$info <- info
  parameters$prior <- prior
  parameters$proposal <- proposal
  
  parameters
}


get_names <- function(state_name, suffix_list, suffix0 = NULL) {
  if (is.null(suffix0)) {
    suffixes <- list()
  } else {
    suffixes <- list(suffix0)
  }
  for (i in seq_along(suffix_list)) {
    nm <- names(suffix_list)[[i]]
    if (length(nm) == 0) {
      nm <- ""
    }
    suffixes <- c(suffixes,
                  list(c("", sprintf("_%s%s", nm,
                                     seq_len(suffix_list[[i]] - 1L)))))
  }
  suffixes <- expand.grid(suffixes)
  nms <- apply(suffixes, 1,
               function(x) sprintf("%s%s",
                                   state_name, paste0(x, collapse = "")))
  nms
}


summarise_trajectory <- function(sample, nm) {
  dates <- sample$trajectories$date[-1]
  state <- t(sample$trajectories$state[nm, , -1])
  
  data.frame(
    date = ZamCovid:::numeric_date_as_date(dates),
    state = nm,
    mean = rowMeans(state),
    lb = matrixStats::rowQuantiles(state, probs = 0.025),
    ub = matrixStats::rowQuantiles(state, probs = 0.975)
  )
}


get_convergence_diagnostic <- function(dat) {
  
  conv_dx <- function(sample) {
    n_full_pars <- nrow(sample$pars_full)
    n_chains <- max(sample$chain)
    
    sample$chain_full <- rep(seq_len(n_chains), each = n_full_pars / n_chains)
    
    chains <- lapply(unname(split(data.frame(sample$pars_full),
                                  sample$chain_full)), coda::as.mcmc)
    
    rhat <- tryCatch(coda::gelman.diag(chains), error = function(e) NULL)
    if (!is.null(rhat)) {
      rhat <- round(max(rhat$psrf[, "Point est."]), 1)
    } else {
      rhat <- NA_real_
    }
    
    ess <- function(p) {
      traces <- matrix(p, ncol = n_chains)
      sum(coda::effectiveSize(coda::as.mcmc(traces)))
    }
    
    pars <- sample$pars_full
    nms <- colnames(pars)
    pars_ess <- lapply(nms, function (nm) {
      ess(pars[, nm])
    })
    pars_ess <- round(min(unlist(pars_ess)))
    
    data.frame(rhat, pars_ess)
  }
  
  conv_dx(dat$fit$samples)
}


process_restart <- function(samples, parameters) {
  pars <- fit_parameters(samples, parameters)
  pars$prior <- fit_process_restart_priors(samples$pars, pars)
  pars$sample <- samples$pars
  class(pars) <- "pars_pmcmc"
  
  pars$base <- parameters$base
  
  list(state = samples$restart,
       info = samples$info,
       pars = pars)
}


fit_process_restart_priors <- function(values, parameters, nms) {
  nms <- colnames(values)
  stopifnot(all(nms %in% parameters$info$name))
  
  info <- parameters$info
  prior <- parameters$prior
  
  if (length(dim(values)) == 3) {
    wrapper_multiregion <- function(nm, region) {
      if (is.na(region)) {
        x <- values[, nm, 1]
        prior <- prior[prior$name == nm & is.na(prior$region), ]
        info <- info[info$name == nm & is.na(info$region), ]
      } else {
        x <- values[, nm, region]
        prior <- prior[prior$name == nm & prior$region == region, ]
        info <- info[info$name == nm & info$region == region, ]
      }
      if (length(unique(x)) > 1) {
        prior <- fit_prior(x, info, prior)
      }
      prior
    }
    
    nms_fixed <- info$name[is.na(info$region)]
    nms_varied <- setdiff(info$name, nms_fixed)
    region <- last(dimnames(values))
    prior_fixed <- lapply(nms_fixed, wrapper_multiregion, NA)
    prior_varied <- unlist(lapply(region, function(r)
      lapply(nms_varied, wrapper_multiregion, r)),
      FALSE, FALSE)
    res <- c(prior_fixed, prior_varied)
  } else {
    wrapper_single <- function(nm) {
      x <- values[, nm]
      info <- info[info$name == nm, ]
      prior <- prior[prior$name == nm, ]
      if (length(unique(x)) > 1) {
        prior <- fit_prior(x, info, prior)
      }
      prior
    }
    res <- lapply(nms, wrapper_single)
  }
  
  do.call("rbind", res)
}


fit_prior <- function(x, info, prior) {
  ## Can't update null priors; these will eventually be replaced in
  ## the underlying priors task, but this requires work from Marc.
  if (prior$type == "null") {
    return(prior)
  }
  
  if (prior$type == "gamma") {
    ll <- function(theta) {
      -suppressWarnings(
        sum(dgamma(x, scale = theta[1], shape = theta[2], log = TRUE)))
    }
    cols <- c("gamma_scale", "gamma_shape")
  } else if (prior$type == "beta") {
    ll <- function(theta) {
      suppressWarnings(
        -sum(dbeta(x, shape1 = theta[1], shape2 = theta[2], log = TRUE)))
    }
    cols <- c("beta_shape1", "beta_shape2")
  } else {
    stop(sprintf("Unsupported prior type '%s'", prior$type))
  }
  
  fit <- optim(unlist(prior[cols]), ll)
  prior[cols] <- fit$par
  prior
}


trajectories_filter_time <- function(trajectories, i) {
  trajectories$time <- trajectories$time[i]
  trajectories$date <- trajectories$date[i]
  trajectories$predicted <- trajectories$predicted[i]
  if (length(dim(trajectories$state)) == 3) {
    trajectories$state <- trajectories$state[, , i, drop = FALSE]
  } else {
    trajectories$state <- trajectories$state[, , , i, drop = FALSE]
  }
  trajectories
}


rt_filter_time <- function(rt, i){
  rt_filter_time1 <- function(x) {
    if (length(dim(x)) == 2) {
      ret <- x[i, , drop = FALSE]
    } else {
      ret <- x[i, , , drop = FALSE]
    }
    ret
  }
  
  ret <- lapply(rt, rt_filter_time1)
  class(ret) <- class(rt)
  ret
}
