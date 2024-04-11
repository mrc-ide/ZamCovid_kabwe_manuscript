
transform_states <- function (state1, info1, info2) {
  fn <- function(state1, info1, info2, d1, d2, i1, i2, ny, 
                 x1, x2, nm) {
    if (length(d1) == 1) {
      x2[1, ] <- x1
    }
    else if (length(d2) == 3) {
      x2[, 1, , ] <- x1
    }
    else if (length(d2) == 4) {
      x2[, 1, , , ] <- x1
    }
    else {
      stop("Unexpected dimension of output")
    }
    x2
  }
  inflate_state(state1, info1, info2, fn)
}


compute_beta_step <- function(fits, control) {
  
  n_scenarios <- length(control$scenarios)
  nm_scenarios <- names(control$scenarios)
  
  dt <- fits$dt
  steps_per_day <- fits$steps_per_day
  step_end <- control$end_date * steps_per_day
  pars <- fits$pars
  
  # Make an array with one row per scenario
  rt_future <- do.call("rbind", replicate(
    n_scenarios, control$rt_scenarios, simplify = FALSE))
  rt_future$scenario <- nm_scenarios
  
  # Match step start date corresponding to each scenario
  scenario_step <- 
    sapply(control$scenarios, function(x) x$start_date) *
    steps_per_day
  rt_future$step_start <-
    scenario_step[match(rt_future$scenario, names(scenario_step))]
  
  # Get sample of fitted beta values, match to Rt values and complete if needed
  beta <- t(vapply(pars, "[[", numeric(length(pars[[1]]$beta_step)),
                   "beta_step"))
  rt <- t(fits$rt$Rt_general)
  eff_rt <- t(fits$rt$eff_Rt_general)
  
  n_beta_add <- step_end - ncol(beta)
  n_rt_add <- max(scenario_step) / steps_per_day - ncol(rt)
  
  if (n_beta_add > 0) {
    beta <- cbind(beta, matrix(rep(beta[, ncol(beta)], n_beta_add), nrow(beta)))
  }
  
  if (n_rt_add > 0) {
    rt <- cbind(rt, matrix(rep(rt[, ncol(rt)], n_rt_add), nrow(rt)))
    eff_rt <- cbind(eff_rt, matrix(rep(eff_rt[, ncol(eff_rt)], n_rt_add), nrow(eff_rt)))
  }

  step_rt <- unique(scenario_step)
  date_rt <- step_rt / steps_per_day
  final_beta_fits <- beta[, step_rt]
  
  final_rt_fits <- list(rt = rt, eff_rt = eff_rt)
  
  ## beta/Rt ratio
  rho <- rt_to_rho(final_rt_fits, final_beta_fits, date_rt)#scenario_step / steps_per_day)
  rt_future <- rt_future %>%
    dplyr::group_by(scenario) %>%
    dplyr::mutate(step_end = c(step_start[-1L] - 1L, step_end)) %>%
    dplyr::ungroup()
  
  control$beta_step <-
    setNames(
      lapply(
        nm_scenarios,
        create_beta_step,
        dt, rt_future, step_end, final_beta_fits, rho, n_beta_add, pars,
        control, beta),
      nm_scenarios
    )
  control
}


## rho = beta / rt
rt_to_rho <- function(fits, final_beta_fits, date_rt) {
  n_particles <- nrow(fits$rt)
  # date_rt <- unique(date_rt)
  n_restart <- length(date_rt)
  f <- function(rt) {
    array(
      matrix(final_beta_fits, c(n_particles, n_restart)) /
        matrix(fits[[rt]][, date_rt], c(n_particles, n_restart)),
      c(n_particles, n_restart),
      list(NULL)
    )
  }
  setNames(lapply(c("rt", "eff_rt"), f), c("rt", "rt_eff"))
}


create_beta_step <- function(key, dt, rt_future, step_end, final_beta_fits,
                             rho, n_beta_add, pars, control_pars, beta) {
  
  message(sprintf("Computing beta_step for scenario '%s'", key))
  scen_event <- rt_future[rt_future$scenario == key, ]
  beta_index <- control_pars$scenarios[[key]]$beta_index
  future_beta <- final_beta_fits[, beta_index]
  beta[, seq(scen_event$step_start, scen_event$step_end)] <- future_beta
  beta
}


update_simulate_pars <- function(fits, control, vaccine_schedule) {
  
  pars <- fits$pars
  pars1 <- pars[[1]]
  baseline <- fits$base
  
  
  shift_doses <- function(schedule, days_to_effect) {
    n_groups <- dim(schedule$doses)[1]
    n_doses <- dim(schedule$doses)[2]
    n_days <- dim(schedule$doses)[3]
    
    shift1 <- function(i, j) {
      c(rep(0, days_to_effect[j]),
        schedule$doses[i,j,],
        rep(0, max(days_to_effect) - days_to_effect[j]))
    }
    
    schedule_effect <-
      vapply(seq_len(n_groups),
             function (i) {
               vapply(seq_len(n_doses), function(j) shift1(i, j),
                      numeric(n_days + max(days_to_effect)))},
             array(0, c(n_days + max(days_to_effect), n_doses)))
    schedule_effect <- aperm(schedule_effect, c(3, 2, 1))
    
    schedule$doses <- schedule_effect
    schedule
    
  }
  
  vaccine_schedule_effect <- lapply(vaccine_schedule, shift_doses,
                                    baseline$vaccine_days_to_effect)

  severity <- baseline$rel_severity
  
  update <- NULL
  for (i in names(vaccine_schedule_effect)) {
    update[[i]] <-
      ZamCovid::ZamCovid_parameters_vaccination(
        pars1$dt,
        pars1$N_tot,
        n_groups = 16L,
        rel_susceptibility = severity$rel_susceptibility,
        rel_p_sympt = severity$rel_p_sympt,
        rel_p_hosp_if_sympt = severity$rel_p_hosp_if_sympt,
        rel_p_death = severity$rel_p_death,
        rel_infectivity = severity$rel_infectivity,
        vaccine_schedule = vaccine_schedule_effect[[i]],
        vaccine_index_dose2 = baseline$vaccine_index_dose2,
        vaccine_progression_rate = baseline$vaccine_progression_rate,
        n_doses = baseline$n_doses)
    update[[i]]$rel_p_G_D <- update[[i]]$rel_p_death
    update[[i]]$rel_p_H_D <- update[[i]]$rel_p_death
  }
  
  pars_update <- NULL
  for (i in seq_along(update)) {
    nm <- names(update[i])
    
    for (j in seq_along(pars)) {
      pars[[j]][names(update[[i]])] <- update[[i]]
    }
    pars_update[[nm]] <- pars
  }
  
  ## Lastly update model parameters for healthcare capability simulation
  pars_update <- create_healthcare_counterfactual(pars_update, fits, control)
  pars_update
}


simulate_local <- function(pars, control, fits, n_threads, seed = 123) {
  
  # Ensure reproducible seed starting from RNG seed when distributing simulations
  # aross multiple nodes.
  seed <- dust::dust_rng_distributed_state(seed, n_nodes = length(pars))
  lapply(seq_along(pars), simulate_scen, pars, fits, control, n_threads, seed)
  
}


simulate_scen <- function(i, pars, fits, control, n_threads, seed) {
  
  el <- pars[[i]]
  st_index <- control$scenarios[[i]]$st_index
  rt_dates <- control$scenarios[[i]]$rt_critical_dates
  st <- fits$state[, , st_index, drop = FALSE]
  stp <- fits$scenario_step[i]
  
  message(sprintf("-----\nRunning scenario %d / %d", i, length(pars)))
  time <- system.time(
    ret <- simulate_scen_one(el, fits, st, rt_dates, stp, n_threads, seed[[i]]))
  message(sprintf("Finished scenario %d in %2.1f s", i, time[["elapsed"]]))
  
  ret
}


simulate_scen_one <- function(pars, fits, st, rt_dates, step_start, n_threads, seed) {
  
  date_start <- step_start / fits$steps_per_day
  date_end <- ZamCovid:::numeric_date(date_end)
  dates <- seq(date_start, date_end)
  steps <- dates * fits$steps_per_day
  step_end <- last(steps)
  pars <- matrix(pars)
  
  message("Creating dust object")
  obj <- ZamCovid::ZamCovid$new(pars, step_start, NULL, pars_multi = TRUE,
                                n_threads = n_threads, seed = seed)
  
  info <- obj$info()[[1]]
  index <- simulate_index(info)
  obj$update_state(state = st)
  obj$set_index(index$run)
  
  message("Simulating!")
  state <- obj$simulate(steps)
  
  
  ret <- list(date = dates,
              state = state[, , 1, ],
              index = index)
  
  message("Calculating Rt")
  idx_s <- grep("^S", names(index$run))[1:80]
  idx_r <- grep("^R", names(index$run))[1:80]
  
  rt <- simulate_rt(steps, ret$state[idx_s, , ], pars,
                    rt_dates, ret$state[idx_r, , ])
  
  message("Calculating vaccination")
  vaccination_summary <- simulate_calculate_vaccination(ret$state, index, pars[[1]])

  message("Calculating severity")
  severity <- simulate_calculate_severity(ret$state, index, pars[[1]])
  
  ret <- c(ret, rt, vaccination_summary, severity)
  
  ret
}



simulate_index <- function(info) {
  
  set_names <- function (x, nms) {
    names(x) <- nms
    x
  }
  index_S <- info$index$S
  names(index_S) <- paste0("S_", seq_along(index_S))
  
  index_D <- info$index$D
  names(index_D) <- paste0("D_", seq_along(index_D))
  
  index_I <- info$index$I_weighted
  names(index_I) <- paste0("I_", seq_along(index_I))
  
  index_infectious <- info$index$I_weighted
  names(index_infectious) <- paste0("I_", seq_along(index_infectious))
  
  index_n_vaccinated <- set_names(
    info$index$cum_n_vaccinated,
    paste0("n_vaccinated_", seq_along(index_S))
  )
  index_R <- info$index$R
  names(index_R) <- paste0("R_", seq_along(index_R))
  
  index <- c(ZamCovid::ZamCovid_index(info)$state,
             index_S, index_D, index_I, index_R,
             index_n_vaccinated, index_infectious)
  
  list(run = index,
       n_vaccinated = index_n_vaccinated,
       S = index_S,
       # I = index_I,
       # infectious = index_infectious,
       R = index_R,
       D = index_D)
}


simulate_rt <- function(steps, S, pars, critical_dates, R, weight_Rt = FALSE,
                        rt_type = c("Rt_general", "eff_Rt_general")) {

  dim_S <- dim(S)
  dim_R <- dim(R)
  x <- list(pars)
  
  rt <- ZamCovid::ZamCovid_Rt_trajectories(
    steps, S, pars, type = rt_type,
    initial_time_from_parameters = FALSE,
    # R = R,
    weight_Rt = weight_Rt)

  ret_rt <- rt[rt_type]
  for (i in rt_type) {
    if (weight_Rt) {
      ret_rt[[i]] <- mcstate::array_reshape(t(rt[[i]]), 1, dim_S[2:3])
      dimnames(ret_rt[[i]]) <- list(NULL, region_names, NULL)
    } else {
      ret_rt[[i]] <- ret_rt[[i]][, 1, ]
    }
  }
  ret_rt
}


simulate_calculate_vaccination <- function(state, index, pars) {
  
  n_groups <- length(pars$N_tot)
  n_strata <- dim(pars$rel_infectivity)[2]
  
  ## Get cumulative numbers in vaccine strata by age and time
  n_vaccinated <- apply(state[names(index$n_vaccinated), , , drop = FALSE],
                        c(1, 3), mean)
  n_vaccinated <-
    mcstate::array_reshape(n_vaccinated, 1L, c(n_groups, n_strata))
  
  ## Get number recovered in vaccine strata by age and time
  # R_raw: [age, vacc_strata, particle, time]
  R_raw <- mcstate::array_reshape(
    state[names(index$R), , , drop = FALSE], 1L, c(n_groups, n_strata))
  
  # take the sum over age
  R <- apply(R_raw, seq(2, 4), sum)
  # take the mean over the particles
  R <- apply(R, c(1, 3), mean)
  
  
  ## Output number of first and second doses
  idx_doses <- c("first_dose_no_effect" = 1, "first_dose_protection" = 2,
                 "second_dose_protection" = 3, "second_dose_waned" = 4)
  doses <- n_vaccinated[, idx_doses, ]
  dimnames(doses)[2] <- list(names(idx_doses))
  
  doses_inc <- aperm(apply(doses, c(1, 2), diff), c(2, 3, 1))
  doses_inc <- mcstate::array_bind(array(NA, c(dim(doses_inc)[-3], 1)),
                                   doses_inc)
  colnames(doses_inc) <- paste0(colnames(doses), "_inc")
  
  n_doses <- abind_quiet(doses, doses_inc, along = 2)
  
  list(n_vaccinated = n_vaccinated,
       n_doses = n_doses)
}

simulate_calculate_severity <- function(state, index, pars) {
  
  n_groups <- length(pars$N_tot)
  n_strata <- dim(pars$rel_infectivity)[2]
  
  ## Get number dead per vaccine strata, age and time
  D_raw <- mcstate::array_reshape(
    state[names(index$D), , , drop = FALSE], 1L, c(n_groups, n_strata))
  
  idx_doses <- c("first_dose_no_effect" = 1, "first_dose_protection" = 2,
                 "second_dose_protection" = 3, "second_dose_waned" = 4)
  D_raw <- D_raw[, idx_doses, , ]
  dimnames(D_raw)[2] <- list(names(idx_doses))
  
  deaths_inc <- state[grep("death", rownames(state), value = TRUE), , ]
  
  ifr <- state[grep("ifr", rownames(state), value = TRUE), , ]
  
  list(D_raw = D_raw,
       deaths_inc = deaths_inc,
       ifr = ifr)
}
