compute_severity <- function(pars, severity_data) {
  dt <- 0.25 # TODO: tidy this up
  ## WARNING: this is a hack
  list2env(pars, environment())
  
  # Time-varying probability of death (mechanistically of dying outside hospital)
  mu_D_date <- ZamCovid:::numeric_date(c("2020-08-15",
                                         "2020-12-01", "2021-03-15",
                                         "2021-05-01"))
  mu_D_vect <- c(1, mu_D_1, mu_D_1, mu_D_2)
  p_G_D_value <- p_G_D * mu_D_vect
  
  # Simplifying assumptions
  p_star_value <- 0.1 # Probability of being admitted with positive PCR
  p_H_value <- 0 # Probability of hospitalisation
  p_H_D_value <- 0 # Probability death in hospital
  
  # e. Simplifying assumption - note sero_sensitivity is 0.831
  ## TODO: this should not be the case; we are saying ALL go from sero_pre into sero_pos
  severity_data[severity_data$Name == "p_sero_pos", 2:17] <- 1
  
  severity <- ZamCovid::ZamCovid_parameters_severity(
    dt,
    severity_data,
    p_G_D = list(value = p_G_D_value, date = mu_D_date),
    p_H = list(value = p_H_value),
    p_H_D = list(value = p_H_D_value),
    p_star = list(value = p_star_value))
  
  severity
}


compute_progression <- function(pars, progression_data, rel_gammas) {
  dt <- 0.25
  ## WARNING: this is a hack
  list2env(pars, environment())
  
  k_parameters <-
    progression_data[grep("^k_", progression_data$parameter), ]
  gammas <-
    progression_data[grep("^gamma_", progression_data$parameter),
                     "value"]
  names(gammas) <-
    progression_data[grep("^gamma_", progression_data$parameter),
                     "parameter"]
  gammas <- as.list(gammas)
  # browser()
  
  gamma_dates <- rel_gammas$date
  
  gamma_E_value <- c(gammas$gamma_E * rel_gammas$wildtype[1], 
                     gammas$gamma_E * rel_gammas$wildtype[1],
                     gammas$gamma_E * rel_gammas$beta[1],
                     gammas$gamma_E * rel_gammas$beta[1],
                     gammas$gamma_E * rel_gammas$delta[1],
                     gammas$gamma_E * rel_gammas$delta[1],
                     gammas$gamma_E * rel_gammas$omicron[1])
  
  gamma_P_value <- c(gammas$gamma_P * rel_gammas$wildtype[1], 
                     gammas$gamma_P * rel_gammas$wildtype[1],
                     gammas$gamma_P * rel_gammas$beta[1],
                     gammas$gamma_P * rel_gammas$beta[1],
                     gammas$gamma_P * rel_gammas$delta[1],
                     gammas$gamma_P * rel_gammas$delta[1],
                     gammas$gamma_P * rel_gammas$omicron[1])
  
  gamma_C_1_value <- c(gammas$gamma_C_1 * rel_gammas$wildtype[1], 
                       gammas$gamma_C_1 * rel_gammas$wildtype[1],
                       gammas$gamma_C_1 * rel_gammas$beta[1],
                       gammas$gamma_C_1 * rel_gammas$beta[1],
                       gammas$gamma_C_1 * rel_gammas$delta[1],
                       gammas$gamma_C_1 * rel_gammas$delta[1],
                       gammas$gamma_C_1 * rel_gammas$omicron[1])
  
  gamma_C_2_value <- c(gammas$gamma_C_2 * rel_gammas$wildtype[2], 
                       gammas$gamma_C_2 * rel_gammas$wildtype[2],
                       gammas$gamma_C_2 * rel_gammas$beta[2],
                       gammas$gamma_C_2 * rel_gammas$beta[2],
                       gammas$gamma_C_2 * rel_gammas$delta[2],
                       gammas$gamma_C_2 * rel_gammas$delta[2],
                       gammas$gamma_C_2 * rel_gammas$omicron[2])
  
  gamma_H_R_value <- gammas$gamma_H_R
  gamma_H_D_value <- gammas$gamma_H_D
  gamma_PCR_pre_value <- 1 / 2
  gamma_PCR_pos_value <- 0.083
  
  gamma_R_value <- gammas$gamma_R
  gamma_sero_pos_value <- gammas$gamma_sero_pos
  
  # Time to diagnosis if admitted without test
  gamma_U_value <- 1 / 3
  
  progression <- ZamCovid::ZamCovid_parameters_progression(
    dt,
    gamma_E = list(value = gamma_E_value, date = gamma_dates),
    gamma_P = list(value = gamma_P_value, date = gamma_dates),
    gamma_C_1 = list(value = gamma_C_1_value, date = gamma_dates),
    gamma_C_2 = list(value = gamma_C_2_value, date = gamma_dates),
    gamma_H_D = list(value = gamma_H_D_value),
    gamma_H_R = list(value = gamma_H_R_value),
    gamma_U = list(value = gamma_U_value),
    gamma_PCR_pre = list(value = gamma_PCR_pre_value),
    gamma_PCR_pos = list(value = gamma_PCR_pos_value),
    gamma_R = list(value = gamma_R_value),
    gamma_sero_pos = list(value = gamma_sero_pos_value)
  )
  progression[k_parameters$parameter] <- k_parameters$value
  
  progression$k_sero_pre <- 1
  progression$gamma_sero_pre <- 1 / 13
  progression$k_PCR_pre <- 1
  progression$k_PCR_pos <- 1
  progression$k_sero_pos <- 1
  
  progression
}


compute_observation <- function(pars) {
  
  ## WARNING: this is a hack
  list2env(pars, environment())
  
  observation <- ZamCovid::ZamCovid_parameters_observation()
  observation$kappa_death_all <- 1 / alpha_D
  
  observation
}


make_transform <- function(baseline) {
  
  expected <- c("date", "region", "population", "restart_date",
                "beta_date", "beta_names",
                "severity_data", "progression_data",
                "sens_and_spec", "seed_size", "seed_pattern",
                "rel_gammas", "base_death_date", "base_death_value",
                "cross_immunity_date", "cross_immunity_value",
                ## Lots of vaccination things
                "rel_severity", "vaccine_eligibility_min_age",
                "vaccine_progression_rate",
                "vaccine_schedule", "vaccine_schedule_effect",
                "vaccine_uptake", "vaccine_mean_days_between_doses",
                "vaccine_index_dose2", "vaccine_days_to_effect",
                "n_doses")
  stopifnot(setequal(expected, names(baseline)))
  expected <- c("start_date", baseline$beta_names, "p_G_D", "alpha_D",
                "mu_D_1", "mu_D_2")
  
  function(pars) {
    
    stopifnot(setequal(expected, names(pars)))
    beta_value <- unname(pars[baseline$beta_names])
    pars <- as.list(pars) # using list access below
    
    progression <- compute_progression(pars, baseline$progression_data, 
                                       baseline$rel_gammas)
    severity <- compute_severity(pars, baseline$severity_data)
    observation <- compute_observation(pars)
    
    vaccine_schedule <- baseline$vaccine_schedule_effect
    
    stage_parameters <- function(vaccine_doses) {
      
      rel_severity <- baseline$rel_severity
      
      if (vaccine_doses == 0) {
        rel_severity <- lapply(
          rel_severity, function(x) x[, 1, drop = FALSE])
        vaccine_progression_rate <- 0
        vaccine_schedule <- NULL
        vaccine_index_dose2 <- NULL
        n_doses <- 2L # this is just the default, but has no effect
      }
      
      if (vaccine_doses == 2) {
        ## TODO: this is a hack! Need to set doses to 0 and instead inflate
        ##      states in simulation task to avoid carrying lots of empty
        ##      compartments
        rel_severity <- rel_severity
        vaccine_progression_rate <- baseline$vaccine_progression_rate
        vaccine_index_dose2 <- baseline$vaccine_index_dose2
        n_doses <- 2L
      }

      ZamCovid::ZamCovid_parameters(
        start_date = pars$start_date,
        population = baseline$population,
        beta_date = baseline$beta_date,
        beta_value = beta_value,
        
        base_death_date = baseline$base_death_date,
        base_death_value = baseline$base_death_value,
        
        cross_immunity_date = baseline$cross_immunity_date,
        cross_immunity_value = baseline$cross_immunity_value,
        
        severity = severity,
        progression = progression,
        observation = observation,
        sens_and_spec = baseline$sens_and_spec,
        
        initial_seed_size = baseline$seed_size,
        initial_seed_pattern = baseline$seed_pattern,
        
        rel_susceptibility = rel_severity$rel_susceptibility,
        rel_p_sympt = rel_severity$rel_p_sympt,
        rel_p_hosp_if_sympt = rel_severity$rel_p_hosp_if_sympt,
        rel_p_death = rel_severity$rel_p_death,
        rel_infectivity = rel_severity$rel_infectivity,
        
        vaccine_progression_rate = vaccine_progression_rate,
        vaccine_schedule = vaccine_schedule,
        vaccine_index_dose2 = vaccine_index_dose2,
        n_doses = n_doses
      )
    }
    stage_parameters(2)
    # stage_parameters(0)
  }
}
