simulate_control <- function(dat, date_end, uptake = 0.7,
                             scale_ifr = c(0.15, 0.30, 0.45)) {
  
  # Restart objects dates for healthcare capability counterfactual scenarios
  healthcare_dates <- dat$restart$state$time[c(2, 3, 5)]
  names(healthcare_dates) <- c("hc_1st_wave", "hc_2nd_wave", "hc_3rd_wave")
  
  # Restart objects dates for vaccination counterfactual scenarios
  vaccine_start_dates <- dat$restart$state$time[c(1, 4, 6)]
  names(vaccine_start_dates) <- c("cepi", "high_income", "zambia_real")
  vaccine_step_start <- vaccine_start_dates * dat$steps_per_day
  end_date <- ZamCovid:::numeric_date(dat$parameters$base$date)
  
  restart_dates <- c(vaccine_start_dates, healthcare_dates)
  
  ## Create a list of weekly dates for simulation scenarios based on restart dates
  weeks <- NULL
  for (i in seq_along(restart_dates)) {
    ds <- ZamCovid:::numeric_date_as_date(restart_dates[i])
    de <- ZamCovid:::numeric_date_as_date(end_date)
    seq_d <- seq.Date(as.Date(ds), as.Date(de), "day")
    weeks[[i]] <- seq_d[lubridate::wday(seq_d, label = TRUE) == "Mon"]
  }
  names(weeks) <- names(restart_dates)

  
  
  ## Create a list of available dose schedules
  N_age <- dat$parameters$base$population$n
  N_tot <- sum(N_age)
  N_eligible <- round(sum(N_age[4:16]) + (N_age[3] / 5) * 3)
  
  
  ## We'll simulate vaccination scenarios based on what was achieved by income
  ## setting vs Zambia (factual). We will shift the schedules, so week one
  ## of data correponds to the roll-out date of counterfactual simulation.
  shift_coverage <- function(x) {
    n <- min(which(!is.na(x)))
    x <- c(x[-(seq(n))], rep(NA, n))
    zoo::na.locf(x)
  }
  
  
  ## Vaccination schedules using officially reported daily doses per-100
  ## population by average of income settings, as per Our World In Data
  schedules_weekly <- read_csv("inputs/owid-covid-data.csv") %>%
    select(date, name = location, value = total_vaccinations_per_hundred) %>%
    filter(name %in% c("Zambia", "High income", "Lower middle income",
                       "Low income", "Upper middle income")) %>%
    mutate(name = gsub("Zambia", "factual", name)) %>%
    mutate(date = lubridate::floor_date(as.Date(date), "week", week_start = 3),
           name = tolower(gsub(" ", "_", name))) %>%
    group_by(date, name) %>%
    summarise(value = mean(value, na.rm = TRUE) * 0.01) %>%
    ungroup() %>%
    pivot_wider() %>%
    filter(date >= as.Date("2021-01-10"), date < as.Date("2022-08-01")) %>%
    mutate(week_number = row_number(), .before = date) %>%
    select(!date) %>%
    mutate(across(!week_number, shift_coverage)) %>%
    filter(week_number <= length(weeks$cepi)) %>%
    pivot_longer(!week_number) %>%
    group_by(name) %>%
    mutate(value = c(diff(value), 0)) %>%
    pivot_wider()

  
  ## Calendar of school term dates. At the moment is not being used, as the
  ## analysis is retrospective and we assume fitted beta values absorbed
  ## changes in transmissibility due to schools being open/closed.
  ## For prospective scenarios beyond available data to fit, this can be
  ## used, making assumptions around how much school term dates affect beta.
  rt_school_dates <- read_csv("inputs/schools_schedule.csv") %>%
    mutate(date = lubridate::ymd(paste(year, month, day, sep = "-"))) %>%
    select(date)
  
  
  scenarios <- NULL
  for (i in seq_along(weeks)) {
    
    nm <- names(weeks)[i]
    start_date <- min(weeks[[i]])
    rt_critical_dates <- rt_school_dates %>%
      filter(date >= start_date, date <= date_end)
    rt_critical_dates <- unname(c(start_date, rt_critical_dates, date_end))
    rt_critical_dates <- sort(ZamCovid:::numeric_date(
      unique(rt_critical_dates[rt_critical_dates > start_date])))
    
    ret <- NULL
    if (nm %in% c("cepi", "high_income", "zambia_real")) {
      for (j in colnames(schedules_weekly)[-1L]) {
        
        # Vaccine scenario dates correspond to restart indexes 1, 4, 6
        st_index <- ifelse(i == 1, i, i * 2)
        nm1 <- paste(nm, j, sep = "_")
        
        times <- seq(1, length(weeks[[i]]))
        lng <- length(times)
        values <- round(c(data.frame(schedules_weekly[, j]))[[1]] * N_tot)
        
        scenarios[[nm1]] <- list(doses = values,
                                 start_date = ZamCovid:::numeric_date(start_date),
                                 st_index = st_index,
                                 beta_index = i,
                                 rt_critical_dates = rt_critical_dates,
                                 n_days = as.integer(max(weeks[[i]]) - min(weeks[[i]])))
      }
      
    } else {
      
      for (j in scale_ifr) {
        
        # Healthcare scenario dates correspond to restart indexes 2, 3, 5
        st_index <- ifelse(i == 4, i - 2, i * 2 - 7)
        beta_index <- ifelse(i == 4, 2, i)
        nm1 <- if (j == scale_ifr[1]) {
          paste(nm, "pessimistic", sep = "_")
        } else if (j == scale_ifr[2]) {
          paste(nm, "central", sep = "_")
        } else if (j == scale_ifr[3]) {
          paste(nm, "optimistic", sep = "_")
        }
        
        scenarios[[nm1]] <- list(start_date = ZamCovid:::numeric_date(start_date),
                                 st_index = st_index,
                                 beta_index = beta_index,
                                 rt_critical_dates = rt_critical_dates,
                                 n_days = as.integer(max(weeks[[i]]) - min(weeks[[i]])),
                                 scale_ifr = j)
      }
    }
  }
  
  ## We assume vaccines are allocated:
      # Anyone aged 12+
      # There is a 70% uptake amongst eligible
      # In decreasing priority by age group
      # With a 12 week gap between doses
  vaccine_eligibility_min_age <- 12
  uptake_by_age <- c(rep(0, 2), # no vaccination in <12
                     (3 / 5) * uptake, # no vaccination in 10-11yo
                     rep(uptake, 13))
  age_priority <- list(16, 15, 14, 13, 12, 11, 9:10, 7:8, 1:6)
  days_between_doses <- 7 * 12
  
  ## TODO: This is not being used at the moment; consider removing
  rt_seasonality_date_peak <- ZamCovid:::numeric_date("2020-07-01")
  rt_seasonality_amplitude <- 0
  rt_sd <- 0.001 # Log-normal standard deviation for generating Rt
  rt_schools_modifier <- 0 # Rt decrease when schools close (0.1 = 10% reduction)
  rt_scenarios <- read_csv("inputs/rt_scenarios.csv")
  rt_scenarios$date <- NULL
  
  
  list(scenarios = scenarios,
       end_date = end_date,
       
       vaccine_start = vaccine_start_dates,
       vaccine_eligibility_min_age = vaccine_eligibility_min_age,
       uptake_by_age = uptake_by_age,
       age_priority = age_priority,
       days_between_doses = days_between_doses,
       
       rt_seasonality_date_peak = rt_seasonality_date_peak,
       rt_seasonality_amplitude = rt_seasonality_amplitude,
       rt_sd = rt_sd,
       rt_schools_modifier = rt_schools_modifier,
       rt_schools_schedule = rt_school_dates,
       rt_scenarios = rt_scenarios)
  
}


fits_simulate_prepare <- function(fits, date_start, n_par) {

  restart <- fits$restart
  info <- restart$info
  
  ## Take a random sample of our parameters without replacement.
  n_par_fits <- nrow(restart$pars$sample)
  if (n_par > n_par_fits) {
    message(sprintf(
      "Reducing n_par from %d to %d as too few available in fits",
      n_par, n_par_fits))
    n_par <- n_par_fits
  }
  
  i <- sort(sample(n_par_fits, n_par, replace = FALSE))
  pars_mcmc <- restart$pars$sample[i, , drop = FALSE]
  state <- restart$state$state[, i, , drop = TRUE]
  
  message("Creating odin parameters from sampled parameters")
  pars <- apply(pars_mcmc, 1, restart$pars$transform)
  
  ## Note Rt objects from the fits have dim = c(time, vacc_classes, n_par)
  ## We will take Rt for n_vacc_class = 1
  rt <- list(
    Rt_general = restart$parent$rt$Rt_general[, 1, i],
    eff_Rt_general = restart$parent$rt$eff_Rt_general[, 1, i]
  )

  ## Our final object that we will use in the simulations
  ret <- fits[c("dt", "steps_per_day")]
  ret$scenario_step <- ZamCovid:::numeric_date(date_start) * ret$steps_per_day
  ret$date <- date_start
  ret$base <- fits$parameters$base
  ret$pars <- pars
  ret$state <- state
  ret$info <- info
  ret$rt <- rt
  ret$index <- i
  ret
}


create_vaccination_counterfactual <- function(control, base) {

  vacc_scenarios <- names(control$scenarios)[1:15]
  ret <- NULL
  for (i in vacc_scenarios) {
    message(paste("Creating counterfactual for scenario", i))
    doses <- control$scenarios[[i]]$doses
    vaccine_start <- control$scenarios[[i]]$start_date
    n_days <- control$scenarios[[i]]$n_days
    
    tmp <- ZamCovid::vaccine_schedule_historic(
      data = NULL,
      uptake = control$uptake_by_age,
      age_priority = control$age_priority,
      pop_to_vaccinate = base$population$n,
      daily_doses_value = doses,
      weekly_doses = TRUE,
      days_between_doses = control$days_between_doses,
      start = vaccine_start,
      dose_waste = 0,
      n_days = n_days)
    
    ret[[i]] <- tmp
  }
  ret
}


create_healthcare_counterfactual <- function(pars_update, fits, control) {
  
  # Fitted parameters
  pars <- fits$pars
  
  # Scale fitted p_G_D piecewise function for each healthcare scenario
  nms <- names(control$scenarios)[16:24]
  
  # Given dates at which we set p_G_D change points for the fits the resulting
  # fitted p_G_D_step trajectory is 1949 steps long. That is, the last change
  # point is on 2021-05-01 and dt is of 1 / 4. We will reconstruct a p_G_D_step
  # trajectory for each healthcare scenario, by scaling values after the
  # appropriate step. Note we assume the improvement in healthcare is gradually
  # implemented over 21 days from the date in each scenario.
  implement_length <- 21 * 4
  p_G_D_length <- dim(fits$pars[[1]]$p_G_D_step)[1]
  
  for (i in seq_along(nms)) {
    
    nm <- nms[i]
    start_date <- control$scenarios[[nm]]$start_date
    step_start <- start_date * 4
    step_end <- step_start + implement_length
    scale_ifr <- control$scenarios[[nm]]$scale_ifr
    
    scale_vect <- c(rep(1, step_start),
                    seq(1, 1 - scale_ifr, length.out = implement_length),
                    rep(1 - scale_ifr, p_G_D_length - step_end))
    
    pars_update[[nm]] <- pars
    for (j in seq_along(pars_update[[nm]])) {
      pars_update[[nm]][[j]]$p_G_D_step <-
        pars_update[[nm]][[j]]$p_G_D_step * scale_vect
    }
  }
  pars_update
}
