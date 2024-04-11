
create_baseline <- function(restart_date, pars, assumptions) {
  
  pars_info <- pars$info
  pars_vcv <- pars$proposal
  
  region <- "kabwe"
  message(sprintf("Creating baseline for '%s'", region))
  pars_info <- pars_info[pars_info$region == region | is.na(pars_info$region), ]
  pars_info <- setNames(pars_info$initial, pars_info$name)
  
  ## 1. Set up basic variables and read in data ----
  date <- as.Date("2021-09-30") # last date of fitting period
  dt <- 0.25
  vaccine_days_to_effect_dose1 <- 0
  vaccine_days_to_effect_dose2 <- 7
  vaccine_days_to_effect <- c(vaccine_days_to_effect_dose1,
                              vaccine_days_to_effect_dose2)
  
  uptake_by_age <- read.csv("data/vaccine_uptake.csv", row.names = "group")
  vaccine_efficacy <- read_csv("data/vaccine_efficacy.csv")
  progression_data <- read_csv("data/progression_data.csv")
  severity_data <- read_csv("data/support_severity.csv")
  population <- read_csv("data/population.csv")[, c("age", region)] %>%
    `colnames<-`(., c("age", "n"))
  
  
  
  ## 2. Set up Kabwe population parameters ----
  
  if (region == "kabwe") {
    
    # We don't have population breakdown from Kabwe.
    # The latest census (https://www.zamstats.gov.zm/2022-census/reports/#)
    # reports the total population of Kabwe district at 299,206 in 2022, compared
    # to 202,360 in 2010, with an average annual growth rate in between at 3.3%.
    # We'll assume the same age-distribution as reported nationally.
    g <-  1 - 0.033
    n_kabwe_2022 <- 299206
    n_kabwe_2020 <- round(n_kabwe_2022 * g^2)
    population$n <- round((population$n / sum(population$n)) * n_kabwe_2020)
    
    
    # Now set assumptions of what proportion of baseline deaths are observed
    # A crude empiric estimation is 0.567 unobserved deaths comparing
    # available timeseries vs linelist in 2020.
    # Total deaths from time series in 2020 and 2021 yield a crude mortality
    # rate of 6.54076 and 7.126761, respectively, which tallies up with
    # World Bank projections for the country
    # https://data.worldbank.org/indicator/SP.DYN.CDRT.IN?locations=ZM
    # We'll therefore treat all-cause mortality data as absolute and only
    # account for non-covid deaths being those predicted by the linear model
    # below, as per Ghafari et al. and McCabe et al.
    
    #### NOTE: this code is commented out give the input data for historic
    ####       monthly deaths as plotted in Figure S2 of the technical
    ####       supplement to this publication are not publicly available.
    ####       Access to these data can be granted upon request by contacting
    ####       the authors at p.perez-guzman@imperial.ac.uk
    ####       We will, instead, read in the calculated monthly baseline death
    ####       rate to input into the model. If running this task with raw
    ####       historic deaths data, uncoment lines 64-76 and comment out lines
    ####       77-86.
    # deaths_observed <- 0.567
    # if (assumptions == "observed_low") {
    #   deaths_observed <- deaths_observed * 0.85
    # } else if (assumptions == "observed_high") {
    #   deaths_observed <- deaths_observed * 1.15
    # }
    # 
    # historic_deaths <- infer_baseline_deaths(historic_deaths, date,
    #                                          inflate = 1 / deaths_observed)
    # base_death_date <-
    #   ZamCovid:::numeric_date(historic_deaths$expected_deaths$date)
    # base_death_date[1] <- 0
    # base_death_value <- round(historic_deaths$expected_deaths$rate)
    baseline_deaths <- read_csv("data/baseline_deaths.csv")
    base_death_date <- ZamCovid:::numeric_date(baseline_deaths$date)
    base_death_date[1] <- 0
    if (assumptions == "observed_low") {
      base_death_value <- baseline_deaths$observed_low
    } else if (assumptions == "observed_high") {
      base_death_value <- baseline_deaths$observed_high
    } else {
      base_death_value <- baseline_deaths$central
    }

    ## Cross-immunity assumptions
    # We have a single-strain model that explicitly accounts for protection vs
    # re-infection after recovery from infection. As a proxy of VOC emergence
    # the value cross_immunity (% protection conferred) can be time-varying.
    # Beta was first detected in Zambia on 2020-12-16, but cases had already
    # been increasing from early December
    # https://www.cdc.gov/mmwr/volumes/70/wr/mm7008e2.htm
    # Gill et al. (https://bmjopen.bmj.com/content/12/12/e066763) estimated
    # Beta detection in Lusaka (mortuaty) peaked in January 2021 and Delta
    # in June 2021.
    # Lastly, no reliable data for Omicron emergence. News feed indicates 
    # first case detected on 2021-12-04
    # We will assume same values of protections as in England Perez-Guzman et al.
    # https://www.medrxiv.org/content/10.1101/2023.02.10.23285516v2
    cross_immunity_date <-
      c(0, ZamCovid:::numeric_date(c("2020-12-01", "2021-01-31",
                                     "2021-03-15", "2021-06-30",
                                     "2021-11-20", "2021-12-31")))
    cross_immunity_value <- c(1, 1, 0.95, 0.95, 0.85, 0.85, 0.25)
    
    # Now, set target p_G_D informed by IFR estimates with seroreversion
    # from Brazeau et al. https://www.nature.com/articles/s43856-022-00106-7
    brazeau_ifr <- c(0, 0.0001, 0.0001, 0.0002, 0.0003, 0.0004, 0.0006, 0.001,
                     0.0016, 0.0024, 0.0038, 0.0059, 0.0092, 0.0143, 0.0223,
                     0.0347, 0.0541, 0.0843, 0.164)
    # Their estimates are for 19 age brackets; we only have 16 so will aggregate
    # 75_plus weighting by age distribution using the breakdown of 75+'s as per
    # population estimates from US IDB (https://www.census.gov)
    n_75_plus_brackets <- c(99000, 53004, 21234, 5617)
    
    ifr_75_plus <- data.frame(n = n_75_plus_brackets) %>%
      mutate(ifr = (tail(brazeau_ifr, 4) * n))
    ifr_75_plus <- sum(ifr_75_plus$ifr) / sum(ifr_75_plus$n)
    
    target_p_G_D <- data.frame(p_C = t(severity_data[1, 2:17])) %>%
      `colnames<-`("p_C") %>%
      mutate(ifr = c(brazeau_ifr[1:15], ifr_75_plus),
             p_G_D = ifr / p_C)
    
    severity_data[severity_data$Name == "p_G_D", 2:length(severity_data)] <- 
      if (assumptions == "ifr_low") {
        target_p_G_D$p_G_D * 0.9
      } else if (assumptions == "ifr_high") {
        target_p_G_D$p_G_D * 1.1
      } else {
        target_p_G_D$p_G_D
      }
    
    ## Lastly, set assumptions of infection-induced immunity waning and
    ## seroreversion/
    
    # Immunity waning assumption based on estimate of 24.7% remaining effectively
    # protected at 12 months in Bobrovitz et al.
    # https://linkinghub.elsevier.com/retrieve/pii/S1473309922008015
    imm_waning <- data.frame(parameter = "gamma_R", value = 1 / (2 * 365))
    if (assumptions == "imm_waning_slow") {
      imm_waning$value <- 1 / (2.5 * 365)
    } else if (assumptions == "imm_waning_fast") {
      imm_waning$value <- 1 / (1.5 * 365)
    }
    progression_data <- rbind(progression_data, imm_waning)
    
    # Krutikov et al. 2022 report median time to sero-rev with Abbott 242.5 days,
    # but their population was CHR and CHW.
    # (https://www.thelancet.com/journals/lanhl/article/PIIS2666-7568(21)00282-8/fulltext)
    serorev <- data.frame(parameter = "gamma_sero_pos", value = 1 / 242.5)
    if (assumptions == "serorev_fast") {
      serorev$value <- 1 / (242.5 * 0.8)
    } else if (assumptions == "serorev_slow") {
      serorev$value <- 1 / (242.5 * 1.2)
    }
    progression_data <- rbind(progression_data, serorev)
  }
  
  
  ## 3. Set-up other basic model parameters ----
  # Beta change points - A vector of date (strings) for the beta parameters.
  # We are agnostic as to the effect of official NPI dates at specific dates
  # This is thus just a vector of monthly dates
  beta_date <- as.character(
    seq.Date(as.Date("2020-03-01"), as.Date(date), "2 weeks"))
  beta_names <- sprintf("beta%d", seq_along(beta_date))
  
  # Set of parameters that will be fitted for each model type
  to_fit_all <- c(
    # direct
    "start_date", beta_names,
    # severity
    "p_G_D", "alpha_D", "mu_D_1", "mu_D_2"
  )
  stopifnot(setequal(to_fit_all, names(pars_info)))
  
  
  ## Initial seeding: seed 10 per million over 1 day (4 steps)
  seed_size <- 20e-6 * sum(population[, "n"])
  seed_pattern <- c(1, 1, 1, 1)
  
  
  # For internal check: the assumed serial interval of Wildtype should be ~ 5.2 days
  # Calculation adapted from Svensson
  # https://www2.math.su.se/matstat/reports/seriea/2005/rep14/report.pdf
  # so SI = E + (P^2 + P·C_1 + C_1^2) / (P + C_1)
  
  # We ensure that mean[T_I_C_1 + T_I_C_2] is unchanged
  # In the abscence of variant data, we will assume a decreasing SI with the
  # sequential replacement of Wildtype by Beta and Beta by Delta. We assume
  # same dates as for imperfect cross-immunity.
  # is replaced by Beta and Delta 
  rel_si_wildtype <- 1
  rel_si_beta <- 0.94
  rel_si_delta <- 0.87
  rel_si_omicron <- 0.75
  
  mean_E <- 0.816
  mean_P <- 1.7808
  mean_C_1 <- 2.2684
  mean_C_2 <- 1.9716
  
  rel_C_2_wildtype <- (mean_C_2 + (1 - rel_si_wildtype) * mean_C_1) / mean_C_2
  rel_C_2_beta <- (mean_C_2 + (1 - rel_si_beta) * mean_C_1) / mean_C_2
  rel_C_2_delta <- (mean_C_2 + (1 - rel_si_delta) * mean_C_1) / mean_C_2
  rel_C_2_omicron <- (mean_C_2 + (1 - rel_si_omicron) * mean_C_1) / mean_C_2
  
  rel_gamma_wildtype <- c(1 / rel_si_wildtype, 1 / rel_C_2_wildtype)
  rel_gamma_beta <- c(1 / rel_si_beta, 1 / rel_C_2_beta)
  rel_gamma_delta <- c(1 / rel_si_delta, 1 / rel_C_2_delta)
  rel_gamma_omicron <- c(1 / rel_si_omicron, 1 / rel_C_2_omicron)
  
  rel_gammas <- list(wildtype = rel_gamma_wildtype, beta = rel_gamma_beta,
                     delta = rel_gamma_delta, omicron = rel_gamma_omicron,
                     date = cross_immunity_date)
  
  mean_SI_wildtype <-
    (2 / mean_E) + (mean_P^2 + mean_P * mean_C_1 + mean_C_2^2) / (mean_P + mean_C_1)
  
  ## Set serology assay sensitivity assumptions
  sens_and_spec <- ZamCovid::ZamCovid_parameters_sens_and_spec()
  if (assumptions == "sero_sens_low") {
    sens_and_spec$sero_sensitivity <- 0.754
  } else if (assumptions == "sero_sens_high") {
    sens_and_spec$sero_sensitivity <- 0.99
  }
  
  
  
  ## 4. Set-up vaccination parameters and assumptions ----
  
  ## TODO: keeping all rel_parameters and n_doses for now, but implementing
  ## schedule with no doses available. Need to instead implement a function to
  ## inflate states in the simulation.
  daily_doses_value <- 0
  vaccine_eligibility_min_age <- 12
  mean_days_between_doses <- 7 * 11 # second dose starts 12 weeks after first
  mean_time_to_waned <- 24 * 7 # assume exponential with mean 24 weeks
  time_to_dose_1_effect <- 7 * 3 # assume exponential with mean 3 weeks
  
  # Compartment progression rates
  # These are progression rates OUT of the specified compartment, where
  # zeroes mean movement is controlled by vaccination schedule rather than
  # a rate parameter and waned is an absorbing state
  vaccine_progression_rate <- c(0,                         # unvaccinated 
                                1 / time_to_dose_1_effect, # first dose no effect
                                0,                         # first dose full effect
                                1/ mean_time_to_waned,     # second dose
                                0)                         # waned
  
  # Proportion of vaccines by type and age
  prop_pfizer_by_age <- get_vacc_doses_age()
  
  # Calculate average vaccine efficacy and rel_severity
  average_vacc_efficacy <-
    calculate_average_vacc_efficacy(vaccine_efficacy, prop_pfizer_by_age)
  rel_severity <- lapply(average_vacc_efficacy, function(e)
    get_vaccine_conditional_prob(e$death, e$severe_disease, e$disease,
                                 e$infection, e$transmission))
  
  # Set VE assumption values
  if (assumptions == "ve_high") {
    rel_severity <- rel_severity$ve_high
  } else if (assumptions == "ve_low") {
    rel_severity <- rel_severity$ve_low
  } else {
    rel_severity <- rel_severity$central
  }
  
  # Note that vaccine_uptake[i, j] is proportional uptake of dose j for group i 
  n_doses <- 2
  # Downstream, we'll simulate vaccination from the first Monday after vaccines
  # were first made available
  dose_start_dates <- c("2020-12-14", "2020-12-14")
  vaccine_index_dose2 <- 3L
  age_priority <- list(16, 15, 14, 13, 12, 11, 9:10, 7:8, 1:6)
  
  vaccine_uptake <- 
    array(uptake_by_age$central, c(length(uptake_by_age$central), n_doses))
  
  ## TODO: this is implementing a vaccination schedule with NO doses at the moment.
  ## this is a hack to allow ejecting restart objects with the right dimension
  ## for downstream simulation. Need to instead implement a function to inflate
  ## model states to include vaccination, at the simulation stages.
  vaccination <-
    vaccination_schedule(date, region, vaccine_uptake,
                         vaccine_days_to_effect, NULL,
                         n_doses, dose_start_dates, age_priority,
                         population, daily_doses_value)
  
  vaccine_schedule_real <- vaccination$schedule
  vaccine_schedule_effect <- shift_doses(vaccine_schedule_real,
                                         vaccine_days_to_effect)
  
  
  ## 5. Create list with all parameters to run the model -----
  ret <- list(
    date = date,
    region = region,
    population = population,
    restart_date = restart_date,
    
    beta_date = ZamCovid:::numeric_date(beta_date),
    beta_names = beta_names,
    severity_data = severity_data,
    progression_data = progression_data,
    sens_and_spec = sens_and_spec,
    seed_size = seed_size,
    seed_pattern = seed_pattern,
    base_death_date = base_death_date,
    base_death_value = unname(base_death_value),
    cross_immunity_date = cross_immunity_date,
    cross_immunity_value = cross_immunity_value,
    
    rel_severity = rel_severity,
    rel_gammas = rel_gammas,
    
    vaccine_eligibility_min_age = vaccine_eligibility_min_age,
    vaccine_progression_rate = vaccine_progression_rate,
    vaccine_schedule = vaccine_schedule_real,
    vaccine_schedule_effect = vaccine_schedule_effect,
    vaccine_uptake = vaccine_uptake,
    vaccine_mean_days_between_doses = mean_days_between_doses,
    vaccine_index_dose2 = vaccine_index_dose2,
    vaccine_days_to_effect = vaccine_days_to_effect,
    n_doses = n_doses)
  
  message("  - Creating transformation function")
  tr <- make_transform(ret)
  message("  - Testing transformation function")
  p <- tr(pars_info)
  message("  - Testing creating model with transformed parameters")
  m <- ZamCovid::ZamCovid$new(p, 0, 1)
  ret
}