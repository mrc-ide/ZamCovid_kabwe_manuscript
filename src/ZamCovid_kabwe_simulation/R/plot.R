# Get model trajectories into long format tibble
traj_to_long <- function(t, tmp, dates, data) {
  data.frame(
    date = dates,
    state = t,
    mean = colMeans(tmp, na.rm = TRUE),
    lb = matrixStats::colQuantiles(tmp, probs = 0.025, na.rm = TRUE),
    ub = matrixStats::colQuantiles(tmp, probs = 0.975, na.rm = TRUE),
    data = data
  )
}

start_labels <- c("100th day roll-out", "HIC roll-out", "Factual roll-out",
                  "1st wave implementation", "2nd wave implementation",
                  "3rd wave implementation")
rate_labels <- c("HIC capacity", "Upper MIC capacity", "Lower MIC capacity", "LIC capacity",
                 "Factual capacity", "Minimal hc. improv.", "Moderate hc. improv.",
                 "Optimal hc. improv.")

vacc_scenario_labels <- function(df) {
  
  st_lab <- start_labels[1:3]
  rt_lab <- rate_labels[1:5]
  
  df <- df %>%
    mutate(scen_start = as.factor(case_when(
      str_detect(scenario, "cepi") ~ st_lab[1],
      str_detect(scenario, "high_income_") ~ st_lab[2],
      str_detect(scenario, "zambia") ~ st_lab[3])),
      
      scen_rate = as.factor(case_when(
        str_detect(scenario, "_high_income") ~ rt_lab[1],
        str_detect(scenario, "_upper_middle_income") ~ rt_lab[2],
        str_detect(scenario, "_lower_middle_income") ~ rt_lab[3],
        str_detect(scenario, "_low_income") ~ rt_lab[4],
        str_detect(scenario, "_factual") ~ rt_lab[5])))
  
  df$scen_start <- factor(df$scen_start, levels = st_lab)
  df$scen_rate <- factor(df$scen_rate, levels = rt_lab)
  
  df
}


health_scenario_labels <- function(df) {
  
  st_lab <- start_labels[4:6]
  rt_lab <- rate_labels[6:8]
  
  df <- df %>%
    mutate(scen_start = as.factor(case_when(
      str_detect(scenario, "1st") ~ st_lab[1],
      str_detect(scenario, "2nd") ~ st_lab[2],
      str_detect(scenario, "3rd") ~ st_lab[3])),
      
      scen_rate = as.factor(case_when(
        str_detect(scenario, "pessimistic") ~ rt_lab[1],
        str_detect(scenario, "central") ~ rt_lab[2],
        str_detect(scenario, "optimistic") ~ rt_lab[3])))
  
  df$scen_start <- factor(df$scen_start, levels = st_lab)
  df$scen_rate <- factor(df$scen_rate, levels = rt_lab)
  
  df
}


values <- c("royalblue", "orange3", "purple2",
            "green4", "yellow4", "red3", "grey30")
breaks <- c(start_labels, "Fitted")
wave_breaks <- c("2020-09-15", "2021-03-15")


plot_scenario_timeseries <- function(fits, sim, outcome = "deaths_covid",
                                     per_wave = FALSE, which = "vaccination") {
  
  if (which == "vaccination") {
    sim <- sim[1:15]
  } else if (which == "healthcare") {
    sim <- sim[16:24]
  }
  
  data <- fits$fitted
  dates_vect <- as.Date(data$date_string)
  
  # Fitted model trajectory to long format
  if (outcome == "deaths_covid") {
    model_fit <- fits$state["deaths_all_inc", , -1L] - fits$state["base_death_inc", , -1L]
  } else {
    model_fit <- fits$state[outcome, , -1L]
  }
  
  df_fit <- traj_to_long("Fitted", model_fit,
                         dates_vect, rep(NA_real_, length(dates_vect))) %>%
    select(date, fit_mean = mean, fit_lb = lb, fit_ub = ub)
  
  # Simulated trajectories 
  df_sim <- NULL
  for (i in names(sim)) {
    tmp <- sim[[i]]
    sim_dates <- ZamCovid:::numeric_date_as_date(tmp$date)
    
    if (outcome == "deaths_covid") {
      model_sim <- tmp$state["deaths_all_inc", , ] - tmp$state["base_death_inc", , ]
    } else {
      model_sim <- tmp$state[outcome, , ]
    }
    
    tmp <- traj_to_long("Counterfactual", model_sim,
                        sim_dates, rep(NA_real_, length(sim_dates))) %>%
      left_join(data.frame(date = dates_vect), ., by = "date") %>%
      mutate(scenario = i) %>%
      left_join(., df_fit)
    
    df_sim <- rbind(df_sim, tmp)
  }
  
  
  if (outcome == "deaths_covid") {
    gg_title <- "all-cause excess deaths"
  } else {
    data$mean <- NA_real_
    gg_title <- "infections"
  }
  
  
  ## Produce plot 
  if (which == "vaccination") {
    df_sim <- vacc_scenario_labels(df_sim) %>%
      filter(scen_rate != "Factual capacity")
  } else if (which == "healthcare") {
    df_sim <- health_scenario_labels(df_sim) %>%
      filter(scen_rate != "Factual capacity")
  }
  
  ggplot(df_sim, aes(x = date)) +
    geom_line(aes(y = mean, col = scen_start), linewidth = 1, linetype = 1) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = scen_start),
                alpha = 0.5, show.legend = FALSE) +
    geom_line(aes(y = fit_mean, col = "Fitted"), linetype = 3) +
    geom_ribbon(aes(ymin = fit_lb, ymax = fit_ub, fill = "Fitted"),
                alpha = 0.3, show.legend = FALSE) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    labs(x = "", y = paste("Daily new", gg_title)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    scale_color_manual(values = values, breaks = breaks) +
    scale_fill_manual(values = values, breaks = breaks) +
    facet_grid(rows = vars(scen_rate), cols = vars(scen_start)) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 10),
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 10))
}



plot_scenario_cumulative <- function(fits, sim, sample, outcome = "deaths_covid",
                                     which = "vaccination") {
  
  if (which == "vaccination") {
    sim <- sim[1:15]
    tit <- "Counterfactual vaccination scenarios"
  } else if (which == "healthcare") {
    sim <- sim[16:24]
    tit <- "Counterfactual healthcare improvement scenarios"
  }
  

  # Fitted model trajectory to long format
  if (outcome == "deaths_covid") {
    model_fit <- fits$state["deaths_cum_hosp", sample, -1L] +
      fits$state["deaths_cum_comm", sample, -1L]
  } else {
    model_fit <- fits$state[outcome, sample, -1L]
  }
  
  df_fit <- as.data.frame(model_fit[, ncol(model_fit)]) %>%
    `colnames<-`("value") %>%
    mutate(scenario = "Fitted") %>%
    select(scenario, value)
  
  # Simulated trajectories 
  df <- NULL
  for (i in names(sim)) {
    tmp <- sim[[i]]
    
    if (outcome == "deaths_covid") {
      model_sim <- tmp$state["deaths_cum_hosp", , ] + tmp$state["deaths_cum_comm", , ]
    } else {
      model_sim <- tmp$state[outcome, , ]
    }
    
    tmp <- as.data.frame(model_sim[, ncol(model_sim)]) %>%
      `colnames<-`("value") %>%
      mutate(scenario = i) %>%
      select(scenario, value)
    
    tmp$value <- tmp$value / df_fit$value
    
    df <- rbind(df, tmp)
  }
  
  if (outcome == "deaths_covid") {
    gg_title <- "deaths"
  } else {
    gg_title <- "infections"
  }
  

  ## Produce plot 
  if (which == "vaccination") {
    df <- vacc_scenario_labels(df)
  } else if (which == "healthcare") {
    df <- health_scenario_labels(df)
  }
  
  df <- df %>%
    filter(scen_rate != "Factual capacity")
  
  if (which == "vaccination") {
    val <- c("red2", "orange2", "yellow2", "green4")
    brk <- c("LIC", "L-MIC", "U-MIC", "HIC")
    
    df$scen_rate <- factor(df$scen_rate,
                           levels = paste(c("LIC", "Lower MIC",
                                            "Upper MIC", "HIC"), "capacity"),
                           labels = brk)
    df$scen_start <- factor(df$scen_start,
                            levels = c("100th day roll-out", "HIC roll-out",
                                       "Factual roll-out"),
                            labels = c("Apr-2020 (100-days mission)",
                                       "Dec-2020 (HIC roll-out)",
                                       "Apr-2021 (Zambia roll-out)"))
    
    legend_tit <- "Vaccination capacity (average of income setting)"
    
  } else {
    val <- c("turquoise3", "blue", "purple")
    brk <- c("Minimal", "Moderate", "Optimal")
    df$scen_rate <- factor(df$scen_rate,
                           levels = c("Minimal hc. improv.", "Moderate hc. improv.",
                                      "Optimal hc. improv."),
                           labels = brk)
    df$scen_start <- factor(df$scen_start,
                            levels = c("1st wave implementation",
                                       "2nd wave implementation",
                                       "3rd wave implementation"),
                            labels = c("Trans-1st wave start",
                                       "Pre-2nd wave start",
                                       "Pre-3rd wave start"))
    
    legend_tit <- "Improvement in healthcare for severe COVID-19"
  }
  
  ggplot(df) +
    geom_violin(aes(x = scen_rate, y = 1 - value, fill = scen_rate),
                alpha = 0.4) +
    geom_boxplot(width = 0.1, aes(x = scen_rate, y = 1 - value)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                       labels = scales::percent_format()) +
    labs(x = "", y = paste("Averted", gg_title, "(% vs fitted)")) +
    scale_fill_manual(values = val, breaks = brk, name = legend_tit) +
    facet_grid(cols = vars(scen_start)) +
    theme_minimal() +
    theme(legend.title = element_text(),
          legend.position = "top",
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 10),
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 10))
}



plot_infection_incidence <- function(fit, sim, pars, xlim = "2020-07-01") {
  
  dates_fit <- ZamCovid:::numeric_date_as_date(fit$date)[-1L]
  dates_sim <- ZamCovid:::numeric_date_as_date(sim$date)
  pop <- sum(pars[[1]]$N_tot_all)
  
  values <- c("darkgoldenrod3", "purple2")
  breaks <- c("Fitted", "Counterfactual")
  
  # Inferred infection incidence rates x 1,000 population
  model_fit <- fit$state["infections_inc", , -1L]
  model_sim <- sim$counterfactual["infections_inc", , ]
  
  df <- rbind(
    traj_to_long("Fitted", model_fit, dates_fit, rep(NA_real_, length(dates_fit))),
    traj_to_long("Counterfactual", model_sim, dates_sim, rep(NA_real_, length(dates_sim)))) %>%
    select(!data) %>%
    filter(date >= as.Date(xlim))
  
  ggplot(df, aes(x = date)) +
    geom_line(aes(y = (mean / pop) * 1e3, col = state)) +
    geom_ribbon(aes(ymin = (lb / pop) * 1e3,
                    ymax = (ub / pop) * 1e3, 
                    fill = state), alpha = 0.3) +
    scale_color_manual(values = values, breaks = breaks) +
    scale_fill_manual(values = values, breaks = breaks) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    labs(y = "Infection incidence (x 1,000 population)", x = "") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "right",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 12))
  
}


get_immunity_levels <- function(fits, sim, pars) {
  N_age <- pars[[1]][[1]]$N_tot
  N_tot <- sum(N_age)
  
  fit_dates <- ZamCovid:::numeric_date_as_date(fits$date)
  immune_states <- c("immune_S_vacc" ,"immune_R_vacc", "immune_R_unvacc")
  factual <- traj_to_long("Inf. immunity (fit period)",
                          fits$state["immune_R_unvacc", , ],
                          fit_dates, rep(NA_real_, length(fit_dates))) %>%
    select(!data)
  
  out <- NULL
  for (i in names(sim)) {
    tmp <- sim[[i]]
    sim_dates <- ZamCovid:::numeric_date_as_date(tmp$date)
    
    ## Susceptible trajectories
    traj <- tmp$state[immune_states, , ]
    ret <- NULL
    for (j in rownames(traj)) {
      
      nm <- if (j == immune_states[1]) {
        "Vacc. immunity"
      } else if (j == immune_states[2]) {
        "Hybrid immunity"
      } else if (j == immune_states[3]) {
        "Inf. immunity"
      }
      
      tmp_fit <- factual %>%
        filter(date <= min(sim_dates))
      
      x <- traj_to_long(nm, traj[j, , ], sim_dates,
                        rep(NA_real_, length(sim_dates))) %>%
        select(!data) %>%
        rbind(tmp_fit, .) %>%
        mutate(scenario = i) %>%
        mutate_if(is.numeric,  ~. / N_tot)
      
      ret <- rbind(ret, x)
    }
    
    out <- rbind(out, ret)
  }
  out
}


get_cum_inf <- function(fits, sim, pars, sample, factual_only = FALSE) {
  
  inf_nm <- grep("infections_inc_age_", rownames(fits$state), value = TRUE)
  reinf_nm <- grep("reinfections_inc_age_", rownames(fits$state), value = TRUE)
  inf_nm <- setdiff(inf_nm, reinf_nm)
  
  N_tot <- sum(pars[[1]][[1]]$N_tot)
  
  fnx <- function(s) {
    a <- apply(s[inf_nm, , ], c(2:3), sum)
    b <- apply(s[reinf_nm, , ], c(2:3), sum)
    t(apply(a - b, 1, cumsum))
  }
  
  fit_dates <- ZamCovid:::numeric_date_as_date(fits$date)[-1L]
  factual <- fits$state[c(inf_nm, reinf_nm), sample, -1L]
  
  if (factual_only) {
    fnx(factual)
  } else {
    out <- NULL
    for (i in names(sim)) {
      # browser()
      tmp <- sim[[i]]
      sim_dates <- tmp$date
      tmp <- tmp$state[c(inf_nm, reinf_nm), , ]
      
      ret <- factual
      ret[, , sim_dates] <- tmp[, , ]
      
      ret <- traj_to_long("Cum. infected", fnx(ret), fit_dates,
                          rep(NA_real_, length(fit_dates))) %>%
        select(!data) %>%
        mutate_if(is.numeric,  ~. / N_tot) %>%
        mutate(scenario = i)
      
      out <- rbind(out, ret)
    }
    out
  }
}


plot_effective_immunity <- function(fits,  sim, pars, which = "vaccination",
                                    scen_plot = NULL, add_ever_infected = FALSE,
                                    sample = NULL) {
  
  
  if (which == "vaccination") {
    sim <- sim[1:15]
  } else if (which == "healthcare") {
    sim <- sim[16:24]
  }

    
  if (!is.null(scen_plot)) {
    stopifnot(all(scen_plot %in% names(sim)))
    sim <- sim[scen_plot]
  }
  
  df <- get_immunity_levels(fits, sim, pars)
  
  ## Produce plot 
  if (which == "vaccination") {
    df <- vacc_scenario_labels(df)
  } else if (which == "healthcare") {
    df <- health_scenario_labels(df)
  }
  
  df$state <- factor(df$state,
                     levels = c("Vacc. immunity", "Hybrid immunity",
                                "Inf. immunity (fit period)", "Inf. immunity"),
                     labels = c("Vacc. immunity", "Hybrid immunity",
                                "Inf. immunity","Inf. immunity"))
  values <- c("darkgreen", "blue4", "orange2", "black")
  
  p <- ggplot(df, aes(x = date)) +
    geom_area(aes(y = mean, fill = state), alpha = 0.4) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                       labels = scales::percent_format()) +
    scale_x_date(limits = c(as.Date("2020-03-15"), NA),
                 date_breaks = "2 month", date_labels = "%b-%y") +
    scale_fill_manual(values = values, breaks = unique(df$state)) +
    facet_grid(rows = vars(scen_rate), cols = vars(scen_start)) +
    labs(y = "Population immune against infection (%)", x = "") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "top",
          axis.line = element_line(),
          panel.spacing.y = unit(1, "lines"),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 12))
  
  if (add_ever_infected) {
    stopifnot(!is.null(sample))
    
    x <- get_cum_inf(fits, sim, pars, sample)
    
    colnames(x)[3:5] <- paste0(colnames(x)[3:5], "_inf")
    x <- x %>%
      select(!state) %>%
      left_join(df, .)
    
    p <- ggplot(x, aes(x = date)) +
      geom_area(aes(y = mean, fill = state), alpha = 0.4) +
      geom_line(aes(y = mean_inf, col = "Cumulative infected"), alpha = 0.4) +
      geom_ribbon(aes(ymin = lb_inf, ymax = ub_inf, fill = "Cumulative infected"),
                  alpha = 0.2) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                         labels = scales::percent_format()) +
      scale_x_date(limits = c(as.Date("2020-03-15"), NA),
                   date_breaks = "2 month", date_labels = "%b-%y") +
      scale_fill_manual(values = values, breaks = unique(x$state)) +
      scale_color_manual(values = "black", breaks = "Cumulative infected") +
      facet_grid(rows = vars(scen_rate), cols = vars(scen_start)) +
      labs(y = "Population immune against infection (%)", x = "") +
      theme_minimal() +
      theme(legend.title = element_blank(),
            legend.position = "top",
            axis.line = element_line(),
            panel.spacing.y = unit(1, "lines"),
            axis.text.x = element_text(angle = 45, vjust = 0.7, size = 12))
  }
  p
}


plot_vaccinated <- function(sim, pars) {
  
  plot_title <- "Population vaccinated"
  N_age <- pars[[1]][[1]]$N_tot
  N_tot <- sum(N_age)
  N_eligible <- sum(N_age[4:16]) + (N_age[3] / 5) * 3
  
  vacc_class <- list(
    `One dose` = 1L, `Two doses` = 3L)
  get_vacc_class <- function(i, x, d, class) {
    nm <- names(class[i])
    i <- class[[i]]
    x <- apply(x[, i, , drop = FALSE], 3, sum)
    data.frame(
      date = as.Date(d),
      dose = nm,
      value = x)
  }
  
  df <- NULL
  for (i in names(sim)) {
    
    tmp <- sim[[i]]$n_vaccinated
    sim_dates <- ZamCovid:::numeric_date_as_date(sim[[i]]$date)

    tmp_sim <- lapply(seq_along(vacc_class),
                      get_vacc_class, tmp, sim_dates, vacc_class)
    
    ret <- purrr::map_df(tmp_sim, ~as.data.frame(.x), .id = "id") %>%
      select(!id) %>%
      mutate(scenario = i)
    
    df <- rbind(df, ret)
  }
  
  df <- df %>%
    pivot_wider(names_from = dose) %>%
    mutate(`One dose` = `One dose` - `Two doses`) %>%
    pivot_longer(!c(date, scenario), names_to = "dose")

  df <- vacc_scenario_labels(df) %>%
    filter(scen_rate != "Factual capacity")
  
  ggplot(df, aes(x = date)) +
    geom_area(aes(y = value / N_eligible, alpha = dose, fill = scen_start)) +
    geom_hline(yintercept = 0.7, linetype = 3, linewidth = 0.7, col = "red") +
    geom_hline(yintercept = 0.4, linetype = 3, linewidth = 0.7, col = "red") +
    scale_color_manual(values = values, breaks = breaks) +
    scale_fill_manual(values = values, breaks = breaks) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                       labels = scales::percent_format()) +
    scale_alpha_discrete(range = c(0.2, 0.7)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    facet_grid(rows = vars(scen_rate), cols = vars(scen_start)) +
    labs(y = "", x = "") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "top",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 12))
  
}


plot_rt <- function(fits, sim, which = "vaccination") {
  
  if (which == "vaccination") {
    sim <- sim[1:15]
  } else if (which == "healthcare") {
    sim <- sim[16:24]
  }
  
  dates_vect <- as.Date(fits$fitted$date_string)
  model_fit <- rbind(traj_to_long("Rt*", t(fits$Rt_general[-1L, 1, ]),
                                  dates_vect, rep(NA_real_, length(dates_vect))),
                     traj_to_long("Effective Rt*", t(fits$eff_Rt_general[-1L, 1, ]),
                                  dates_vect, rep(NA_real_, length(dates_vect)))) %>%
    select(!data) %>% 
    mutate(source = "Fitted",
           scenario = NA) %>%
    filter(date >= as.Date("2020-03-14"))
  
  model_sim <- NULL
  for (i in names(sim)) {
    
    tmp <- sim[[i]]
    sim_dates <- ZamCovid:::numeric_date_as_date(tmp$date)
    tmp_rt <- t(tmp[["Rt_general"]])
    tmp_eff_rt <- t(tmp[["eff_Rt_general"]])# / 2.33
    tmp <- rbind(traj_to_long("Rt*", tmp_rt, sim_dates, rep(NA_real_, length(sim_dates))),
                 traj_to_long("Effective Rt*", tmp_eff_rt, sim_dates, rep(NA_real_, length(sim_dates)))) %>%
      select(!data) %>%
      mutate(source = "Simulated",
             scenario = i) %>%
      left_join(., model_fit) %>%
      filter(date >= as.Date("2020-03-14"))
    
    tmp_fit <- model_fit %>%
      filter(date < min(tmp$date)) %>%
      mutate(scenario = i)
    
    tmp <- rbind(tmp_fit, tmp)
    
    model_sim <- rbind(model_sim, tmp)
  }
  
  ## Produce plot 
  if (which == "vaccination") {
    model_sim <- vacc_scenario_labels(model_sim) %>%
      filter(scen_rate != "Factual capacity")
  } else if (which == "healthcare") {
    model_sim <- health_scenario_labels(model_sim) %>%
      filter(scen_rate != "Factual capacity")
  }
  
  model_sim$source <- factor(model_sim$source, levels = c("Simulated", "Fitted"))
  
  ggplot(model_sim, aes(x = date)) +
    geom_line(aes(y = mean, col = state, linetype = source)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = state),
                alpha = 0.3) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    facet_grid(rows = vars(scen_rate), cols = vars(scen_start)) +
    geom_hline(yintercept = 1, linetype = 3, col = "red") +
    labs(x = "", y = "Reproduction number (Rt)") +
    theme_minimal() +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.line = element_line(),
          axis.text.x = element_text(size = 10, angle = 45, vjust = 0.7),
          panel.spacing = unit(2, "lines"),
          strip.text = element_text(size = 10))
}


plot_severity <- function(fits, sim, which = "vaccination") {
  
  if (which == "vaccination") {
    sim <- sim[1:15]
  } else if (which == "healthcare") {
    sim <- sim[16:24]
  }
  
  dates_vect <- as.Date(fits$fitted$date_string)
  model_fit <- traj_to_long("Effective IFR", t(fits$ifr[-1L, ]),
                            dates_vect, rep(NA_real_, length(dates_vect))) %>%
    select(!data) %>% 
    mutate(source = "Fitted",
           scenario = NA) %>%
    filter(date >= as.Date("2020-03-14"))
  
  
  model_sim <- NULL
  for (i in names(sim)) {
    tmp <- sim[[i]]
    sim_dates <- ZamCovid:::numeric_date_as_date(tmp$date)
    tmp <- tmp$ifr["ifr", , ]
    
    tmp <- traj_to_long("Effective IFR", tmp, sim_dates, rep(NA_real_, length(sim_dates))) %>%
      select(!data) %>%
      mutate(source = "Simulated",
             scenario = i) %>%
      # left_join(data.frame(date = dates_vect), ., by = "date") %>%
      left_join(., model_fit) %>%
      filter(date >= as.Date("2020-03-14"))
    
    tmp_fit <- model_fit %>%
      filter(date < min(tmp$date)) %>%
      mutate(scenario = i)
    
    tmp <- rbind(tmp_fit, tmp)
    
    model_sim <- rbind(model_sim, tmp)
    
  }
  
  ## Produce plot 
  if (which == "vaccination") {
    model_sim <- vacc_scenario_labels(model_sim) %>%
      filter(scen_rate != "Factual capacity")
  } else if (which == "healthcare") {
    model_sim <- health_scenario_labels(model_sim) %>%
      filter(scen_rate != "Factual capacity")
  }
  
  model_sim <- model_sim %>%
    mutate(source2 = case_when(
      source == "Fitted" ~ source,
      TRUE ~ scen_start))
  
  model_sim$source <- factor(model_sim$source, levels = c("Simulated", "Fitted"))
  
  ggplot(model_sim, aes(x = date)) +
    geom_line(aes(y = mean, col = source2)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = source2),
                alpha = 0.3) +
    geom_hline(yintercept = 0.0025, col = "red", linetype = 3) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.02),
                       labels = scales::percent_format()) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    scale_color_manual(values = values, breaks = breaks) +
    scale_fill_manual(values = values, breaks = breaks) +
    facet_grid(rows = vars(scen_rate), cols = vars(scen_start)) +
    geom_hline(yintercept = 1, linetype = 3, col = "red") +
    labs(x = "", y = "Effective IFR") +
    theme_minimal() +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.line = element_line(),
          axis.text.x = element_text(size = 10, angle = 45, vjust = 0.7),
          panel.spacing = unit(2, "lines"),
          strip.text = element_text(size = 10))
}

