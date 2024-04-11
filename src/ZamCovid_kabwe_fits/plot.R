# Get model trajectories into long format tibble
traj_to_long <- function(t, tmp, dates, data) {
  data.frame(
    date = dates,
    state = t,
    mean = colMeans(tmp),
    lb = matrixStats::colQuantiles(tmp, probs = 0.025),
    ub = matrixStats::colQuantiles(tmp, probs = 0.975),
    data = data
  )
}


calc_seropos <- function(tmp, dates, name) {
  
  pos <- tmp[, grep("pos", colnames(tmp))]
  tot <- tmp[, grep("tot", colnames(tmp))]
  pos[is.na(pos)] <- 0
  tot[is.na(tot)] <- 0
  
  data.frame(date = dates) %>%
    cbind(., as.data.frame(Hmisc::binconf(pos, tot))) %>%
    `colnames<-`(c("date", "data_mean", "data_lb", "data_ub")) %>%
    mutate(state = name)
}


plot_fit_traces <- function(samples) {
  
  if (is.null(samples$chain)) {
    n_chains <- 1L
  } else {
    n_chains <- length(unique(samples$chain))
  }
  
  cols <- rev(viridisLite::viridis(n_chains))
  
  reorder_beta <- function(nms) {
    i <- grep("^beta[0-9]+$", nms)
    j <- order(as.integer(sub("^beta", "", nms[i])))
    k <- seq_along(nms)
    k[i] <- i[j]
    k
  }
  
  i <- reorder_beta(colnames(samples$pars_full))
  pars <- samples$pars_full[, i]
  nms <- colnames(pars)
  probs <- samples$probabilities_full
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  new_grid <- function(n, title) {
    par(mfrow = rep(ceiling(sqrt(n + 1)), 2),
        mar = c(3, 3, 2, 1),
        mgp = c(2, 0.5, 0),
        oma = c(1, 1, 1 + as.integer(title), 1))
  }
  
  plot_traces1 <- function(p, name) {
    traces <- matrix(p, ncol = n_chains)
    ess <- coda::effectiveSize(coda::as.mcmc(traces))
    
    if (name == "log_likelihood") {
      main <- ""
    } else {
      main <- sprintf("ess = %s", round(sum(ess)))
    }
    matplot(traces, type = "l", lty = 1,
            xlab = "Iteration", bty = "n",
            ylab = name, col = cols,
            main = main,
            font.main = 1)
  }
  
  new_grid(length(nms), FALSE)
  for (nm in nms) {
    plot_traces1(samples$pars_full[, nm], nm)
  }
  plot_traces1(probs[, "log_likelihood"], "log_likelihood")
}


plot_serology <- function(dat, which = NULL, labels = NULL, immunity = FALSE) {
  
  data <- dat$fit$data$full
  if (is.null(which)) {
    which <- c("all", "over15", "15_19", "20_29",
               "30_39", "40_49", "50_plus")
    labels <- c("All", "≥15-yo.", "15 to 19",
                "20 to 29", "30 to 39", "40 to 49",
                "50 and older")
  }
  
  sero_cols <- grep("sero", colnames(data), value = TRUE)
  dates <- ZamCovid:::numeric_date_as_date(dat$fit$samples$trajectories$date)[-1]
  data <- data[, c("date", sero_cols)] %>%
    filter(as.Date(date) %in% dates)
  
  values <- c("cadetblue", "darkviolet","darkgoldenrod3", "chartreuse4", "orange")
  breaks <- c("Model", "Ever infected", "Data (fitted)", "Data", 
              "Immune")

  
  data_plot <- NULL
  for (w in which) {
    keep <- grep(w, colnames(data), value = TRUE)
    tmp <- calc_seropos(data[, keep], dates, w)
    
    data_plot <- rbind(data_plot, tmp)
  }
  
  states <- dat$fit$samples$trajectories$state
  nms <- rownames(states)
  
  p <- dat$fit$samples$predict$transform(dat$fit$samples$pars[1, ])
  sens <- p$sero_sensitivity
  spec <- p$sero_specificity
  
  df <- NULL
  for (t in which) {
    n <- p[[paste0("N_tot_", t)]]
    sero_pos <- states[paste0("sero_pos_", t), , -1]
    
    tmp_sero <- ((sens *  sero_pos + (1 - spec) * (n - sero_pos)) / n) - 0.001
    cum_inf <- get_cum_inf_age(states, t) / n
    
    ret_sero <- traj_to_long(t, tmp_sero, dates, rep(NA_real_, dim(tmp_sero)[2])) %>%
      select(!data)
    ret_cum_inf <- traj_to_long(t, cum_inf, dates, rep(NA_real_, dim(cum_inf)[2])) %>%
      select(!data) %>%
      rename(cum_inf_mean = mean, cum_inf_lb = lb, cum_inf_ub = ub)
    
    ret <- left_join(ret_sero, ret_cum_inf)
    df <- rbind(df, ret)
  }
  
  df <- left_join(df, data_plot) %>%
    mutate(col = case_when(state == "over15" ~ "Data",
                           TRUE ~ "Data (fitted)"))
  df$state <- factor(df$state, levels = unique(df$state), labels = labels)
  
  p <- df %>% 
    filter(state != "All") %>%
    ggplot(., aes(x = date)) +
    geom_line(aes(y = mean, col = "Model"), linewidth = 0.5) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = "Model"), alpha = 0.4, show.legend = FALSE) +
    geom_line(aes(y = cum_inf_mean, col = "Ever infected"), alpha = 0.7, linewidth = 0.5) +
    geom_ribbon(aes(ymin = cum_inf_lb, ymax = cum_inf_ub, fill = "Ever infected"),
                alpha = 0.4, show.legend = FALSE) +
    geom_point(aes(y = data_mean, col = col), size = 1, alpha = 0.5) +
    geom_errorbar(aes(ymin = data_lb, ymax = data_ub, 
                      col = col), linewidth = 0.5, alpha = 0.5) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 1),
                       labels = scales::percent_format()) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y",
                 limits = c(as.Date("2020-02-15"), as.Date("2021-09-30"))) +
    facet_wrap(~state) +
    labs(x = "", y = "Seropositive (% age group)") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "top",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
  if (length(unique(df$col)) == 2) {
    
    p <- p +
      scale_color_manual(values = values, breaks = breaks,
                         guide =
                           guide_legend(override.aes = list(linetype = rep(1, 4),
                                                            size = rep(2, 4),
                                                            shape = c(NA, NA, 16, 16)))) +
      scale_fill_manual(values = values, breaks = breaks)
    
  } else {
    p <- p +
      scale_color_manual(values = values, breaks = breaks) +
      scale_fill_manual(values = values, breaks = breaks)
    
  }
    
  
  if ("all" %in% which) {
    
    p_all <- df %>% 
      filter(state == "All") %>%
      ggplot(., aes(x = date)) +
      geom_line(aes(y = mean, col = "Model")) +
      geom_ribbon(aes(ymin = lb, ymax = ub, fill = "Model"), alpha = 0.4, show.legend = FALSE) +
      geom_line(aes(y = cum_inf_mean, col = "Ever infected")) +
      geom_ribbon(aes(ymin = cum_inf_lb, ymax = cum_inf_ub, fill = "Ever infected"), alpha = 0.4, show.legend = FALSE) +
      geom_point(aes(y = data_mean, col = col), size = 1.5, alpha = 1) +
      geom_errorbar(aes(ymin = data_lb, ymax = data_ub, 
                        col = col), linewidth = 0.3, alpha = 0.5) +
      scale_y_continuous(expand = c(0, 0),
                         limits = c(0, 0.2),
                         labels = scales::percent_format()) +
      scale_x_date(date_breaks = "1 month", date_labels = "%b-%y",
                   limits = c(as.Date("2020-03-15"), as.Date("2020-08-15"))) +
      scale_color_manual(values = values, breaks = breaks) +
      scale_fill_manual(values = values, breaks = breaks) +
      labs(x = "", y = "Seropositive (% all population)") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.line = element_line(),
            axis.text.x = element_text(angle = 45, vjust = 0.7, size = 10))
  }
  
  if (length(which) == 1) {
    
    if (immunity) {
      N_tot <- sum(dat$fit$parameters$base$population$n)
      immune <- traj_to_long("Immune",
                             samples$trajectories$state["immune_R_unvacc", , -1L] / N_tot,
                             dates, rep(NA_real_, length(dates))) %>%
        select(!data)
      
      p <- p +
        geom_line(aes(y = immune$mean, col = "Immune"), alpha = 0.7, linewidth = 0.5) +
        geom_ribbon(aes(ymin = immune$lb, ymax = immune$ub, fill = "Immune"),
                    alpha = 0.4, show.legend = FALSE)
      
    }
    
    p <- add_voc_labs(p, labs = TRUE, ymax = 1)
    p + 
      labs(y = paste0("Sero-positive (% ", labels,")"), title = "") +
      theme(strip.text = element_blank(),
            legend.title = element_blank(),
            legend.position = "right")
  } else {
    p_all / p
  }
}


add_voc_labs <- function(p, ymax, labs = FALSE) {
  
  dates <- as.Date(c("2020-12-16", "2021-06-01"))
  names(dates) <- paste(c("Beta", "Delta"), "det.")
  xmax <- max(p$data$date)
  lab_pos_1 <- dates[1] + (dates[2] - dates[1]) / 2
  lab_pos_2 <- dates[2] + (xmax - dates[2]) / 2
  
  vacc_dates <- as.Date(c("2020-12-08", "2021-04-20"))
  
  p <- p +
    annotate(geom = "rect", xmin = dates[1], xmax = dates[2],
             ymin = 0, ymax = ymax, alpha = 0.1, fill = "#723e8a") +
    annotate(geom = "rect", xmin = dates[2], xmax = xmax,
             ymin = 0, ymax = ymax, alpha = 0.1, fill = "#f76002") +
    geom_vline(xintercept = vacc_dates, linetype = 3, linewidth = 0.7)
  
  if (labs) {
    p +
      annotate("label", x = c(lab_pos_1, lab_pos_2),
               y = ymax * 0.95, size = 2.8, label = names(dates)) +
      annotate("text", x = vacc_dates, y = ymax * c(0.7, 0.82), size = 2.6,
               label = c("1st vax World", "1st vax Zambia"))
  } else {
    p
  }
}


plot_deaths <- function(dat, week_only = FALSE) {
  
  data <- dat$fit$data$fitted
  which <- paste0("deaths_", c("all", "hosp"))
  values <- c("cadetblue", "darkgoldenrod3", "chartreuse4")
  breaks <- c("Model", "Data (fitted)", "Data")
  
  data_plot <- data[, c("date", "date_string", which)]
  data_plot$deaths_comm_inc <- data_plot$deaths_all - data_plot$deaths_hosp
  colnames(data_plot) <- c("date", "date_string", "deaths_all_inc",
                           "deaths_hosp_inc", "deaths_comm_inc")
  
  states <- dat$fit$samples$trajectories$state
  nms <- rownames(states)
  
  traj <- c("deaths_all_inc")
  dates <- as.Date(data_plot$date_string)
  df <- NULL
  for (t in traj) {
    ret <- traj_to_long(t, states[t, , -1], dates, data_plot[, t])
    df <- rbind(df, ret)
  }
  
  df$state <- factor(df$state, levels = unique(df$state))
  df <- df %>%
    filter(!is.na(data))
  
  ylim <- max(df$ub) * 1.25
  
  daily <- ggplot(df, aes(x = date)) +
    geom_line(aes(y = mean, col = "Model")) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = "Model"), alpha = 0.4, show.legend = FALSE) +
    geom_point(aes(y = data, col = "Data (fitted)"), size = 0.7, alpha = 0.9) +
    geom_point(aes(y = NA_real_, col = "Data"), size = 0.7, alpha = 0.9) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y",
                 limits = c(as.Date("2020-02-15"), NA)) +
    scale_color_manual(values = values, breaks = breaks,
                       guide =
                         guide_legend(override.aes = list(linetype = c(1, NA, NA),
                                                          size = rep(2, 3),
                                                          shape = c(NA, 16, 16)))) +
    scale_fill_manual(values = values, breaks = breaks) +
    labs(x = "", y = "Daily deaths (all-cause)") +
    theme_minimal() +
    theme(legend.position = "top",
          legend.title = element_blank(),
          axis.line = element_line(),
          axis.text.x = element_blank())
  
  weekly <- df %>%
    mutate(date = lubridate::floor_date(date, "week", week_start = 3)) %>%
    select(!state) %>%
    pivot_longer(!c(date)) %>%
    group_by(date, name) %>%
    summarise(value = sum(value, na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider()
  weekly <- weekly[-c(1, nrow(weekly)), ]
  
  ymax <- max(weekly$ub, weekly$data) * 1.05
  weekly <- ggplot(weekly, aes(x = date)) +
    geom_line(aes(y = mean, col = "Model"), linewidth = 0.5) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = "Model"), alpha = 0.4, show.legend = FALSE) +
    geom_point(aes(y = data, col = "Data"), size = 1, alpha = 0.9) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y",
                 limits = c(as.Date("2020-02-15"), NA)) +
    scale_color_manual(values = values, breaks = breaks,
                       guide =
                         guide_legend(override.aes = list(linetype = c(1, NA),
                                                          size = rep(2, 2),
                                                          shape = c(NA, 16)))) +
    scale_fill_manual(values = values, breaks = breaks) +
    labs(x = "", y = "Weekly deaths (all-cause)") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 10))
  
  
  if (week_only) {
    weekly <- add_voc_labs(weekly, ymax = ymax)
    weekly + labs(y = "Weekly deaths (all-cause)", title = "") + #,
                  # caption = "* Seropositivity data fitted to age-specific timeseries; all-cause deaths fitted to daily timeseries.
                  # For clarity of visualisation, aggregated data are presented.
                  # Shaded area around model trajectory represents 95% Credible Interval") +
      theme(legend.position = "left",
            legend.title = element_blank())
  } else {
    daily / weekly
  }
}


plot_deaths_disag <- function(dat, plot_age = TRUE) {
  
  which <- c("deaths_all_inc", "base_death_inc", "deaths_covid_inc",
             paste0("D_", seq(0, 75, 5)))
  dates <- ZamCovid:::numeric_date_as_date(dat$fit$samples$trajectories$date)[-1]
  data <- rep(NA_real_, length(dates))
  
  states <- dat$fit$samples$trajectories$state
  
  df <- NULL
  for (t in which) {
    
    if (t == "deaths_covid_inc") {
      tmp <- states["deaths_comm_inc", , -1] + states["deaths_hosp_inc", , -1]
    } else {
      tmp <- states[t, , -1]
    }
    ret <- traj_to_long(t, tmp, dates, data) %>% select(!data)
    df <- rbind(df, ret)
  }
  
  
  deaths_agg <- df %>%
    filter(state %in% c("base_death_inc", "deaths_covid_inc")) %>%
    select(date, state, mean) %>%
    group_by(date, state) %>%
    summarise(n = sum(mean)) %>%
    mutate(percent = n / sum(n)) %>%
    filter(state == "deaths_covid_inc")
  
  agg <- ggplot(deaths_agg, aes(x = date)) +
    geom_area(aes(y = percent, fill = state), alpha = 0.6) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                       labels = scales::percent_format()) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    scale_fill_manual(values = "brown4") +
    labs(x = "", y = "Excess deaths (% all deaths)") +
    theme_minimal() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 12))
  
  deaths_age <- df %>%
    filter(!state %in% c("deaths_all_inc", "base_death_inc", "deaths_covid_inc")) %>%
    select(date, state, mean) %>%
    group_by(date, state) %>%
    summarise(n = sum(mean)) %>%
    mutate(percent = n / sum(n))
  
  deaths_age$state <- factor(deaths_age$state,
                             levels = paste0("D_", seq(0, 75, 5)),
                             labels = c(paste0(seq(0, 70, 5),
                                               " to ", seq(4, 74, 5)),
                                        "75+"))
  
  
  age <- ggplot(deaths_age, aes(x = date)) +
    geom_area(aes(y = percent, fill = state), alpha = 0.75) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                       labels = scales::percent_format()) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y",
                 limits = c(as.Date("2020-01-01"), NA)) +
    scale_fill_viridis_d(direction = 1) +
    labs(x = "", y = "Deaths by age (% daily Covid-19 deaths)") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 10))
  
  if (plot_age) {
    ((agg + theme(axis.text.x = element_blank())) / age) +
      plot_layout(heights = c(0.5, 1))
  } else {
    add_voc_labs(agg, ymax = 1)
  }
}



plot_rt <- function(dat, vocs = FALSE, which = "eff_rt") {
  
  stopifnot(which %in% c("eff_rt", "rt", "all"))

  rt <- dat$fit$rt$Rt_general[, 1, ]
  eff_rt <- dat$fit$rt$eff_Rt_general[, 1, ]
  dates <- ZamCovid:::numeric_date_as_date(dat$fit$rt$time[, 1] / 4)
  
  df_to_long <- function(df, which) {
    data.frame(
      date = dates,
      mean = rowMeans(df),
      lb = matrixStats::rowQuantiles(df, probs = 0.025),
      ub = matrixStats::rowQuantiles(df, probs = 0.975),
      type = which
    )
  }
  
  if (which == "all") {
    df <- rbind(df_to_long(rt, "Rt"),
                df_to_long(eff_rt, "Effective Rt"))
  } else if (which == "rt") {
    df <- df_to_long(rt, "Rt") %>%
      mutate(type = "Model    ")
  } else if (which == "eff_rt") {
    df <- df_to_long(eff_rt, "Effective Rt") %>%
      mutate(type = "Model    ")
  }
  
  values <- c("red4", "lightblue4", "cadetblue")
  breaks <- c("Rt", "Effective Rt", "Model    ")
  
  p <- ggplot(df, aes(date)) +
    geom_line(aes(y = mean, col = type), linewidth = 0.5) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = type), alpha = 0.3) +
    scale_color_manual(values = values, breaks = breaks) +
    scale_fill_manual(values = values, breaks = breaks) +
    geom_hline(yintercept = 1, linewidth = 0.5, col = "red", linetype = 3) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 10, 0.5), limits = c(0, 3.25)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    labs(x = "", y = expression(R["t"])) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
  if (vocs) {
    add_voc_labs(p, ymax = 3.25)
  } else {
    p
  }
}


plot_severity <- function(dat, age = TRUE, xmin = "2020-04-01") {
  
  ifr <- get_severity(dat$fit$severity, "ifr")
  ifr[ifr$date < as.Date(xmin), c("mean", "lb", "ub")] <- NA_real_
  
  p1 <- ggplot(ifr, aes(date, mean)) +
    geom_line(col = "cadetblue", linewidth = 0.5) +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.3, fill = "cadetblue") +
    geom_hline(yintercept = 0.0012, linetype = 3, col = "red", linewidth = 0.5) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.015),
                       breaks = seq(0, 1, 0.0025),
                       labels = scales::percent_format()) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    labs(y = expression(IFR["t"]), x = "") +
    theme_minimal() +
    theme(axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
  ifr_age <- get_severity(dat$fit$severity, "ifr", TRUE) %>%
    filter(date >= as.Date(xmin)) %>%
    pivot_longer(!c(date, age), names_to = "estimate") %>%
    group_by(age, estimate) %>%
    summarise(value = mean(value)) %>%
    pivot_wider(names_from = estimate) %>%
    mutate(age = factor(age, levels = seq(0, 75, 5)))
  
  p2 <- ggplot(ifr_age, aes(y = age, col = age)) +
    geom_pointrange(aes(x = mean, xmin = lb, xmax = ub)) +
    scale_x_continuous(labels = scales::percent_format(),
                       trans = "log",
                       breaks = c(0, 0.001, 0.01, 0.05, 0.1, 0.2)) +
    labs(x = "Effective IFR", y = "Age group") +
    theme_minimal() +
    theme(axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7),
          legend.position = "none")
  
  if (age) {
    p1 + p2
  } else {
    add_voc_labs(p1, ymax = 0.015)
  }
}


plot_infection_incidence <- function(dat, vocs = FALSE) {
  
  sample <- dat$fit$samples$trajectories$state
  dates <- ZamCovid:::numeric_date_as_date(dat$fit$samples$trajectories$date)
  data <- rep(NA_real_, length(dates))
  
  inf_traj <- c( "infections_inc", "reinfections_inc")
  
  inf <- traj_to_long("infections_inc", 
                      sample["infections_inc", ,], 
                      dates, data)
  
  reinf <- traj_to_long("reinfections_inc", 
                        sample["reinfections_inc", ,], 
                        dates, data)
  
  pop <- sum(dat$fit$parameters$base$population$n)
  ymax <- max((inf$ub / pop)) * 1e3 * 1.05
  
  p <- ggplot(inf, aes(x = date)) +
    geom_line(aes(y = (mean / pop) * 1e3, col = "Mean (95%CrI)"), show.legend = FALSE) +
    geom_ribbon(aes(ymin = (lb / pop) * 1e3,
                    ymax = (ub / pop) * 1e3, fill = "Mean (95%CrI)"),
                alpha = 0.3, show.legend = FALSE) +
    geom_ribbon(aes(ymax = (mean / pop) * 1e3, ymin = (reinf$mean / pop) * 1e3,
                    fill = "1st infections"), alpha = 0.3) +
    geom_ribbon(aes(ymax = (reinf$mean / pop) * 1e3, ymin = 0,
                    fill = "Re-infections"), alpha = 0.5) +
    scale_color_manual(values = c("black", "lightblue4", "blue4"),
                       breaks = c("Mean (95%CrI)", "1st infections",
                                  "Re-infections")) +
    scale_fill_manual(values = c("black", "lightblue4", "blue4"),
                      breaks = c("Mean (95%CrI)", "1st infections",
                                 "Re-infections")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    labs(y = "New infections (x 1,000 population)", x = "") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 12))
  
  if (vocs) {
    add_voc_labs(p, ymax = ymax)
  } else {
    p
  }
}


get_severity <- function(severity, what, by_age = FALSE) {
  
  age_bands <- ZamCovid:::model_age_bins()$start
  dates <- ZamCovid:::numeric_date_as_date(severity$date[-1])
  
  if (by_age) {
    ret <- NULL
    get <- paste0(what, "_age_", age_bands)
  } else {
    ret <- data.frame(date = dates)
    get <- what
  }
  
  for (g in get) {
    sev <- severity[[g]][-1, ]
    
    tmp <- data.frame(t(apply(sev, 1, quantile, c(0.5, 0.025, 0.975),
                              na.rm = TRUE))) %>%
      setNames(., c("mean", "lb", "ub"))
    
    if (by_age) {
     
      tmp <- tmp %>%
       mutate(date = dates,
              age = gsub("ifr_age_", "", g)) %>%
        select(date, age, mean, lb, ub)
      
     ret <- rbind(ret, tmp)
     
    } else {
      ret <- cbind(ret, tmp)
    }
  }
  
  ret
}


get_new_pars <- function(samples, priors) {
  
  i <- which.max(samples$probabilities[, "log_posterior"])
  initial <- samples$pars[i, ]
  
  vcv <- cov(samples$pars)
  dimnames(vcv) <- NULL
  
  for (i in 1:length(priors)) {
    nm <- priors[[i]]$name
    priors[[i]]$initial <- initial[[nm]]
  }
  
  list(vcv = vcv,
       priors = priors)
}


plot_forest <- function(dat, plot_type = "betas") {
  
  labels <- get_par_labels(dat)
  stopifnot(plot_type %in% c("all", "betas", "non_betas"))
  
  samples <- dat$fit$samples
  date <- dat$fit$parameters$base$date
  par_names <- colnames(samples$pars)
  
  beta_date <- ZamCovid:::numeric_date_as_date(samples$info$beta_date)
  beta_names <- par_names[substr(par_names, 1, 4) == "beta"]
  beta_names <- beta_names[order(as.numeric(gsub("beta", "", beta_names)))]
  
  hps <- dat$fit$parameters$prior
  pars_info <- dat$fit$parameters$info
  
  ## Can set the xmax for specific parameters here. For any left as NA, it will
  ## instead just use the maximum from the parameters info
  par_max <- rep(NA, length(par_names))
  names(par_max) <- par_names
  par_max[beta_names] <- 0.25
  
  extract_sample <- function(par_name) {
    lapply(samples, function(x) as.numeric(x$pars[, par_name]))
  }
  
  plot_par <- function(nm, type = NULL, plot = FALSE) {
    
    numeric_par <- as.numeric(samples$pars[, nm])
    numeric_par <- data.frame(
      name = nm,
      mean = mean(numeric_par),
      lb = quantile(numeric_par, 0.025),
      ub = quantile(numeric_par, 0.975)
    )
    
    xlims <- as.vector(pars_info[pars_info$name == nm, c("min", "max")])
      
    if (type == "beta") {
      prior <- hps[hps$name == nm, c("beta_shape1", "beta_shape2")]
      shape1 <- prior$beta_shape1
      shape2 <- prior$beta_shape2
      
      prior <- mapply(qbeta,
                      shape1 = shape1,
                      shape2 = shape2,
                      MoreArgs = list(p = c(0.025, 0.975)),
                      SIMPLIFY = TRUE)
      
      prior <- as.list(prior)
      names(prior) <- c("min", "max")
      
    } else if (type == "gamma") {
      prior <- hps[hps$name == nm, c("gamma_shape", "gamma_scale")]
      shape <- prior$gamma_shape
      scale <- prior$gamma_scale
      
      prior <- mapply(qgamma,
                      shape = shape,
                      scale = scale,
                      MoreArgs = list(p = c(0.025, 0.975)),
                      SIMPLIFY = TRUE)
      
      prior <- as.list(prior)
      names(prior) <- c("min", "max")
      
    } else {
      prior <- xlims
    }
    
    if (plot) {
      ggplot(numeric_par) +
        geom_pointrange(aes(y = NA, x = mean, xmin = lb, xmax = ub)) +
        labs(y = labels[nm], x = "") +
        scale_x_continuous(limits = c(xlims$min, xlims$max)) +
        geom_vline(xintercept = c(prior$min, prior$max), linetype = 3, col = "red") +
        theme_minimal() +
        theme(axis.text.y = element_blank(),
              axis.text.x = element_text(size = 10),
              axis.title = element_text(size = 16, face = "bold"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank())
    } else {
      prior <- data.frame(min = prior$min, max = prior$max)
      cbind(numeric_par, prior)
    }
  }
  
  ## Plot start date
  numeric_start_date <- as.numeric(samples$pars[, "start_date"])
  numeric_start_date <- data.frame(
    name = "start_date",
    mean = ZamCovid:::numeric_date_as_date(mean(numeric_start_date)),
    lb = ZamCovid:::numeric_date_as_date(quantile(numeric_start_date, 0.025)),
    ub = ZamCovid:::numeric_date_as_date(quantile(numeric_start_date, 0.975))
  )

  start_date <- ggplot(numeric_start_date) +
    geom_pointrange(aes(y = NA, x = mean, xmin = lb, xmax = ub)) +
    labs(y = labels$start_date, x = "") +
    scale_x_date(limits = c(as.Date("2020-02-15"), as.Date("2020-04-15")),
                 date_labels = "%d-%b-%Y") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 10, angle = 45, vjust = 0.7),
          axis.title = element_text(size = 16, face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank())
  
  
  if (plot_type == "betas") {
    
    ## Plot betas
    beta_plot <- NULL
    for (i in grep("beta", names(labels), value = TRUE)) {
      beta_plot <- rbind(beta_plot, plot_par(i, "gamma", plot = FALSE))
    }
    
    levels <- paste0("beta", seq(1:nrow(beta_plot)))
    beta_plot$name <- factor(beta_plot$name, levels = levels)
    
    ggplot(beta_plot, aes(group = name)) +
      geom_pointrange(aes(y = NA, x = mean, xmin = lb, xmax = ub)) +
      geom_vline(aes(xintercept = min), linetype = 3, col = "red") +
      geom_vline(aes(xintercept = max), linetype = 3, col = "red") +
      facet_wrap(~name, scales = "free_x") +
      labs(y = "", x = "") +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(size = 10, angle = 45, vjust = 0.7),
            axis.title = element_text(size = 16, face = "bold"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank())
    
    
  } else if (plot_type == "non_betas") {
    ## Plot severity parameters
    sev_plots <- NULL
    for (i in grep("_D", names(labels), value = TRUE)) {
      sev_plots[[i]] <- plot_par(i, "beta", plot = TRUE)
    }
    
    sev <- sev_plots[[1]] / sev_plots[[2]] / sev_plots[[3]] / sev_plots[[4]]
    
    (start_date / sev) +
      plot_layout(heights = c(0.25, 1))
  }
  
}


plot_ci_bar <- function(res, at, width = 1,
                        min = 0.025, max = 0.975, col = "grey20",
                        segments = FALSE, pt_col = NULL, horiz = TRUE, ...) {
  cols  <- c("grey80", col)
  qs <- quantile(res,
                 probs = seq(min, max, by = 0.005),
                 na.rm = TRUE)
  
  palette <- grDevices::colorRampPalette(cols)
  if (segments) {
    segments(y0 = at, x0 = min(res), x1 = max(res), col = cols[2])
    points(y = rep(at, 2), x = range(res), col = cols[2], pch = "-")
  }
  ci_bands(quantiles = cbind(qs, qs),
           y = at + c(-1, 1) * width,
           palette = palette, leg = FALSE, horiz = horiz)
  
  if (is.null(pt_col)) pt_col <- col
  if (horiz) {
    points(y = at, x = mean(res), col = "white", pch = 23, bg = pt_col, ...)
  } else {
    points(x = at, y = mean(res), col = "white", pch = 23, bg = pt_col, ...)
  }
}


ci_bands <- function (quantiles, y, palette = NULL, cols = NULL, leg = TRUE, 
                      leg_y = 0, leg_x = 1, horiz = TRUE, ...) 
{
  yy <- c(y, rev(y))
  yy <- c(yy, yy[1])
  n_bands <- (nrow(quantiles) - 1)/2 + 1
  if (!is.null(palette)) {
    cols <- do.call(what = palette, args = list(n = n_bands))
  }
  for (band in seq_len(n_bands)) {
    x1 <- quantiles[band, ]
    x2 <- quantiles[nrow(quantiles) + 1 - band, ]
    x2 <- rev(x2)
    x2 <- c(x2, x1[1])
    if (horiz) {
      polygon(y = yy, x = c(x1, x2), col = cols[band], 
              border = NA)
    }
    else {
      polygon(x = yy, y = c(x1, x2), col = cols[band], 
              border = NA)
    }
  }
  if (leg) {
    leg_cols <- which(row.names(quantiles) %in% leg)
    leg <- c(row.names(quantiles)[1], seq(5, 50, 5), "%")
    leg[seq(2, 10, 2)] <- ""
    legend(y = leg_y, x = leg_x, pch = 15, col = cols[leg_cols], 
           legend = leg, border = NA, bty = "n", ...)
  }
}


get_par_labels <- function(dat) {
  
  par_names <- unique(dat$fit$parameters$proposal$name)
  beta_date <- ZamCovid:::numeric_date_as_date(dat$fit$samples$info$beta_date)
  
  n_betas <- length(beta_date)
  
  labels <- lapply(seq_len(n_betas), 
                   function(x) {
                     date_x <- as.character(beta_date[x])
                     bquote(beta[.(x)] ~ (.(date_x)))
                   })
  names(labels) <- paste0("beta", seq_len(n_betas))
  
  plus <- "+"
  labels <- c(labels,
              alpha_D = expression(alpha["D"]),
              mu_D_1 = expression(mu["D,1"]),
              mu_D_2 = expression(mu["D,2"]),
              p_G_D = expression(p[paste(G[D])]^max),
              start_date = expression(t["seed"])
  )
  labels
}


get_cum_inf_age <- function(states, t) {
  
  stopifnot(t %in% c("all", "over15", "15_19", "20_29", "30_39",
                     "40_49", "50_plus"))
  
  inf_nm <- grep("infections_inc_age_", rownames(states), value = TRUE)
  reinf_nm <- grep("reinfections_inc_age_", rownames(states), value = TRUE)
  inf_nm <- setdiff(inf_nm, reinf_nm)
  
  if (t == "over15") {
    inf_nm <- inf_nm[4:length(inf_nm)]
    reinf_nm <- reinf_nm[4:length(reinf_nm)]
  } else if (t == "15_19") {
    inf_nm <- inf_nm[4]
    reinf_nm <- reinf_nm[4]
  } else if (t == "20_29") {
    inf_nm <- inf_nm[5:6]
    reinf_nm <- reinf_nm[5:6]
  } else if (t == "30_39") {
    inf_nm <- inf_nm[7:8]
    reinf_nm <- reinf_nm[7:8]
  } else if (t == "40_49") {
    inf_nm <- inf_nm[9:10]
    reinf_nm <- reinf_nm[9:10]
  } else if (t == "50_plus") {
    inf_nm <- inf_nm[11:length(inf_nm)]
    reinf_nm <- reinf_nm[11:length(reinf_nm)]
  }
  
  inf <- states[inf_nm, , -1]
  reinf <- states[reinf_nm, , -1]
  
  if (length(dim(inf)) == 3) {
    inf <- apply(inf, c(2:3), sum)
    reinf <- apply(reinf, c(2:3), sum)
  }
  
  out <- inf - reinf
  t(apply(out, 1, cumsum))
}
