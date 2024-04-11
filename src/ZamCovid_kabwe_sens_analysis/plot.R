ggtheme <- theme_minimal() +
  theme(legend.position = "none",
        axis.line = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 0.7))
theme_set(ggtheme)


plot_trajectory_scenarios <- function(scenarios, data, which = "deaths_all_inc") {
  
  # Unlist objects of interest into tibbles for plotting
  traj <- unlist_scenarios(scenarios, "trajectories")

  # Some graphic parameters
  cols <- c("grey30", "goldenrod", "darkolivegreen3")
  col_names <- c('Central fit', 'Scenario', 'Data')
  
  # Plot trajectories
  if (which == "deaths_all_inc") {
    ymax <- 100
    dta <- data %>%
      select(date, data_mean = deaths_all) %>%
      mutate(data_lb = NA_real_, data_ub = NA_real_)
  } else if (which == "sero_pos_over15") {
    ymax <- 1
    dta <- data %>%
      select(date, data_mean = sero_mean, data_lb = sero_lb, data_ub = sero_ub)
  }
  
  tmp <- traj %>%
    filter(state == which)
  
  ctr <- tmp %>%
    filter(scenario == "central") %>%
    mutate(source = "Central fit") %>%
    left_join(., dta)
  
  tmp <- tmp %>%
    filter(scenario != "central") %>%
    mutate(source = "Scenario") %>%
    left_join(., dta) %>%
    rbind(., ctr)
  
  if (which == "deaths_all_inc") {
    tmp <- tmp %>%
      mutate(date = lubridate::floor_date(date, "week", week_start = 3)) %>%
      pivot_longer(!c(date, scenario, state, source)) %>%
      group_by(date, scenario, state, source, name) %>%
      summarise(value = sum(value, na.rm = TRUE)) %>%
      ungroup() %>%
      pivot_wider() %>%
      filter(date < as.Date("2021-09-29"))
    
    p <- ggplot(tmp, aes(date)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
      labs(y = "Weekly all-cause deaths", x = "")
    
    
  } else if (which == "sero_pos_over15") {
    p <- ggplot(tmp, aes(date)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, ymax),
                         labels = scales::percent_format()) +
      labs(y = "Seropositive (15+)", x = "")
  }

  p +
    geom_pointrange(aes(y = data_mean, ymin = data_lb,
                        ymax = data_ub, col = "Data"), alpha = 0.8, size = 0.15) +
    geom_line(aes(y = mean, col = source)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = source),
                alpha = 0.2, show.legend = FALSE) +
    facet_wrap(~scenario) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    scale_color_manual(values = cols, breaks = col_names) +
    scale_fill_manual(values = cols, breaks = col_names) +
    theme(legend.position = "top",
          legend.title = element_blank())
}


plot_ifr_scenarios <- function(scenarios) {
  
  ifr <- unlist_scenarios(scenarios, "severity") %>%
    mutate(source = case_when(scenario == "central" ~ "Central fit",
                              TRUE ~ "Scenario")) %>%
    filter(date > as.Date("2020-03-11"))
  
  # Some graphic parameters
  cols <- c("grey30", "goldenrod")
  col_names <- c('Central fit', 'Scenario')
  
  ggplot(ifr, aes(date)) +
    geom_hline(yintercept = 0.0012, col = "red", linetype = 3) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.015),
                       labels = scales::percent_format()) +
    labs(y = "Infection fatality ratio", x = "") +
    geom_line(aes(y = mean, col = source)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = source),
                alpha = 0.2, show.legend = FALSE) +
    facet_wrap(~scenario) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    scale_color_manual(values = cols, breaks = col_names) +
    scale_fill_manual(values = cols, breaks = col_names) +
    theme(legend.position = "top",
          legend.title = element_blank())
}


plot_rt_scenarios <- function(scenarios) {
  
  rt <- unlist_scenarios(scenarios, "rt") %>%
    mutate(source = case_when(scenario == "central" ~ "Central fit",
                              TRUE ~ "Scenario")) %>%
    filter(date > as.Date("2020-03-11"))
  
  # Some graphic parameters
  ggplot(rt, aes(date)) +
    geom_hline(yintercept = 1, col = "red", linetype = 3) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    labs(y = "Reproduction number", x = "") +
    geom_line(aes(y = mean, col = state)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = state),
                alpha = 0.2, show.legend = FALSE) +
    facet_wrap(~scenario) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    theme(legend.position = "top",
          legend.title = element_blank())
}
