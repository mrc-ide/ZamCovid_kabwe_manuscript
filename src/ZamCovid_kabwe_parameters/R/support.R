load_mcmc_parameters <- function(assumptions, deterministic) {
  
  message(sprintf("Assumptions are with '%s' parameters", assumptions))
  
  deterministic <- ifelse(deterministic, "deterministic", "stochastic")
  path_pars <- file.path("pars", assumptions, deterministic)
  
  info <- read_csv(file.path(path_pars, "info.csv"))
  proposal <- read_csv(file.path(path_pars, "proposal.csv"))
  
  prior <- create_priors(info)
  proposal <- update_proposal(info, proposal)
  
  ret <- list(info = info,
              proposal = proposal,
              prior = prior)
  
  ret
}


update_proposal <- function(info, proposal) {
  msg <- setdiff(info$name, proposal$name)
  if (length(msg) > 0) {
    message(sprintf("Adding %d parameters to proposal matrix", length(msg)))
    extra <- expand.grid(region = unique(proposal$region), name = msg)
    extra[names(proposal)[-(1:2)]] <- 0
    proposal <- rbind(proposal, extra)
    proposal[msg] <- 0
    proposal <- proposal[order(proposal$region, proposal$name), ]
    col_order <- c(1, 2, 2 + order(names(proposal)[-c(1, 2)]))
    proposal <- proposal[, col_order, drop = FALSE]
  }
  proposal
}


infer_baseline_deaths <- function(historic, date, alpha = 0.95,
                                  inflate = 1 / 0.516) {
  
  expected_model <- historic_deaths %>%
    mutate(deaths = floor(deaths * inflate)) %>%
    nest_by(site) %>%
    mutate(model = list(lm(deaths ~ year + month, data = data))) %>%
    mutate(sigma = broom::glance(model)$sigma,
           nobs = broom::glance(model)$nobs,
           df = model$df,
           t_df = abs(qt((1 - alpha) / 2, df)))
  
  expected_deaths <- expected_model %>%
    summarise(
      df = model$df,
      broom::augment(model,
                     newdata = tibble(month = rep(c(1:12), 2) %>% as.factor(),
                                      year = rep(c(2020, 2021), each = 12)),
                     interval = "prediction",
                     se_fit = TRUE)) %>%
    mutate(month_num = as.numeric(month)) %>%
    rename(expected = .fitted,
           lb = .lower,
           ub = .upper) %>%
    mutate(t_df = abs(qt((1 - alpha) / 2, df)),
           se_fit_pred = (ub - expected) / t_df)
  
  ## Output linear model diagnostics for scrutiny
  mod <- expected_model$model[[1]]
  dat <- expected_model$data[[1]]
  
  names(mod$coefficients)[3:13] <- month.abb[-1]
  
  p1 <- jtools::plot_summs(mod, plot.distributions = TRUE)
  p2 <- ggplot(dat, aes(x = year, y = deaths)) +
    geom_point() +
    scale_x_continuous(breaks = c(2017:2019)) +
    scale_y_continuous(expand = c(0, 0), limits = c(100, 250)) +
    stat_smooth(method = "lm") +
    theme_minimal() +
    theme(axis.line = element_line())
  
  # Extract only death rate from model for replacement function
  expected_deaths <- baseline_death_rate(expected_deaths, date)
  
  list(
    expected_deaths = expected_deaths,
    model = mod,
    plot = p1 + p2
  )
}


baseline_death_rate <- function(expected_deaths, d) {
  
  out <- expected_deaths %>%
    ungroup() %>%
    mutate(date = as.Date(paste0(year, "-", month, "-01")),
           days = lubridate::days_in_month(date),
           rate = round(expected / days)) %>%
    select(date, rate)
  
  if (out$date[nrow(out)] > d) {
    out[out$date <= d, ]
  } else {
    out
  }
}


prior_rename <- function (x, from, to, name = deparse(substitute(x))) {
  verify_names(x, required = from, allow_extra = TRUE, name = name)
  i <- match(from, names(x))
  names(x)[i] <- to
  x
}


verify_names <- function (x, required = NULL, optional = NULL, allow_extra = FALSE, 
                          name = deparse(substitute(x))) {
  nms <- names(x)
  if (anyDuplicated(nms)) {
    dups <- unique(nms[duplicated(nms)])
    stop(sprintf("Duplicate element names in '%s': %s", name, 
                 paste(squote(dups), collapse = ", ")))
  }
  if (!allow_extra) {
    extra <- setdiff(nms, c(required, optional))
    if (length(extra) > 0) {
      stop(sprintf("Extra elements in '%s': %s", name, 
                   paste(squote(extra), collapse = ", ")))
    }
  }
  msg <- setdiff(required, nms)
  if (length(msg) > 0) {
    stop(sprintf("Elements missing from '%s': %s", name, 
                 paste(squote(msg), collapse = ", ")))
  }
  invisible(x)
}


squote <- function (x) {
  sprintf("'%s'", x)
}
