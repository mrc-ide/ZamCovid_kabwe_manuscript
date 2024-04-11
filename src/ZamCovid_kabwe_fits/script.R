source("util.R")
version_check("ZamCovid", "0.1.3")


## 1. Prepare elements of particle filter
pars <- fit_pars_load("inputs", "kabwe", assumptions,
                      short_run, deterministic)

restart_date <- readRDS("inputs/base.rds")[[1]]$restart_date
control <- set_control(short_run, deterministic,
                       restart_date = restart_date)

data_full <- read_csv("data/data_timeseries.csv")
data_fit <- parse_data(data_full, fit_sero = TRUE, sero_by_age = TRUE,
                       fit_deaths = TRUE)

## 2. Build particle filter and run pMCMC
filter <- ZamCovid_particle_filter(data_fit, pars$mcmc,
                                   control$particle_filter,
                                   deterministic = deterministic)

samples <- fit_run(pars, filter, control$pmcmc)


## 3. Post-processing of model fits ----
dat <- ZamCovid_fit_process(samples, pars, data_full, data_fit)

dir.create("outputs", FALSE, TRUE)
saveRDS(dat, "outputs/fit.rds")
write_csv(dat$fit$parameters$info, "outputs/info.csv")
write_csv(dat$fit$parameters$proposal, "outputs/proposal.csv")

dir.create("plots", FALSE, TRUE)
message("Creating plots")

png("plots/si_pmcmc_traceplots.png", units = "in", width = 16, height = 12, res = 300)
plot_fit_traces(dat$fit$samples)
dev.off()

png("plots/si_fits_serology.png", units = "in", width = 8, height = 8, res = 300)
plot_serology(dat)
dev.off()

png("plots/si_fits_deaths.png", units = "in", width = 10, height = 10, res = 300)
plot_deaths(dat) | plot_deaths_disag(dat)
dev.off()

png("plots/si_forest_betas.png", units = "in", width = 14, height = 12, res = 300)
plot_forest(dat, plot_type = "betas")
dev.off()

png("plots/si_forest_non_betas.png", units = "in", width = 6, height = 10, res = 300)
plot_forest(dat, plot_type = "non_betas")
dev.off()


# Figure 1 for the manuscript
col1 <- (
  (plot_serology(dat, which = "over15",
                 labels = "≥15-yo.", immunity = TRUE) +
     theme(legend.position = "right",
           axis.text.x = element_blank())) /
    (plot_rt(dat, vocs = TRUE, "eff_rt") + theme(axis.text.x = element_blank())) /
    plot_infection_incidence(dat, vocs = TRUE)) +
  plot_layout(heights = c(1, 1, 1))

col2 <- ((plot_deaths(dat, week_only = TRUE) +
  theme(axis.text.x = element_blank(),
        legend.position = "none")) /
  (plot_severity(dat, age = FALSE) +
     theme(axis.text.x = element_blank())) / 
  plot_deaths_disag(dat, plot_age = FALSE)) +
  plot_layout(heights = c(1, 1, 1))

png("plots/figure_1.png", units = "in", width = 10, height = 8, res = 300)
(col1 | col2) +
  plot_annotation(tag_levels = list(c("A", "C", "E", "B", "D", "F")))
dev.off()


rmarkdown::render("fit_results.Rmd")
