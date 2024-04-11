source("util.R")

version_check("ZamCovid", "0.1.3")
n_threads <- as.integer(Sys.getenv("CONTEXT_CORES",
                                   Sys.getenv("MC_CORES",
                                              getOption("mc.cores", 1))))

dat <- readRDS("inputs/fits.rds")$fit
dat$steps_per_day <- 4
dat$dt <- 1 / dat$steps_per_day
# This is a sanity check - 0 == 100% vaccine efficacy; 1 == no VE
# for (i in seq_along(dat$parameters$base$rel_severity)) {
#   tmp <- dat$parameters$base$rel_severity[[i]]
#   tmp[] <- 1
#   dat$parameters$base$rel_severity[[i]] <- tmp
# }

## Prepare control parameters for general simulation
date_end <- as.Date(fit_date)
control <- simulate_control(dat, date_end)
date_start <- as.Date(
  ZamCovid:::numeric_date_as_date(sapply(control$scenarios, function(x) x$start_date)))


## Prepare fits for simulation
fits_prepared <- fits_simulate_prepare(dat, date_start, n_par)
fits_real <- c(dat$samples$trajectories[c("date", "state")],
               dat$rt[c("Rt_general", "eff_Rt_general")],
               dat$data,
               dat$severity["ifr"])
dat <- NULL


## Create counterfactual vaccination schedules
# real_vaccine_schedule <- fits_prepared$base$vaccine_schedule
counterfactual_vaccine_schedule <-
  create_vaccination_counterfactual(control, fits_prepared$base)

## Now prepare future beta_step for simulation based on Rt scenario(s)
control <- compute_beta_step(fits_prepared, control)

## Update model parameters for vaccine simulation
pars <- update_simulate_pars(fits_prepared, control, counterfactual_vaccine_schedule)
counterfactual_vaccine_schedule <- NULL

## Running the simulation!
res <- simulate_local(pars, control, fits_prepared, n_threads, seed = 4321)
names(res) <- names(control$scenarios)


## Save outputs 
saveRDS(list(simulation = res,
             date_start = date_start,
             fits_prepared = fits_prepared,
             fits_real = fits_real,
             parameters = pars),
        "simulation_results.rds")


## Plot simulation results
dir.create("plots", FALSE, TRUE)


## Main manuscript plots
png("plots/figure_2.png", units = "in", width = 10, height = 8, res = 300)
(plot_scenario_cumulative(fits_real, res, sample = fits_prepared$index,
                         outcome = "deaths_covid") /
  plot_scenario_cumulative(fits_real, res, sample = fits_prepared$index,
                           outcome = "deaths_covid", which = "healthcare")) +
  plot_annotation(tag_levels = "A")
dev.off()


scen_nms <- names(res)[c(2, 5, 4, 3, 7, 10, 9, 8, 12, 15, 14, 13)]
png("plots/figure_3.png", units = "in", width = 10, height = 8, res = 300)
plot_effective_immunity(fits_real, res, pars, scen_plot = scen_nms,
                        add_ever_infected = TRUE,
                        sample = fits_prepared$index)
dev.off()


## SI plots
png("plots/si_inf_averted.png", units = "in", width = 10, height = 8, res = 300)
plot_scenario_cumulative(fits_real, res, sample = fits_prepared$index,
                         outcome = "inf_cum_all") /
  plot_scenario_cumulative(fits_real, res, sample = fits_prepared$index,
                           outcome = "inf_cum_all", which = "healthcare")
dev.off()

png("plots/si_vacc_infections.png", units = "in", width = 10, height = 8, res = 300)
plot_scenario_timeseries(fits_real, res, outcome = "infections_inc")
dev.off()

png("plots/si_vacc_deaths.png", units = "in", width = 10, height = 8, res = 300)
plot_scenario_timeseries(fits_real, res, outcome = "deaths_covid")
dev.off()

png("plots/si_vacc_prop_vaccinated.png", units = "in", width = 10, height = 8, res = 300)
plot_vaccinated(res[1:15], pars)
dev.off()

png("plots/si_vacc_prop_immune.png", units = "in", width = 10, height = 8, res = 300)
plot_effective_immunity(fits_real, res, pars)
dev.off()

png("plots/si_vacc_Rt.png", units = "in", width = 10, height = 8, res = 300)
plot_rt(fits_real, res)
dev.off()

png("plots/si_vacc_ifr.png", units = "in", width = 10, height = 8, res = 300)
plot_severity(fits_real, res)
dev.off()

png("plots/si_health_infections.png", units = "in", width = 10, height = 8, res = 300)
plot_scenario_timeseries(fits_real, res, outcome = "infections_inc",
                             which = "healthcare")
dev.off()

png("plots/si_health_deaths.png", units = "in", width = 10, height = 8, res = 300)
plot_scenario_timeseries(fits_real, res, outcome = "deaths_covid",
                             which = "healthcare")
dev.off()

png("plots/si_health_prop_immune.png", units = "in", width = 10, height = 8, res = 300)
plot_effective_immunity(fits_real, res, pars, which = "healthcare")
dev.off()

png("plots/si_health_Rt.png", units = "in", width = 10, height = 8, res = 300)
plot_rt(fits_real, res, which = "healthcare")
dev.off()

png("plots/si_health_ifr.png", units = "in", width = 10, height = 8, res = 300)
plot_severity(fits_real, res, which = "healthcare")
dev.off()


## Render rmd
rmarkdown::render("paper_numbers.Rmd")
