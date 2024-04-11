source("util.R")
version_check("ZamCovid", "0.1.0")

# Read scenario inputs
scenarios <- load_scenarios("fits/")
data <- load_data()
pars <- readRDS("fits/central_fits.rds")$fit$parameters

dir.create("plots", TRUE, FALSE)
# Outputs comparative plots
png("plots/sa_deaths.png", units = "in", width = 8, height = 8, res = 300)
plot_trajectory_scenarios(scenarios, data, "deaths_all_inc")
dev.off()

png("plots/sa_seropositive.png", units = "in", width = 8, height = 8, res = 300)
plot_trajectory_scenarios(scenarios, data, "sero_pos_over15")
dev.off()

png("plots/sa_ifr.png", units = "in", width = 8, height = 8, res = 300)
plot_ifr_scenarios(scenarios)
dev.off()

png("plots/sa_rt.png", units = "in", width = 8, height = 8, res = 300)
plot_rt_scenarios(scenarios)
dev.off()
