# ZamCovid_kabwe_manuscript

Public repository to the manuscript

> Pandemic burden in low-income settings and impact of limited
and delayed interventions: a granular modelling analysis of COVID-19 in Kabwe, Zambia

This is an [`orderly`](https://github.com/vimc/orderly) project. The directories are:

* `src`: create new reports here
* `global`: R scripts with some basic functions used across reports
* `draft`: versioned results of running your report
* `archive`: committed reports from draft

## Running the analysis

A sequence of tasks needs to be run to infer (fit) the model parameters to reproduce all results and plots in the manuscript. This is sketched out as a workflow in the [`run.R`](run.R) script, it is also described here for documentation purposes.

For the publication manuscript, the fitting tasks were run over several days on an HPC server. However, shorter MCMC (deterministic) and pMCMC (stochastic) chains can be run on a local computer as follows:

1. Run the `ZamCovid_kabwe_data` task to prepare the data for fitting.
2. Run the `ZamCovid_kabwe_parameters` task.
3. Run the `ZamCovid_kabwe_fits` task for each of the central or sensitivity analysis scenarios.
4. Run the `ZamCovid_kabwe_sens_analysis` task.
5. Run the `ZamCovid_kabwe_simulation` task.

Note the `ZamCovid_kabwe_sens_analysis` task will require a matching set of parameters (2) and fitting (3) tasks with all of central and other assumptions used for sensitivity analysis, as specified in `run.R`. Also note that tasks 2 to 4 can be run either with the deterministic equivalent or the full stochastic model. The `ZamCovid_kabwe_simulation` task, however, depends on the `ZamCovid_kabwe_fits` task with the stochastic model and central parameters.

## Requirements

The core requirement is our [ZamCovid](https://mrc-ide.github.io/ZamCovid/) package and its dependencies. Please ensure you install the versions we used for preparation of the manuscript:

```r
remotes::install_github(c(
  "mrc-ide/dust@v0.15.1",
  "mrc-ide/mcstate@v0.9.20",
  "mrc-ide/ZamCovid@v0.1.3"))
```

You will also need a recent [orderly](https://www.vaccineimpact.org/orderly/) which can be installed with

```r
drat:::add("vimc")
install.packages("orderly")
```

## License

MIT © Imperial College of Science, Technology and Medicine
