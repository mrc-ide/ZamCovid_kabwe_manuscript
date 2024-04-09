setwd(orderly::orderly_config()$root)
source("global/util.R")

add_parameter <- function(name, initial, min, max, proposal,
                          assumptions = "central", model = "deterministic",
                          multidistrict = FALSE, integer = FALSE, include = TRUE) {
  
  if (multidistrict) {
    dir <- "src/ZamCovid_parameters_multidistrict/pars"
  } else {
    dir <- "src/ZamCovid_kabwe_parameters/pars"
  }
  
  info_filename <-
    paste(dir, assumptions, model,
          "info.csv", sep = "/")
  proposal_filename <- 
    paste(dir, assumptions, model,
          "proposal.csv", sep = "/")
  
  parameters_info <- read_csv(info_filename) 
  
  districts <- unique(parameters_info$district)
  districts <- districts[!is.na(districts)]
  
  new_par <- data.frame(
    district = districts,
    name = name,
    initial = initial,
    min = min,
    max = max,
    integer = integer,
    include = include
  )
  
  parameters_info <- rbind(parameters_info, new_par)
  parameters_info <- dplyr::arrange(parameters_info,
                                    district, name)

  write.csv(parameters_info, info_filename, row.names = FALSE)
  
  parameters_proposal <- read_csv(proposal_filename)
  
  parameters_proposal[[name]] <- 0
  k <- parameters_proposal$name == parameters_proposal$name[1]
  
  new_prop <- parameters_proposal[k, ]
  new_prop$name <- name
  new_prop[, -c(1, 2)] <- 0
  
  if (!is.numeric(proposal)) {
    avg_prop <- NULL
    
    for (i in regions) {
      
      avg_prop[i] <- parameters_proposal[
        parameters_proposal$district == i &
          parameters_proposal$name == proposal, proposal]
    }
    proposal = abs(mean(avg_prop))
    stopifnot(is.numeric(proposal))
  }
  new_prop[[name]] <- proposal
  
  parameters_proposal <- rbind(parameters_proposal, new_prop)
  parameters_proposal <- dplyr::arrange(parameters_proposal, district, name)
  col_order <- c(1, 2, 2 + order(names(parameters_proposal)[-c(1,2)]))
  parameters_proposal <- parameters_proposal[, col_order]
  
  write.csv(parameters_proposal, proposal_filename, row.names = FALSE)
}


remove_parameter <- function(name, assumptions = "central", model = "deterministic",
                             multidistrict = FALSE) {
  
  if (multidistrict) {
    dir <- "src/ZamCovid_parameters_multidistrict/pars"
  } else {
    dir <- "src/ZamCovid_kabwe_parameters/pars"
  }
  
  info_filename <-
    paste(dir, assumptions, model,
          "info.csv", sep = "/")
  proposal_filename <- 
    paste(dir,assumptions, model,
          "proposal.csv", sep = "/")
  
  parameters_info <- read_csv(info_filename)
  parameters_info <- parameters_info[parameters_info$name != name, ]
  write.csv(parameters_info, info_filename, row.names = FALSE)
  
  
  parameters_proposal <- read_csv(proposal_filename)
  parameters_proposal[[name]] <- NULL
  parameters_proposal <- parameters_proposal[parameters_proposal$name != name, ]
  
  write.csv(parameters_proposal, proposal_filename, row.names = FALSE)
  
}


rename_parameter <- function(old_name, new_name, multidistrict = FALSE,
                             assumptions = "central", model = "deterministic") {
  
  if (multidistrict) {
    dir <- "src/ZamCovid_parameters_multidistrict/pars"
  } else {
    dir <- "src/ZamCovid_kabwe_parameters/pars"
  }
  
  info_filename <-
    paste(dir, assumptions, model,
          "info.csv", sep = "/")
  proposal_filename <- 
    paste(dir, assumptions, model,
          "proposal.csv", sep = "/")
  
  
  parameters_info <- read_csv(info_filename)
  parameters_info$name[parameters_info$name == old_name] <- new_name
  parameters_info <- dplyr::arrange(parameters_info,
                                    region, name)

  write.csv(parameters_info, info_filename, row.names = FALSE)
  
  parameters_proposal <- read_csv(proposal_filename)
  parameters_proposal$name[parameters_proposal$name == old_name] <- new_name
  parameters_proposal <- dplyr::arrange(parameters_proposal, region, name)
  names(parameters_proposal)[names(parameters_proposal) == old_name] <- new_name
  col_order <- c(1, 2, 2 + order(names(parameters_proposal)[-c(1,2)]))
  parameters_proposal <- parameters_proposal[, col_order]

  write.csv(parameters_proposal, proposal_filename, row.names = FALSE)
}


add_beta <- function(beta_name, beta_initial, min = 0, max = 1,
                     proposal = NULL, factor = NULL, multidistrict = FALSE,
                     assumptions = "central", model = "deterministic") {
  
  if (multidistrict) {
    dir <- "src/ZamCovid_parameters_multidistrict/pars"
  } else {
    dir <- "src/ZamCovid_kabwe_parameters/pars"
  }
  
  info_filename <-
    paste(dir, assumptions, model,
          "info.csv", sep = "/")
  proposal_filename <- 
    paste(dir, assumptions, model,
          "proposal.csv", sep = "/")
  
  parameters_info <- read_csv(info_filename)
  
  regions <- unique(parameters_info$region)
  regions <- regions[!is.na(regions)]
  
  if (!is.null(factor)) {
    f <- factor
  } else {
    f <- 1
  }
  
  new_par <- data.frame(
    region = regions,
    name = beta_name,
    initial = parameters_info$initial[parameters_info$name == beta_initial] * f,
    min = min,
    max = max,
    integer = FALSE,
    include = TRUE
  )
  
  parameters_info <- rbind(parameters_info, new_par)
  parameters_info <- dplyr::arrange(parameters_info,
                                    region, name)
  
  write.csv(parameters_info, info_filename, row.names = FALSE)
  
  
  parameters_proposal <- read_csv(proposal_filename)
  
  if (is.null(proposal)) {
    avg_prop <- NULL
    regions <- unique(parameters_proposal$region)
    for (i in regions[!is.na(regions)]) {
      
      avg_prop[i] <- parameters_proposal[
        parameters_proposal$region == i &
          parameters_proposal$name == beta_initial, beta_initial]
    }
    proposal = abs(mean(avg_prop))
  }
  
  parameters_proposal[[beta_name]] <- 0
  k <- parameters_proposal$name == "beta1"
  new_prop <- parameters_proposal[k, ]
  new_prop$name <- beta_name
  new_prop[, -c(1, 2)] <- 0
  new_prop[[beta_name]] <- proposal
  
  parameters_proposal <- rbind(parameters_proposal, new_prop)
  parameters_proposal <- dplyr::arrange(parameters_proposal, region, name)
  col_order <- c(1, 2, 2 + order(names(parameters_proposal)[-c(1,2)]))
  parameters_proposal <- parameters_proposal[, col_order]
  
  write.csv(parameters_proposal, proposal_filename, row.names = FALSE)
}



nudge_info <- function(region, pars, factor = NULL, value = NULL,
                       multidistrict = FALSE,
                       assumptions = "central", deterministic = TRUE) {
  
  if (multidistrict) {
    dir <- "src/ZamCovid_parameters_multidistrict/pars"
  } else {
    dir <- "src/ZamCovid_kabwe_parameters/pars"
  }
  
  if (deterministic) {
    deterministic <- "deterministic"
  } else {
    deterministic <- "stochastic"
  }
  
  info_filename <-
    paste(dir, assumptions, deterministic,
          "info.csv", sep = "/")
  
  info <- read_csv(info_filename)
  
  for (p in pars) {
    
    init <- info[info$region == region & info$name == p, "initial"]
    new <- NULL
    if (!is.null(factor)) {
      new <- init * factor
      info[info$region == region & info$name == p, "initial"] <- new
    } else if (!is.null(value)) {
      new <- value
      info[info$region == region & info$name == p, "initial"] <- new
    } else {
      stop("Need to specify either a scalar factor or value")
    }
    print(paste(p, "for", region, "changed from", init, "to", new))
  }
  
  write.csv(info, info_filename, row.names = FALSE)
  
}


add_district <- function(name, source, assumptions = "central",
                         model = "deterministic") {
  
  dir <- "src/ZamCovid_parameters_multidistrict/pars"
  
  info_filename <-
    paste(dir, assumptions, model,
          "info.csv", sep = "/")
  proposal_filename <- 
    paste(dir, assumptions, model,
          "proposal.csv", sep = "/")
  
  parameters_info <- read_csv(info_filename) 
  
  new_info <- parameters_info[parameters_info$district == source, ]
  new_info$district <- name
  
  parameters_info <- rbind(parameters_info, new_info)
  write.csv(parameters_info, info_filename, row.names = FALSE)
  
  
  parameters_proposal <- read_csv(proposal_filename)
  
  new_prop <- parameters_proposal[parameters_proposal$district == source, ]
  new_prop$district <- name
  
  parameters_proposal <- rbind(parameters_proposal, new_prop)
  write.csv(parameters_proposal, proposal_filename, row.names = FALSE)
}
