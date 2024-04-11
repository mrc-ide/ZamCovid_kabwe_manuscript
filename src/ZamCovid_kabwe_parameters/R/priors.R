create_priors <- function(pars_info) {
  
  ### Set a beta prior for p_G_D
  # We will zero the probabilities of going into hospital for now as there are
  # no reliable sources of hospitalisation data. Fitted p_G_D will reflect
  # overall probability of death by mechanistically assuming all deaths
  # happen outside of hospital. This simplifies things, as we do not have
  # to make any wild assumptions around p_H or p_H_D in relation to p_G_D!
  regional_ps <- data.frame(
    param = "p_G_D",
    mean = 0.0916337063, # Max p_G_D approximated from Brazeau et al. IFRs
    lower = 0.05 # mean - lower = lower range of prior
  )
  
  p_hps <- signif(mapply(FUN = fit_beta,
                         mean = regional_ps$mean,
                         lower = regional_ps$mean - regional_ps$lower,
                         ci = 0.95))
  p_hps <- as.data.frame(t(p_hps))
  colnames(p_hps) <- c("shape1", "shape2")
  
  regional_ps <- cbind(regional_ps, p_hps)
  
  ### set beta_priors
  beta_hps <- data.frame(
    scale = rep(NA, 2),
    shape = rep(NA, 2)
  )
  row.names(beta_hps) <- c("beta1", "beta2")
  
  ## beta value that would give R0 = 1
  ## assumes Wildtype R0 was ~ 2.58 in Zambia
  R0_fac <- 0.0367
  
  ## beta1 aim for 95% CI of [2.5, 3.5]
  beta_hps["beta1", ] <- fit_gamma(mean_D = 2.58,
                                   lower_D = 2.025,
                                   upper_D = 3.125,
                                   ci = 0.95)
  beta_hps["beta1", "scale"] <- beta_hps["beta1", "scale"] * R0_fac
  ## beta2 aim for 95% CI of [1.125, 3.125]
  beta_hps["beta2", ] <- fit_gamma(mean_D = 1.756,
                                   lower_D = 1.125,
                                   upper_D = 3.125,
                                   ci = 0.95)
  beta_hps["beta2", "scale"] <- beta_hps["beta2", "scale"] * R0_fac

  beta_hps <- beta_hps[c(rep("beta1", 10), rep("beta2", 32)), ]
  rownames(beta_hps) <- paste0("beta", seq_len(nrow(beta_hps)))
  beta_names <- rownames(beta_hps)
  
  pars <- c(beta_names, unique(regional_ps$param))
  
  ## Make data frame of all parameters to fit
  hps <- matrix(NA, nrow = length(pars), ncol = 7,
                dimnames = list(pars, c("par", "region", "scale", "shape",
                                        "shape1", "shape2", "correlation")))
  hps <- as.data.frame(hps)
  hps$par <- pars
  hps$region <- "kabwe"
  hps[beta_names, colnames(beta_hps)] <- beta_hps
  hps[unique(regional_ps$param), c("shape1", "shape2")] <- as.matrix(p_hps)
  
  ret <- priors_wide_to_long(hps)
  
  par <- c("alpha_D", "mu_D_1", "mu_D_2")
  
  extra_uniform <-
    expand.grid(region = "kabwe",
                type = "null",
                name = par,
                gamma_scale = NA_real_,
                gamma_shape = NA_real_,
                beta_shape1 = NA_real_,
                beta_shape2 = NA_real_,
                stringsAsFactors = FALSE)
  ret <- rbind(ret, extra_uniform)
  
  nms_expected <- unique(pars_info$name)
  nms_found <- unique(ret$name)
  msg <- setdiff(nms_expected, nms_found)
  if (length(msg) > 0) {
    stop(sprintf("Missing parameters, update priors (missing %s)",
                 paste(msg, collapse = ", ")))
  }
  extra <- setdiff(nms_found, nms_expected)
  if (length(extra)) {
    message(sprintf("Dropping %d unused priors: %s",
                    length(extra), paste(extra, collapse = ", ")))
    ret <- ret[ret$name %in% nms_expected, ]
  }
  rownames(ret) <- NULL
  
  invisible(ret)
}


## specify gamma in terms of mean and variance rather than shape and scale

# convert mean and variance of gamma to shape and scale
mv_to_ss <- function(mean, var) {
  scale <- var / mean
  shape <- mean / scale
  list(shape = shape, scale = scale)
}

# gamma dist functions
qgammamv <- function(p, mean, var) {
  X <- mv_to_ss(mean, var)
  qgamma(p = p, shape = X$shape, scale = X$scale)
}

dgammav <- function(x, mean, var) {
  X <- mv_to_ss(mean, var)
  dgamma(x = x, shape = X$shape, scale = X$scale)
}


## fitting function by least-squares based on mean and CIs
fit_gamma <- function(mean_D, lower_D, upper_D, ci = 0.99) {
  
  alpha <- (1 - ci)/2
  p <- c(alpha , 1 - alpha)
  
  f <- function(v) {
    x <- qgammamv(p = p, mean = mean_D, var = v)
    sum((x[1] - c(lower_D))^2)
  }
  
  var_D <- optimise(f = f, interval = c(0,10), maximum = FALSE)$minimum
  X <- mv_to_ss(mean_D, var_D)
  message(paste(c("fitted qs =", round(qgamma(p, shape = X$shape, scale = X$scale),3)), collapse = " "))
  message(paste(c("target qs =", c(lower_D, upper_D)), collapse = " "))
  message(paste("fitted var =", round(var_D, 3)))
  c(scale = X$scale, shape = X$shape)
  
}


fit_beta <- function(mean, lower, ci = 0.99) {
  a <- (1 - ci)/2
  p <- c(a , 1 - a)
  
  f <- function(alpha) {
    beta <- alpha*(1 - mean) / mean
    x <- qbeta(p = p,shape1 = alpha, shape2 = beta)
    sum((x[1] - c(lower))^2)
  }
  
  alpha <- optimise(f = f, interval = c(0,1e3), maximum = FALSE)$minimum
  beta <- alpha*(1 - mean) / mean
  
  message(paste(c("fitted qs =", round(qbeta(p, shape1 = alpha, shape2 =  beta),3)), collapse = " "))
  
  c(alpha = alpha, beta = beta[[1]])
  
}


priors_wide_to_long <- function(d) {
  stopifnot(all(xor(is.na(d$shape1), is.na(d$shape))))
  d$type <- ifelse(is.na(d$shape1), "gamma", "beta")
  
  tr <- c(par = "name",
          scale = "gamma_scale",
          shape = "gamma_shape",
          shape1 = "beta_shape1",
          shape2 = "beta_shape2")
  d <- prior_rename(d, names(tr), tr)
  d <- d[c("region", "type", tr)]
  
  extra <- data.frame(
    region = "kabwe",
    type = "null",
    name = "start_date",
    gamma_scale = NA_real_,
    gamma_shape = NA_real_,
    beta_shape1 = NA_real_,
    beta_shape2 = NA_real_,
    stringsAsFactors = FALSE)
  
  d <- rbind(d, extra)
  
  ## Not all parameters are region-specific, let's fix that too:
  f <- function(p) {
    s <- d[d$name == p, ]
    msg <- setdiff(d$region, s$region)
    if (length(msg) > 0) {
      i <- match("england", s$region)
      extra <- s[rep(i, length(msg)), ]
      extra$region <- msg
      s <- rbind(s, extra)
    }
    s
  }
  
  res <- do.call(rbind, lapply(unique(d$name), f))
  
  ## We must have all parameters for all regions, and no doubles
  stopifnot(all(table(res$region, res$name) == 1))
  
  res <- res[order(res$region, res$name), ]
  rownames(res) <- NULL
  res
}
