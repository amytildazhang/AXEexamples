##### helper functions for MCV
.collate_folds <- function(objs) {
  obj_names <- names(objs[[1]])

  loops <- purrr::map_lgl(objs, ~ length(.$idx) > 1)
  if (!all(loops)) {
    i <- 1
  } else {
    i <- which(loops)[1]
  }

  purrr::map(obj_names, function(obnm) {
    ex <- objs[[i]][[obnm]]

    if (checkmate::test_atomic_vector(ex)) {
      op <- "c"
    } else {
      op <- "rbind"
    }

    bind <- do.call(op, purrr::map(objs, ~ .[[obnm]]))
    names(bind) <- NULL

    bind
  }) %>% purrr::set_names(obj_names)
}

.posterior_mode <- function(vec) {
  density_estimate <- density(vec)
  density_estimate$x[which.max(density_estimate$y)]
}


posterior_mode <- function(samps) {
  if (is.null(dim(samps))) {
    .posterior_mode(samps)
  } else {
    apply(samps, 2, function(parsamps) {
      .posterior_mode(parsamps)
    })
  }
}


###### functions within each fold

crossval_eight <- function(loop, data, loops, X, seed, cmdmod, ...) {
  a <- Sys.time()
  test_mask <- loop == loops
  cv_dat <- data

  cv_dat$train_idx <- which(!test_mask)
  cv_dat$test_idx <- as.array(which(test_mask))
  cv_dat$H <- sum(test_mask)

  fit_cv <- cmdmod$sample(data = cv_dat,
                          refresh = 0,
                          seed = seed,
                          ...)



  cv_betas <- fit_cv$draws(variables = c("mu", "theta"),
                           format = "draws_matrix")
  cv_samps <- X %*% t(cv_betas)
  b <- Sys.time()
  cv_taus <- apply(fit_cv$draws(variables = "tau"), 3, c)


  stansummary <- fit_cv$summary()

  map_idx <- which.max(fit_cv$draws(variables = "lp__",
                                    format = "draws_matrix"))

  yhat_sij <- cv_betas[map_idx, "mu"]
  list(
    idx = which(test_mask), cv_yhat = mean(cv_samps[test_mask, ]),
    sighat_mcv = mean(cv_taus),
    time_cv = as.numeric(difftime(b, a, units = "secs")),
    max_rhat = max(stansummary$rhat),
    mean_rhat = mean(stansummary$rhat),
    yhat_cv_sij = as.vector(yhat_sij),
    min_ess_bulk = min(stansummary$ess_bulk),
    mean_ess_bulk = mean(stansummary$ess_bulk),
    min_ess_tail = min(stansummary$ess_tail),
    mean_ess_tail = mean(stansummary$ess_tail)

  )
}

crossval_radon <- function(loop, data, loops, X, formula, seed, ...) {
  a <- Sys.time()
  test_mask <- loop == loops
  rfit <- rstanarm::stan_glmer(
    formula,
    data = dplyr::filter(data, !test_mask),
    seed = seed, refresh = 0, ...
  )
  yhat <- rstanarm::posterior_predict(rfit,
                                      newdata = dplyr::filter(data, test_mask),
                                      draws = 1e3
  )
  b <- Sys.time()
  stopifnot(ncol(yhat) == sum(test_mask))
  samples <- rstan::extract(rfit$stanfit)

  tau <- mean(samples$aux)
  sig2 <- mean(samples$theta_L)

  Y_t <- data$log_radon[test_mask]
  if ("beta" %in% names(samples)) {
    fixef <- cbind(samples$alpha, samples$beta)
  } else {
    fixef <- samples$alpha
  }
  coeff_samps <- cbind(fixef, samples$b)

  # since test county is not provided, it is present in
  # samples as last column in samples$b
  test_col <- which(apply(X, 2, function(col) {
    all(col == as.numeric(test_mask))
  }))
  ind <- X[, test_col]
  X <- cbind(X[, -test_col], ind)

  mean <- coeff_samps %*% t(X)
  elpd <- sapply(1:sum(test_mask), function(i) {
    log(mean(dnorm(Y_t[i], mean = mean[, i], sd = samples$aux)))
  })

  map_idx <- which.max(samples$lp__)
  yhat_sij <- posterior_mode(yhat)
  # posterior mode
  data.frame(
    idx = which(test_mask),
    yhat_cv = apply(yhat, 2, mean),
    yhat_cv_sij = yhat_sij,
    tauhat_cv = tau,
    time_cv = as.double(difftime(b, a, units = "secs")),
    sighat_cv = sqrt(sig2),
    elpd_cv = elpd
  )
}

crossval_radon_simul <- function(mod_no, dat, n, perc,
                                 subset_idx, info, seed, ...) {
  mod <- info$models_mcv[[mod_no]]

  test_mask <- dat$county == info$c_county
  train_mask <- !test_mask

  # fit to just training data
  start_time <- Sys.time()
  cv_fit <- rstanarm::stan_glmer(
    mod,
    data = dat[train_mask, ],
    seed = seed,
    refresh = 0, iter = 1e4, thin = 5, ...
  )
  yhat_cv <- rstanarm::posterior_predict(
    cv_fit,
    newdata = dplyr::filter(dat, test_mask),
    draws = 1e3
  )

  end_time <- Sys.time()


  cv_samps <- rstan::extract(cv_fit$stanfit)
  sf <- radon_simul_samples(perc, n, subset_idx, mod_no, cv = TRUE)
  readr::write_rds(cv_samps, sf)


  samps <- as.array(cv_fit$stanfit)

  # log lik
  elpd <- rstanarm::log_lik(cv_fit, newdata = dplyr::filter(dat, test_mask)) %>%
    exp() %>%
    apply(2, mean) %>%
    log()



  X <- mm_radon(mod_no, dat, info$models_axe[[mod_no]])

  county_col <- which(apply(X, 2, function(col) {
    all(col == as.numeric(test_mask))
  }))
  X <- cbind(X[, -county_col], X[county_col])
  colnames(X) <- c(colnames(X)[-ncol(X)], sprintf("county%s", info$c_county))
  coeffs <- coeffs_radon(cv_samps, cv = TRUE)

  tau_samps <- cv_samps$aux
  lpd <- do.call("cbind", purrr::map((1:nrow(coeffs)), function(mit) {
    dnorm(
      dat$log_radon[test_mask],
      mean = X[test_mask, ] %*% t(coeffs[mit, , drop = FALSE]),
      sd = tau_samps[mit]
    )
  })) %>%
    apply(1, mean) %>%
    log()

  # cond_lp

  sigma2 <- mean(cv_samps$theta_L)
  sig_inv <- rep(0, ncol(X))
  sig_inv[stringr::str_detect(colnames(X), "county")] <- 1 / sigma2
  cond_lpd <- axe_ll(diag(sig_inv), mean(tau_samps),
                     X, test_mask, dat$log_radon)
  cv_ycond <- axe(diag(sig_inv), mean(tau_samps)^2,
                  X, dat$county == info$c_county, dat$log_radon)

  data.frame(
    y = dat$log_radon[test_mask],
    rst_elpd = elpd,
    elpd_cv = lpd,
    elpd_cv_cond = cond_lpd,
    yhat_cv_sij = posterior_mode(yhat_cv),
    cv_yhat = apply(yhat_cv, 2, mean),
    cv_ycond = cv_ycond,
    tau2_cv = mean(tau_samps)^2,
    sig2_cv = sigma2,
    max_rhat = max(cv_fit$stan_summary[, "Rhat"]),
    min_rhat = min(cv_fit$stan_summary[, "Rhat"])
  ) %>%
    dplyr::mutate(
      time_cv = as.double(difftime(end_time, start_time, units = "secs"))
    )
}




crossval_lol <- function(loop, data, loops, seed, X, ...) {
  test_mask <- loop == loops
  a <- Sys.time()
  cv_fit <- rstanarm::stan_glmer(
    kills ~ position + team + (1 | player) + (1 | champion) +
      log_dpm + log_egpm,
    family = "poisson", data = dplyr::filter(data, !test_mask), refresh = 0,
    seed = seed, refresh = 0, ...
  )


  yhats <- rstanarm::posterior_predict(
    cv_fit,
    newdata = dplyr::filter(data, test_mask)
  )
  b <- Sys.time()

  data.frame(
    time_cv = as.double(difftime(b, a, units = "secs")),
    loop = loop,
    idx = which(test_mask),
    yhat_cv = apply(yhats, 2, mean),
    yhat_cv_sij = posterior_mode(yhats)
  )
}


crossval_slc <- function(loop, info, loops, seed, cmdmod, ...) {
  message(loop)
  test_mask <- loop == loops
  a <- Sys.time()

  D <- diag(rowSums(info$W))
  W_rr <- info$W[!test_mask, !test_mask]
  W_ir <- info$W[test_mask, !test_mask]

  X_cv <- rbind(info$X[!test_mask, ], info$X[test_mask, , drop = FALSE])
  W_cv <- cbind(
    rbind(W_rr, W_ir),
    c((W_ir), 0)
  )
  loff_cv <- c(info$log_offset[!test_mask], info$log_offset[test_mask])
  sdat_cv <- list(
    n = info$n,
    p = info$p, X = X_cv, y_train = info$y[!test_mask],
    log_offset = loff_cv, W = W_cv, test_idx = info$n,
    D_sparse = rowSums(W_cv), W_n = info$W_n
  )
  cv_fit <- cmdmod$sample(data = sdat_cv, iter_warmup = 5e3,
                          refresh = 0,
                          iter_sampling = 5e3, thin = 5)

  yhats <- apply(cv_fit$draws(), 3, c)[, "y_test"]

  b <- Sys.time()




  logliks <- sapply(yhats, function(y) {
    dpois(info$y[test_mask], yhats)
  }) %>%
    mean() %>%
    log()


  data.frame(
    time_cv = as.double(difftime(b, a, units = "secs")),
    loop = loop,
    idx = which(test_mask),
    yhat_cv = mean(yhats),
    yhat_cv_sij = posterior_mode(yhats),
    elpd_cv = logliks
  )
}



crossval_air <- function(loop, info, loops, seed, formula, ...) {
  message(loop)
  test_mask <- loop == loops
  a <- Sys.time()

  set.seed(seed)
  cvdat <- info$df
  cvdat$observed[test_mask] <- NA

  chain1 <- CARBayesST::ST.CARar(
    formula = formula, family = "poisson",
    data = cvdat, W = info$W, burnin = 20000, n.sample = 220000,
    thin = 100, AR = 1
  )
  chain2 <- CARBayesST::ST.CARar(
    formula = formula, family = "poisson",
    data = cvdat, W = info$W, burnin = 20000, n.sample = 220000,
    thin = 100, AR = 1
  )
  chain3 <- CARBayesST::ST.CARar(
    formula = formula, family = "poisson",
    data = cvdat, W = info$W, burnin = 20000, n.sample = 220000,
    thin = 100, AR = 1
  )

  samples <- purrr::map(names(chain1$samples), function(param) {
    rbind(
      chain1$samples[[param]],
      chain2$samples[[param]],
      chain3$samples[[param]]
    )
  }) %>% purrr::set_names(names(chain1$samples))




  yhats <- samples$Y

  logliks <- sapply(1:sum(test_mask), function(test_idx) {
    log(mean(dpois(info$df$observed[test_mask][test_idx], yhats[, test_idx])))
  })
  b <- Sys.time()

  data.frame(
    time_cv = as.double(difftime(b, a, units = "secs")),
    loop = loop,
    idx = which(test_mask),
    yhat_cv = apply(yhats, 2, mean),
    yhat_cv_sij = posterior_mode(yhats),
    elpd_cv = logliks
  )
}


gp_cmvn <- function(a, alpha, rho, f, X, test_mask) {
  sigma <- cov_exp_quad(X, alpha, rho)

  sig_train_inv <- solve(sigma[!test_mask, !test_mask])
  sig12 <- sigma[test_mask, !test_mask, drop = FALSE]
  sig11 <- sigma[test_mask, test_mask]


  # draw fi
  mean <- a + sig12 %*% sig_train_inv %*% f
  sd <- sqrt(sig11 - sig12 %*% sig_train_inv %*% t(sig12))

  list(mean = mean, sd = sd)

}

predict_boston <- function(a, alpha, rho, f, X, test_mask, sigma) {
  pars <- gp_cmvn(a, alpha, rho, f, X, test_mask)
  # draw fi
  fi <- rnorm(1, mean = pars$mean, sd = pars$sd)

  rnorm(1, fi, sigma)

}

predict_sonar <- function(a, alpha, rho, f, X, test_mask) {
  # sigma
  pars <- gp_cmvn(a, alpha, rho, f, X, test_mask)

  # draw fi
  fi <- rnorm(1, mean = pars$mean, sd = pars$sd)

  # draw from poisson
  rbinom(1, 1, pnorm(fi))

}

crossval_sonar <- function(loop, info, loops, seed, cmdmod, ...) {
  message(loop)
  tryCatch({

    test_mask <- loop == loops

    X <- t(standardize_x_sonar(t(info$X), test_mask))
    set.seed(seed)
    a <- Sys.time()

    cvinfo <- list(
      X = t(X)[!test_mask, ],
      y = info$y[!test_mask],
      N = info$N - sum(test_mask),
      D = info$D
    )
    cvfit <- cmdmod$sample(data = cvinfo,
                           refresh = 0,
                           seed = seed,
                           ...)

    samples <- cvfit$draws(variables = c("a", "alpha", "rho"))
    as <- as.vector(samples[, , 1])
    alphas <- as.vector(samples[, , 2])
    rhos <- as.vector(samples[, , 3])
    f <-  apply(cvfit$draws(variables = "f"), 3, c)


    yhats <- sapply(seq_along(as), function(s) {
      predict_sonar(as[s], alphas[s], rhos[s], f[s, ], X, test_mask)
    })


    b <- Sys.time()


    data.frame(
      time_cv = as.double(difftime(b, a, units = "secs")),
      loop = loop,
      idx = which(test_mask),
      yhat_cv = mean(yhats)
    )
  },
  error = function(cond) {
    message(cond)
    return(NA)
  })

}

crossval_boston <- function(loop, info, loops, seed, cmdmod, ...) {
  message(loop)
  tryCatch({

    test_mask <- loop == loops

    X <- t(standardize_x_sonar(t(info$X), test_mask))
    set.seed(seed)
    a <- Sys.time()

    cvinfo <- list(
      X = t(X)[!test_mask, ],
      y = info$y[!test_mask],
      N = info$N - sum(test_mask),
      D = info$D
    )
    cvfit <- cmdmod$sample(data = cvinfo,
                           refresh = 0,
                           seed = seed,
                           ...)

    samples <- apply(cvfit$draws(variables = c("alpha", "rho", "sigma", "f")),
                     3, c)

    yhats <- sapply(seq_along(as), function(s) {
      predict_boston(0, samples[s, "alpha"], samples[s, "rho"],
                     samples[s, 4:ncol(samples)], X, test_mask,
                     samples[s, "sigma"])
    })


    b <- Sys.time()


    data.frame(
      time_cv = as.double(difftime(b, a, units = "secs")),
      loop = loop,
      idx = which(test_mask),
      yhat_cv = mean(yhats)
    )
  },
  error = function(cond) {
    message(cond)
    return(NA)
  })

}

#####

# "..." captures parameters passed to underlying STAN model.
#' Manual cross-validation for AXE paper examples
#'
#' `mcv_*()` are high-level functions which obtain posterior predictive
#' samples Y_j | Y_-j for each cross-validation loop
#'
#' @param info All. List of data to use for CV folds; assumes list is created by
#'     `prep_*()`.
#' @param n_cores All. Number of cores to use, if parallelizing. Default value
#'   varies.
#' @param seed Random seed generator passed to model sampling method.
#' @param ... Additional arguments passed to modeling method.
#'
#' @return Dataframe or list of dataframes with ground-truth MCV values.
#'
#'  @export
mcv_eight <- function(info = eight, n_cores = 7,
                      data_scale = seq(0.1, 4, by = 0.1),
                      seed = 238972, ...) {
  X <- model.matrix(
    ~school,
    data = info$df,
    contrasts.arg = list(school = contrasts(info$df$school, contrasts = FALSE))
  )
  cmdmod <- cmdstanr::cmdstan_model("data-raw/stan/eightschools_sim.stan")
  loops <- 1:8
  purrr::map(data_scale, function(scl) {
    scld_dat <- scale_dat(info$stan_list, scl)
    if (any(round(scl -  c(2.2, 2.4, 2.8, 4), digits = 1) == 0)) {
      n_warmup <- 4000
      n_iter <- 2000
      thin <- 2
    } else {
      n_warmup <- 1000
      n_iter <- 1000
      thin <- 1
    }
    if (n_cores > 1) {
      cl <- snow::makeCluster(n_cores, outfile = "")
      snow::clusterExport(cl, list("posterior_mode", ".posterior_mode"),
                          envir = environment())

      yhats <- snow::parLapply(
        cl, loops, crossval_eight,
        scld_dat, loops, X, seed, cmdmod, iter_warmup = n_warmup,
        iter_sampling = n_iter, thin = thin
      )

      snow::stopCluster(cl)
    } else {
      yhats <- lapply(
        loops, crossval_eight,
        scld_dat, loops, X, seed, cmdmod, iter_warmup = n_warmup,
        iter_sampling = n_iter, thin = thin)
    }
    objs <- .collate_folds(yhats)
    c(objs, data_scale = scl)
  })
}



#' @describeIn mcv_eight
#'
#'  @export
mcv_radon_full <- function(info = radon_1, models = radon_1$models_mcv,
                           seed = 32897, n_cores = 7, ...) {
  xdat <- info$data %>%
    dplyr::mutate(
      floor = factor(floor), county = factor(county)
    )
  loops <- as.character(info$data$county)
  y <- info$data$log_radon

  purrr::map(seq_along(models), function(mod_no) {
    print(mod_no)
    X <- model.matrix(info$models_axe[[mod_no]],
                      data = xdat,
                      contrasts = contrasts_for_pooling(xdat)
    )

    if (n_cores > 1) {
      cl <- snow::makeCluster(n_cores, outfile = "")
      snow::clusterExport(cl, list("X", "%>%",
                                   "posterior_mode", ".posterior_mode",
                                   "radon2_prefix"),
                          envir = environment()
      )

      objs <- snow::parLapply(
        cl, unique(loops), crossval_radon,
        info$data, loops,
        X, models[[mod_no]], seed, ...
      )

      snow::stopCluster(cl)
    } else {
      objs <- lapply(
        unique(loops), crossval_radon,
        info$data, loops,
        X, models[[mod_no]], seed
      )
    }
    objs <- dplyr::bind_rows(objs)

    objs %>%
      dplyr::arrange(idx) %>%
      dplyr::mutate(
        loop = as.numeric(factor(loops))
      )
  })
}


#' @describeIn mcv_eight
#'
#'  @export
mcv_radon_simul <- function(info = radon_2,
                            seed = 9871, n_cores = 5, ...) {
  loop_over_radon_simul(
    info = info, seed = seed, n_cores = n_cores,
    export = list("mm_radon", "coeffs_radon", "axe_ll", "axe",
                  "axe_v", "radon_simul_samples",
                  "posterior_mode", ".posterior_mode",
                  "radon2_prefix"),
    expenv = environment(),
    FUN = crossval_radon_simul,

    ...
  )
}


######




#' @describeIn mcv_eight
#'
#'  @export
mcv_lol <- function(info = lol,
                    seed = 23982, n_cores = 7, ...) {
  loops <- info$loops

  if (n_cores > 1) {
    cl <- snow::makeCluster(n_cores, outfile = "")
    snow::clusterExport(cl, list("%>%", "posterior_mode", ".posterior_mode"),
                        envir = environment()
    )

    objs <- snow::parLapply(
      cl, unique(loops), crossval_lol,
      info$data, loops, seed, info$X, ...
    )

    snow::stopCluster(cl)
  } else {
    objs <- lapply(
      unique(loops), crossval_lol,
      info$data, loops, seed, info$X, ...
    )
  }

  objs
}


#' @describeIn mcv_eight
#'
#'  @export
mcv_slc <- function(info = slc, seed = 29873, n_cores = 3, ...) {
  loops <- info$data$loops

  cmdmod <- cmdstanr::cmdstan_model("data-raw/stan/sparsecar_slp_cv.stan")
  if (n_cores > 1) {
    cl <- snow::makeCluster(n_cores, outfile = "")
    snow::clusterExport(cl, list("stanmodels", "posterior_mode",
                                 ".posterior_mode"),
                        envir = environment()
    )

    objs <- snow::parLapply(
      cl, unique(loops), crossval_slc,
      info$data, loops, seed, cmdmod, ...
    )

    snow::stopCluster(cl)
  } else {
    objs <- lapply(
      unique(loops), crossval_slc,
      info$data, loops, seed, cmdmod, ...
    )
  }

  dplyr::bind_rows(objs)
}



#' @describeIn mcv_eight
#'
#'  @export
mcv_air <- function(info = air, seed = 69872, n_cores = 3, ...) {
  loops <- info$data$loops

  if (n_cores > 1) {
    cl <- snow::makeCluster(n_cores, outfile = "")
    snow::clusterExport(cl, list("posterior_mode", ".posterior_mode"),
                        envir = environment()
    )

    objs <- snow::parLapply(
      cl, unique(loops), crossval_air,
      info$data, loops, seed, info$data$formula, ...
    )

    snow::stopCluster(cl)
  } else {
    objs <- lapply(
      unique(loops), crossval_air,
      info$data, loops, seed, info$data$formula, ...
    )
  }

  dplyr::bind_rows(objs)
}
