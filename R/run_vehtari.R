vehtari.loop <- function(mu, X_r, reffs, Y, log_offset, test_mask,
                         Sig_inv, family, tau_inv = NULL) {
  if (is.null(log_offset)) {
    log_offset <- rep(0, length(Y))
  }
  fixef <- diag(Sig_inv) == 0
  X_test <- X_r[test_mask, , drop = FALSE]
  if (length(tau_inv) > 1) {
    tau_test <- tau_inv[!test_mask]
  } else {
    tau_test <- tau_inv
  }
  lpl_fun <- switch(family,
                    "poisson" = lpl_approx_pois,
                    "normal" = lpl_approx_norm)
  lpl <- lpl_fun(Y[!test_mask], X = X_r[!test_mask, ],
                 beta = reffs, Sigma_inv = Sig_inv[!fixef, !fixef],
                 mu = log_offset[!test_mask] + mu[!test_mask], dist = family,
                 tau_inv = tau_test)
  yhat <- as.vector(X_test %*% lpl$beta + mu[test_mask] + log_offset[test_mask])
  if (family == "poisson") {
    yhat <- exp(yhat)
  }
  data.frame(
    yhat = yhat
  )
}

vehtari <- function(loops, mu_samps, X, Y, sig_samps, tau_samps = NULL,
                    P = NULL, loop = NULL, n_cores = 0,
                    Sigma = c("diagonal", "car"),
                    Sigma_objs = list(), family = "poisson",
                    log_offset = NULL, beta_samps = NULL, glmm = TRUE,
                    fixef_prior_scale = NULL) {
  if (Sigma == "diagonal") {
    xmatch <- match_to_x(names(sig_samps), X)
    reff_mask <- !is.na(xmatch)
  } else if (Sigma == "carar") {
    reff_mask <- rep(TRUE, ncol(X))
    reff_mask[1:4] <- FALSE
  } else if (Sigma == "car") {
    reff_mask <- rep(TRUE, ncol(X))
    reff_mask[1:(ncol(X) - ncol(Sigma_objs$W))] <- FALSE
  }
  reff_samps <- beta_samps[, reff_mask]

  S <- seq_along(sig_samps[[1]]) # number of MCMC iterations
  N <- length(Y)

  start_iis <- Sys.time()
  logliks_y_cp <- purrr::map_df(S, function(s) { # for each MC sample
    if (is.null(P)) {
      P <- rep(tau_samps[s], N)
    }

    if (is.null(tau_samps)) {
      tau_inv <- 1 / P
    } else {
      tau_inv <- 1 / tau_samps[s]
    }
    # create Sigma^-1
    if (Sigma == "diagonal") {
      sigi_diag <- rep(0, ncol(X))
      sigs <- purrr::map_dbl(sig_samps, ~ 1 / .[s])
      sigi_diag[reff_mask] <- sigs[xmatch[reff_mask]]
      Sig_inv <- diag(sigi_diag)
    } else if (Sigma == "car") {
      prec <- (sig_samps$tau[s]) *
        (Sigma_objs$D - sig_samps$alpha[s] * Sigma_objs$W)
      fixed <- ncol(X) - ncol(prec)
      Sig_inv <- matrix(0, ncol(X), ncol(X))
      Sig_inv[-c(1:fixed), -c(1:fixed)] <- prec
    } else if (Sigma == "carar") {

      Sig_inv <- stcarar_sigmainv(
        Sigma_objs$rho[s, "rho.S"],
        Sigma_objs$rho[s, "rho.T"],
        Sigma_objs$W, (Sigma_objs$tau2[s])
      )
    }

    if (is.null(loop)) {
      purrr::map_df(unique(loops), \(loop) {
        test_mask <- loop == loops
        vehtari.loop(mu_samps[s, ], X[, reff_mask], reff_samps[s, ],
                     Y, log_offset, test_mask, Sig_inv, family = family,
                     tau_inv = tau_inv) |>
          mutate(
            idx = which(test_mask),
            loop = loop
          )
      })
    } else {
      test_mask <- loop == loops
      vehtari.loop(mu_samps[s, ], X[, reff_mask], reff_samps[s, ],
                   Y, log_offset, test_mask, Sig_inv, family = family,
                   tau_inv = tau_inv) |>
        mutate(
          idx = which(test_mask),
          loop = loop
        )
    }
  })

  res <- logliks_y_cp |> group_by(loop, idx) |>
    summarise(yhat_vt = mean(yhat))
  end_iis <- Sys.time()

  select(res, -loop) |>
    mutate(
      time_vt = difftime(end_iis, start_iis, units = "secs") /
        length(unique(loops))
    ) |>
    arrange(idx)
}
