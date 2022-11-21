# Contains LCO method functions for all examples except SRD

# calculations done within each loop for iIS
iis.loglik <- function(test_mask, Y, mu, P, sig, log = TRUE,
                          Sigma = c("diagonal", "car"), Sigma_objs = list(),
                          glmm = FALSE, dim_theta = 1, X_eff = NULL,
                          Sigma_inv = NULL) {
  Y_i <- Y[test_mask]
  mu_i <- mu[test_mask]
  P_i <- P[test_mask]
  n_i <- length(Y_i)
  dif <- Y_i - mu_i
  if (!glmm && Sigma == "diagonal") {
    if (n_i == 1) {
      ll <- dnorm(Y_i, mean = mu_i, sd = sqrt(unique(P_i) + sig), log = log)
    } else {
      ones <- rep(1, n_i)

      ll <-  -n_i / 2 * log(2 * pi) - sum(P_i) / 2 -
        log(sig) / 2 - log((1 / sig + sum(1 / P_i))) / 2 -
        t(dif) %*% solve(diag(P_i) + sig^2 * ones %*% t(ones)) %*% dif / 2
    }
    yhat <- mu_i
  } else if (glmm) {
    js <- apply(X_eff, 2, \(col) any(col[test_mask] != 0))
    X_eff <- X_eff[, js, drop = FALSE]
    if (is.null(Sigma_inv)) {
      Sigma_inv <- 1 / sig^2
    }
    approx <- lpl_approx_pois(
      Y[!test_mask], X = X_eff[!test_mask, , drop = FALSE],
      beta = rep(0, dim_theta), Sigma_inv = Sigma_inv,
      mu = mu[!test_mask], dist = "poisson"
    )
    beta0 <- approx$beta
    hess_beta <- hess_pois(
      X_eff[!test_mask,  , drop = FALSE], beta0, Y[!test_mask], mu[!test_mask]
    ) - 2 * Sigma_inv

    # second laplce approximation
    hess <- hess_pois(
      X_eff[test_mask, , drop = FALSE], beta0, Y[test_mask], mu[test_mask]
    ) + hess_beta

    ll <- log(2 * pi) - log(diag(-hess)) +
      ll_pois(X_eff[test_mask, , drop = FALSE],
              beta0, Y[test_mask], mu[test_mask]) +
      ll_norm(beta0, beta0, -hess_beta)
    yhat <- as.numeric(exp(mu[test_mask] +
                             X_eff[test_mask, , drop = FALSE] %*% beta0))

  }
  data.frame(
    yhat = yhat,
    loglik = ll
  )
}
# function which takes samples, performs calculations a la iis.loop to get
# log ratios, then inputs them to loo::loo to obtain yhats, DIC, WAIC
iis <- function(loops, mu_samps, X_b, Y, sig_samps, tau_samps = NULL,
                   P = NULL, Sigma = c("diagonal", "car"), glmm = FALSE,
                   X = NULL, dim_theta = 1) {
  # X is without the b_i
  # b_samps is an MCMC x data point matrix, with the b_i for each data point
  # all values are provided as variances
  # either tau_samps or P must be null

  N <- length(Y)
  start_iis <- Sys.time()

  logliks <- purrr::map_df((1:nrow(mu_samps)), function(s) {
    if (is.null(P)) {
      P_t <- rep(tau_samps[s], N)
    } else {
      P_t <- P
    }
    purrr::map_df(unique(loops), function(loop) {
      test_mask <- loop == loops
      iis.loglik(test_mask,
                    Y, mu_samps[s, ],
                    (P_t), sig_samps[s],
                    Sigma = Sigma,
                    glmm = glmm, dim_theta = dim_theta, X_eff = X_b
      ) |>
        mutate(
          i = which(test_mask),
          iter = s,
          loop = loop
        )
    })
  }) |> dplyr::arrange(i, iter)


  # yhat psiis
  logliks_iis <- matrix(logliks$loglik, nrow = nrow(mu_samps),
                        ncol = ncol(mu_samps))
  mu_iis <- matrix(logliks$yhat, nrow = nrow(mu_samps), ncol = ncol(mu_samps))
  r_eff_iis <- loo::relative_eff(exp(-logliks_iis),
                                 chain_id = rep(1:4, each = nrow(mu_samps) / 4))
  psis_obj_iis <- loo::psis(-logliks_iis, r_eff = r_eff_iis, cores = 2)

  yhat <- loo::E_loo(mu_iis, psis_obj_iis,
                     type = "mean",
                     log_ratios = -logliks_iis
  )
  end_iis <- Sys.time()

  # yhat iis
  # draw posterior predictive samples
  yhat_iis <- logliks |>
    dplyr::group_by(i) |>
    dplyr::summarise(
      yhat_mean = mean(yhat * exp(-loglik - log(mean(exp(-loglik)))))
    ) |>
    dplyr::pull(yhat_mean)


  time <- as.double(
    difftime(end_iis, start_iis, units = "secs")) / length(unique(loops))

  df_means <- data.frame(
    loop = loops,
    yhat_psiis_fm = yhat$value,
    yhat_iis_fm = yhat_iis,
    pk_iis_fm = sighat$pareto_k,
    time_iis_fm = time
  )

  df_means |> select(-loop)
}
