
#' Obtain LCO approximations for AXE paper examples
#'
#' `pfit_*()` fits models for the full data of each example in the
#'     AXE paper and obtains LCO approximations based on posterior samples.
#'
#'  @export
pfit_eight <- function(info = eight, seed = 238972, nodes = 6, ...) {
  mod <- get_cmdstan_obj("eight")
  purrr::map_df((info$data_scale), function(scl) {
    message(scl)
    scld_dat <- scale_dat(info$stan_list, scl)

    # fit to full data
    pfit <- mod$sample(data = scld_dat, seed = seed, refresh = 0, ...)
    samples <- pfit$draws(format = "draws_matrix")

    X <- model.matrix(
      ~school, data = info$df,
      contrasts.arg = list(school = contrasts(info$df$school,
                                              contrasts = FALSE))
    )


    #
    # Integrated Importance sampling
    # (Li, Qiu, Zhang, Feng, 2016, Statistics and Computing)
    #
    iis_og <- iis(
      1:8, do.call("cbind", lapply(1:8, function(x) samples[, "mu"])),
      X_b = X[, -1], Y = scld_dat$y,
      sig_samps = samples[, "tau"]^2, P = info$df$sd^2
    )


    #
    #  Vehtari's method (Vehtari et al, 2016, JMLR)
    #
    beta_cols <- stringr::str_detect(colnames(samples), "theta") |
      stringr::str_detect(colnames(samples), "mu")

    beta_samps <- unclass(samples[, beta_cols])
    Xbeta_samps <- t(info$X %*% t(beta_samps))
    vt <- vehtari(
      1:8, mu_samps = Xbeta_samps, X, Y = scld_dat$y,
      sig_samps = list("school" = unclass(samples[, "tau"]^2)),
      P = info$df$sd^2, Sigma = "diagonal", family = "normal",
      log_offset = NULL, beta_samps = beta_samps, glmm = FALSE,
      fixef_prior_scale = NULL
    ) |> select(-idx, -loop)


    #
    # PSIS-LOO elpd and yhats (Vehtari, Gelman, Gabry 2017)
    #
    set.seed(seed)
    start <- Sys.time()
    logliks <- sapply((1:nrow(info$df)), function(j) {
      dnorm(scld_dat$y[j], mean = Xbeta_samps[, j], sd = info$df$sd[j],
            log = TRUE)
    })
    r_eff <- loo::relative_eff(exp(-logliks), chain_id = rep(1:4, each = 1000))
    psis_obj <- loo::psis(-logliks, r_eff = r_eff, cores = 2)
    psis_yhats <- loo::E_loo(Xbeta_samps, psis_obj,
                             type = "mean",
                             log_ratios = -logliks
    )
    end <- Sys.time()


    #
    # Ghosting
    #
    ghst <- ghosting_c(
      seed, 1:8, do.call("cbind", lapply(1:8, function(x) samples[, "mu"])),
      samples[, "tau"]^2, X[, -1], scld_dat$y,
      P = info$df$sd^2
    )


    # save data
    data.frame(
      tau_hat = (mean(samples[, "tau"]^2)),
      yhat_psis_c = psis_yhats$value,
      data_scale = scl,
      school_idx = 1:8
    ) |>
      dplyr::bind_cols(
        iis_og, vt, ghst
      ) |>
      dplyr::mutate(
        i = 1:dplyr::n(),
        loop = 1:8
      )
  })
}



#' @describeIn pfit_eight
#' @export
pfit_radon_full <- function(info = radon_1, seed = 32897, ...) {
  xdat <- info$data |>
    dplyr::mutate(floor = factor(floor), county = factor(county))

  purrr::map(seq_along(info$models_mcv), function(mod_no) {
    message(sprintf("Model %s", mod_no))
    mod <- info$models_mcv[[mod_no]]

    rfit <- rstanarm::stan_glmer(mod,
                                 data = info$data, seed = seed,
                                 refresh = 0, ...
    )

    X <- model.matrix(info$models_axe[[mod_no]],
                      data = xdat,
                      contrasts = contrasts_for_pooling(xdat)
    )

    samples <- rstan::extract(rfit$stanfit)
    coeffs <- coeffs_radon(samples)
    loops <- info$data$county

    fixef <- !stringr::str_detect(colnames(X), "county")
    mu_samps <- coeffs[, fixef] %*% t(X[, fixef])


    #
    # iIS
    #
    message(sprintf("Running iIS"))
    iis_og <- iis(
      loops, mu_samps,
      X_b = X[, !fixef], Y = info$data$log_radon,
      sig_samps = samples$theta_L, tau_samps = samples$aux^2, Sigma = "diagonal"
    )

    #
    #  Vehtari
    #
    Xbeta_samps <- coeffs %*% t(X)
    message("Running Vehtari")
    vt <- vehtari(loops, Xbeta_samps,
                  X = X, Y = info$data$log_radon,
                  sig_samps = list("county" = samples$theta_L),
                  tau_samps = samples$aux^2, n_cores = 4, family = "normal",
                  Sigma = "diagonal", glmm = FALSE, beta_samps = coeffs
    )


    message("Ghosting")
    #
    #  Ghosting
    #
    ghst <- ghosting_c(
      seed, loops,
      mu_samps = mu_samps,
      sig_samps = samples$theta_L, X_b = X[, !fixef],
      Y = info$data$log_radon, tau_samps = samples$aux^2, Sigma = "diagonal"
    )

    # posterior samples
    data.frame(
      i = (1:nrow(X)),
      tau_hat = (mean(samples$aux)),
      sighat = mean(sqrt(samples$theta_L)),
      tau_map = as.numeric(bayestestR::map_estimate(samples$aux)),
      sig2_map = as.numeric(bayestestR::map_estimate(samples$theta_L)),
      model = mod_no,
      loop = loops
    ) |>
      dplyr::bind_cols(
        iis_og, vt |> select(-idx, -loop), ghst
      )
  })
}


#' @describeIn pfit_eight
#' @export
pfit_radon_simul <- function(info = radon_2,
                             seed = 9871, n_cores = 5, ...) {
  # Unlike all other pfit_ function, the radon simulated subsets example is so
  # large that this function only saves full-data posterior information
  loop_over_radon_simul(
    info, seed, n_cores,
    export = list("use_saved", "mm_radon", "coeffs_radon", "info",
                  "radon_simul_samples", "savefolder"),
    expenv = environment(),
    FUN = function(mod_no, dat, n, perc, subset_idx, info, seed, ...) {
      rfit <- rstanarm::stan_glmer(
        info$models_mcv[[mod_no]],
        data = dat, seed = seed, refresh = 0, iter = 4e3, ...)
      samples <- rstan::extract(rfit$stanfit)
      readr::write_rds(samples, sf)

      # get fitted values
      X <- mm_radon(mod_no, dat, info$models_axe[[mod_no]])
      coeffs <- coeffs_radon(samples)
      fv <- apply(coeffs %*% t(X), 2, mean)

      # save X
      dat <- dplyr::mutate(dat, floor = factor(floor), county = factor(county))
      X <- model.matrix(info$models_axe[[mod_no]],
                        data = dat,
                        contrasts = contrasts_for_pooling(dat)
      )

      test_mask <- dat$county == info$c_county
      cv_loop <- rep(1, sum(test_mask))


      X_test <- X[test_mask, , drop = FALSE]
      Y_test <- dat$log_radon[test_mask]
      Xbeta_samps <- coeffs %*% t(X_test)
      data.frame(
        yhat_post = as.vector(fv[test_mask]),
        tau2hat_post = mean(samples$aux)^2,
        sig2hat_post = mean(samples$theta_L),
        n_tot = nrow(dat),
        tau2hat_map = as.vector(bayestestR::map_estimate(samples$aux)^2),
        sig2hat_map = as.vector(bayestestR::map_estimate(samples$theta_L))
      )
    })
}


#' @describeIn pfit_eight
#'
#'  @export
pfit_lol <- function(info = lol, seed = 18973, ...) {
  pfit <- rstanarm::stan_glmer(
    kills ~ position + team + log_dpm + log_egpm +
      (1 | player) + (1 | champion),
    family = "poisson", data = info$data, seed = seed
  )
  set.seed(seed)
  phats <- apply(2, rstanarm::posterior_predict(pfit), mean)
  P <- 1 / phats

  samples <- rstan::extract(pfit$stanfit)
  sigs <- samples$theta_L |> apply(2, mean)

  coeffs <- cbind(samples$alpha, samples$beta, samples$b)
  coeffs <- coeffs[, -c(109, ncol(samples$b))] # remove "new" effs
  Xbeta_samps <- coeffs %*% t(info$X)

  X_effs <- info$X[, stringr::str_detect(colnames(info$X), "player")]
  mu_samps <- Xbeta_samps -
    samples$b[, -c(1:109, ncol(samples$b))] %*% t(X_effs)

  message("iIS")
  iis_og <- iis(
    info$loops, mu_samps,
    X_b = X_effs, Y = info$data$kills,
    sig_samps = samples$theta_L[, 2], P = P, Sigma = "diagonal", glmm = TRUE
  )



  vt <- vehtari(
    info$loops, mu_samps = mu_samps, info$X, Y = info$data$kills,
    sig_samps =  list(
      "champion" = samples$theta_L[, 1],
      "player" = samples$theta_L[, 2]
    ), P = P, Sigma = "diagonal",
    family = "poisson",
    log_offset = NULL, beta_samps = coeffs, glmm = TRUE,
    fixef_prior_scale = NULL)

  list(
    pfit = pfit, phats = phats, P = P, sig_champion = sigs[1],
    sig_player = sigs[2],
    df = data.frame(
      idx = seq_along(Y),
      yhat_post = phats
    ) |> dplyr::bind_cols(iis_og, vt)
  )
}

car_inv <- function(D, W, alpha, tau) {
  tau^2 * (D - alpha * W)
}

#' @export
sparsecar_slp <- function(...) {
  rstan::sampling(stanmodels$sparsecar_slp, ...)
}

#' @describeIn pfit_eight
#'
#'  @export
pfit_slc <- function(info = slc$data, seed = 29873, ...) {
  sdat <- list(
    n = info$n, p = info$p, X = info$X, y = info$y,
    log_offset = info$log_offset, W_n = info$W_n, W = info$W,
    D_sparse = rowSums(info$W)
  )

  cmdmod <- get_cmdstan_obj("slc")
  pfit <- cmdmod$sample(
    data = sdat, refresh = 0, seed = seed,
    iter_warmup = 2e3, iter_sampling = 2e3, ...
  )
  set.seed(seed)
  beta_samps <- as.matrix(pfit$draws(variables = "beta",
                                     format = "draws_matrix"))

  D <- diag(rowSums(info$W))
  Xbeta <- info$X %*% apply(beta_samps, 2, mean)

  phi_samps <- as.matrix(pfit$draws(variables = "phi",
                                    format = "draws_matrix"))
  phats <- Xbeta + apply(phi_samps, 2, mean)
  muhats <- exp(phats + info$log_offset)


  P <- as.vector(1 / muhats)
  n <- length(info$y)


  X <- cbind(info$X, diag(rep(1, length(Y))))
  coeffs <- cbind(beta_samps, phi_samps)
  Xbeta_samps <- coeffs %*% t(X)

  alpha_samps <- as.matrix(pfit$draws(variables = "alpha",
                                      format = "draws_matrix"))
  a <- Sys.time()
  spat_mu <- sapply(seq_along(info$y), function(i) {
    Nbhd_i <- which(info$W[i, ] != 0)

    alpha_samps * rowSums(phi_samps[, Nbhd_i, drop = FALSE]) / length(Nbhd_i)
  })

  mu_samps <- apply(beta_samps %*% t(info$X) + spat_mu, 1, function(samp) {
    samp + info$log_offset
  }) |> t()

  b <- Sys.time()

  tau_samps <- as.vector(pfit$draws(variables = "tau",
                                    format = "draws_matrix"))

  message("iIS")
  iis_mc <- iis(
    info$loops, mu_samps,
    X_b = diag(rep(1, length(info$y))),
    Y = info$y, sig_samps = 1 / tau_samps, P = P, Sigma = "car", glmm = TRUE
  )


  #
  #  Integrated Importance sampling (Vehtari et al, 2016, JMLR)
  #

  message("VT")
  vt <- vehtari(
    info$loops, Xbeta_samps, X = X, Y = info$y,
    sig_samps = list("alpha" = alpha_samps, "tau" = tau_samps^2),
    P = P, Sigma = "car", glmm = TRUE,
    Sigma_objs = list(D = diag(rowSums(info$W)), W = info$W),
    beta_samps = coeffs, log_offset = info$log_offset,
    fixef_prior_scale = 1,
    family = "poisson")

  list(
    P = P,
    tauhat_post = mean(tau_samps),
    alpha_post = mean(alpha_samps),
    muhat = Xbeta + apply(phi_samps, 2, mean),
    mutime = difftime(b, a, units = "secs"),
    df = data.frame(
      idx = seq_along(info$y),
      yhat_post = muhats
    ) |> dplyr::bind_cols(iis_mc, vt)
  )
}

#' @describeIn pfit_eight
#'
#'  @export
pfit_air <- function(info = air$data, seed = 69872, ...) {
  set.seed(seed)
  formula <- info$formula

  chain1 <- CARBayesST::ST.CARar(
    formula = formula, family = "poisson",
    data = info$df, W = info$W, burnin = 20000, n.sample = 220000,
    thin = 200, AR = 1
  )
  chain2 <- CARBayesST::ST.CARar(
    formula = formula, family = "poisson",
    data = info$df, W = info$W, burnin = 20000, n.sample = 220000,
    thin = 200, AR = 1
  )
  chain3 <- CARBayesST::ST.CARar(
    formula = formula, family = "poisson",
    data = info$df, W = info$W, burnin = 20000, n.sample = 220000,
    thin = 200, AR = 1
  )
  chain4 <- CARBayesST::ST.CARar(
    formula = formula, family = "poisson",
    data = info$df, W = info$W, burnin = 20000, n.sample = 220000,
    thin = 200, AR = 1
  )

  samples <- purrr::map(names(chain1$samples), function(param) {
    rbind(
      chain1$samples[[param]],
      chain2$samples[[param]],
      chain3$samples[[param]],
      chain4$samples[[param]]
    )
  }) |> purrr::set_names(names(chain1$samples))

  samples <- rstan::extract(pfit)
  D <- diag(rowSums(info$W))
  Xbeta <- info$X %*% apply(samples$beta, 2, mean)

  phats <- Xbeta + apply(samples$phi, 2, mean)
  muhats <- exp(phats + info$log_offset)



  rho <- apply(samples$rho, 2, mean)
  tau2 <- mean(samples$tau2)
  P <- as.vector(1 / muhats)
  loops <- info$loops
  J <- length(unique(loops))
  n <- length(loops)

  Xbeta_samps <- samples$beta %*% t(info$X)
  X <- info$X

  sample_idx <- sample((1:nrow(Xbeta_samps)), 1000)
  cl <- snow::makeCluster(6)
  snow::clusterExport(
    cl,
    list("X", "stcarar_sigmainv", "samples", "info", "air",
         "loops", "J", "Xbeta_samps", "iis.loglik", "vehtari.loop"),
    envir = environment()
  )
  snow::clusterCall(cl, function(x) {
    for (file in list.files("R")) {
      source(file.path("R", file))
    }

  })

  logliks <- do.call("rbind", snow::parLapply(cl, sample_idx, function(s) {

    a <- Sys.time()
    Sig_inv <- stcarar_sigmainv(
      samples$rho[s, "rho.S"],
      samples$rho[s, "rho.T"],
      info$W, (samples$tau2[s])
    )
    b <- Sys.time()
    Sigi <- Sig_inv[-c(1:4), -c(1:4)]

    shared <- as.vector(difftime(b, a, units = "secs")) / J
    purrr::map_df(unique(loops), function(loop) {
      check_a <- Sys.time()
      test_mask <- loop == loops

      j <- which(test_mask)

      # ghost
      A <- Sigi[j, j]
      B <- Sigi[j, -j]
      D <- Sigi[-j, -j]
      D_inv <- solve(D)
      sig_dd <- D - t(B) %*% solve(A) %*% B
      sig_aa <-  solve(A - B %*% D_inv %*% t(B))
      exp_beta_m <- as.vector(sig_aa %*% B %*% D_inv %*% sig_dd %*%
                                matrix(samples$phi[s, -j],  ncol = 1))

      # iIS
      exp_beta_s <- samples$phi[s, ]
      exp_beta_s[j] <- exp_beta_m
      val_iisc <- iis.loglik(
        test_mask, info$df$observed,
        mu = X %*% samples$beta[s, ] + info$log_offset + exp_beta_s,
        P = NULL, sig = NULL, log = TRUE, glmm = TRUE, dim_theta = length(j),
        X_eff = diag(1, nrow(X)), Sigma_inv = A)
      yhat_iisc <- val_iisc$yhat
      ll_iisc <- val_iisc$loglik
      check_b <- Sys.time()

      # vt 254 seconds for s = 10, loop = loops[1], so excluded
      # val_vt <- vehtari.loop(
      #   mu =  X %*% samples$beta[s,] + info$log_offset + exp_beta_s,
      #   X_r = diag(1, nrow(X)),
      #   reffs = samples$phi[s, ], Y = info$df$observed,
      #   log_offset = info$log_offset, test_mask = test_mask,
      #   Sig_inv = Sigi, family = "poisson")
      # yhat_vt <- val_vt$yhat

      data.frame(
        i = which(test_mask),
        iter = s,
        yhat_iisc = yhat_iisc,
        loglik = ll_iisc,
        time_iisc =  as.vector(difftime(check_b, check_a, units = "secs")),
        loop = loop,
        phi_draw = exp_beta_m
      )
    })
  })) |> dplyr::arrange(i, iter)

  # yhat psiis
  logliks_iis <- matrix(logliks$loglik,
                        nrow = length(sample_idx), ncol = n)
  mu_iis <- matrix(logliks$yhat_iisc,
                   nrow = length(sample_idx), ncol = n)
  r_eff_iis <- loo::relative_eff(
    exp(-logliks_iis),
    chain_id = rep(1:4, each = length(sample_idx) / 4))
  psis_obj_iis <- loo::psis(-logliks_iis, r_eff = r_eff_iis, cores = 2)

  yhat <- loo::E_loo(mu_iis, psis_obj_iis,
                     type = "mean",
                     log_ratios = -logliks_iis
  )
  snow::stopCluster(cl)

  # yhat iis
  # draw posterior predictive samples
  yhat_iis <- logliks |>
    dplyr::group_by(i) |>
    dplyr::summarise(
      yhat_mean = mean(yhat_iisc * exp(-loglik - log(mean(exp(-loglik)))))
    ) |>
    dplyr::arrange(i) |>
    dplyr::pull(yhat_mean)

  df_iis <- data.frame(
    loop = loops,
    yhat_psiis_fm = yhat$value,
    yhat_iis_fm = yhat_iis,
    pk_iis_fm = yhat$pareto_k
  ) |>
    dplyr::bind_cols(
      dplyr::group_by(logliks, i) |>
        dplyr::summarise(time_psiis_fm = sum(time_iisc))
    )

  list(
    P = P,
    muhat = Xbeta + apply(samples$phi, 2, mean),
    rho_s = rho["rho.S"],
    rho_t = rho["rho.T"],
    tau2 = mean(samples$tau2),
    df = data.frame(
      idx = seq_along(Xbeta),
      yhat_post = muhats,
      loop = loops
    ) |> dplyr::full_join(df_iis, by = c("idx" = "i"))
  )
}
