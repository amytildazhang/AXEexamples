# functions to fit methods to compare to

idxs <- function(start = 7, max_iter = 100) {
  M <- rep(NA, max_iter)
  M[1] <- start
  for (iter in 2:max_iter) {
    M[iter] <- (M[iter - 1] * 1.5)
    if (floor(M[iter]) %% 2 == 1) {
      M[iter] <- floor(M[iter])
    } else {
      M[iter] <- ceiling(M[iter])
    }

  }
  M
}

# f: data, tau^2, a_jm, si
ly_sigtau <- function(sig_samps, tau_samps, mu_samps, b_samps, Y, j_idxs,
                      M_start = 7,
                      tol = 0.01, max_iter = 15) {
  mu <- apply(b_samps, 2, mean)
  phi <- apply(b_samps, 2, sd)
  # go 4 standard deviations away

  M <- M_start
  fmj <- matrix(NA, nrow = nrow(b_samps), ncol = ncol(b_samps))
  # fmy <-  matrix(NA, nrow = nrow(b_samps), ncol = length(Y))
  fm_0 <- NA
  for (iter in 1:max_iter) {
    if (iter %% 5 == 0)
      message(sprintf("iter %s", iter))
    gh_rule <- fastGHQuad::gaussHermiteData(M)
    a_jm <- mu +  phi %*% t(gh_rule$x)

    for (j in 1:ncol(b_samps)) {
      j_mask <- 1:length(Y) %in% j_idxs[[j]]
      for (s in 1:nrow(b_samps)) {


        g_m <- dnorm(a_jm[j, ], sd = sig_samps[s], log = T)
        w_js <- log(2*pi)/2 + log(phi[j]) + gh_rule$x^2/2  + g_m + log(gh_rule$w)


        if (!is.null(dim(tau_samps))) {
          tau_s <- tau_samps[s,j_mask]
        } else {
          tau_s <- tau_samps[s]
        }
        fc <- sapply(1:M, function(m) {
          sum(dnorm(
            Y[j_mask],
            mean = mu_samps[s, j_mask] + a_jm[j, m],
            sd = tau_s, log = T))

        }) # M x n_j

        fmj[s, j] <- log(sum(exp(w_js + fc)))

      }

    }
    fm <- apply(fmj, 1, function(sampval) {
      (sum((sampval)))
    })
    max_dif <- max(abs(fm - fm_0))
    if ((is.finite(max_dif) && max_dif <= tol)) break;

    M <- (M * 1.5)
    if (floor(M) %% 2 == 1) {
      M <- floor(M)
    } else {
      M <- ceiling(M)
    }

    fm_0 <- fm
    if (iter == max_iter) message(sprintf("Tolerance not reached, max dif is %s.", max_dif))
  }

  return(list(M = M, fm = (fmj)))
}
# fit posteriors



aq_psis_waic_m <- function(seed, sig_samps, tau_samps, mu_samps, b_samps, Y, j_idxs, yhats, ...) {
  set.seed(seed)
  start_lco <- Sys.time()

  logliks_lco <- ly_sigtau(
    sig_samps = sig_samps, tau_samps = tau_samps,  mu_samps = mu_samps,
    b_samps =b_samps, Y = Y, j_idxs = j_idxs, ...
  )
  M <- logliks_lco$M
  message(sprintf("Adaptive quadrature done with %s nodes", M))
  logliks_lco <- logliks_lco$fm
  infs <- is.infinite(logliks_lco)
  if (any(infs))
    message(sprintf("%s infinite logliks in LCO determined, omitting.", sum(infs)))
  fin_mask <- apply(logliks_lco, 1, function(samp) all(is.finite(samp)))
  r_eff_lco <- loo::relative_eff(exp(-logliks_lco),
                                 chain_id = rep(1:4, each = 1000))
  psis_obj_lco <- loo::psis(-logliks_lco[fin_mask, ], r_eff = r_eff_lco, cores = 2 )


  # yhats
  # draw posterior predictive samples
  psis_yhat_lco <- loo::E_loo(yhats[fin_mask, ], psis_obj_lco, type = 'mean',
                              log_ratios = -logliks_lco[fin_mask, ])

  end_lco <- Sys.time()


  list(
    time = end_lco - start_lco,
    elpd_psis_m = loo::loo(logliks_lco[fin_mask, ])$pointwise[, "elpd_loo"],
    elpd_waic_m = loo::waic(logliks_lco[fin_mask, ])$pointwise[, "elpd_waic"],
    yhat_psis_m = psis_yhat_lco
  )
}

pfit_eight <- function(info = eight, seed = 238972, nodes = 6, ...) {

  # cl <- snow::makeCluster(nodes, outfile = "")
  # snow::clusterExport(cl, list("stanmodels", "info", "seed", "%>%", "ly_sigtau"))

  purrr::map((info$data_scale), function(scl) {
    message(scl)
    scld_dat <- scale_dat(info$stan_list, scl)
    # fit to full data
    fit_sch8 <- rstan::sampling(stanmodels$eightschools_sim, data = scld_dat,
                                seed = seed, refresh = 0, ...)


    samples <- rstan::extract(fit_sch8)
    beta_samps <- cbind(samples$mu, samples$theta)
    Xbeta_samps <- t(info$X %*% t(beta_samps))
    y.pp <- sapply(1:nrow(info$df), function(j) {
      rnorm(n = nrow(Xbeta_samps), mean = Xbeta_samps[, j],
            sd = info$df$sd[j])
    })
    X <- model.matrix(~school, data = info$df,
                      contrasts.arg =
                        list(school = contrasts(info$df$school, contrasts = F)))


    #
    # Integrated Importance sampling (Li, Qiu, Zhang, Feng, 2016, Statistics and Computing)
    #
    iis_og <- iis_fm(1:8, do.call("cbind", lapply(1:8, function(x) samples$mu)),
                     X_b = X[, -1], Y = scld_dat$y,
                     sig_samps = samples$tau^2, P = info$df$sd^2)


    #
    #  Integrated Importance sampling (Vehtari et al, 2016, JMLR)
    #

    iis_vt <- iis_im(1:8, Xbeta_samps, X = X, Y = scld_dat$y,
                     sig_samps =list("school" = samples$tau^2), P = info$df$sd^2,
                     Sigma = "diagonal")


    #
    # PSIS-LOO elpd and yhats (Vehtari, Gelman, Gabry 2017)
    #




    set.seed(seed)
    start <- Sys.time()
    logliks <- sapply(1:nrow(info$df), function(j) {
      dnorm(scld_dat$y[j], mean = Xbeta_samps[, j], sd = info$df$sd[j], log = T)
    })
    r_eff <- loo::relative_eff(exp(-logliks), chain_id = rep(1:4, each = 1000))
    psis_obj <- loo::psis(-logliks,  r_eff = r_eff, cores = 2)
    psis_yhats <- loo::E_loo(Xbeta_samps, psis_obj, type = 'mean',
                             log_ratios = -logliks)
    end <- Sys.time()





    ghst <- ghosting_c(
      seed, 1:8, do.call("cbind", lapply(1:8, function(x) samples$mu)),
      samples$tau^2, X[, -1], scld_dat$y, P = info$df$sd^2
    )

    analy_var <- 9/8 * (info$df$sd^2 + mean(samples$tau)^2)
    #
    # PSIS-LCO elpd and yhats (Marginal log-likelihood as in Merkle et al 2018)
    #
    # lco <- aq_psis_waic_m(seed = seed, sig_samps = samples$tau,
    #                       tau_samps = do.call("rbind", lapply(1:4000, function(x) scld_dat$sigma)),
    #                       mu_samps = do.call("cbind", lapply(1:8, function(x) samples$mu)),
    #                       b_samps = samples$theta,
    #                       Y = scld_dat$y, j_idxs = as.list(1:8), yhats = Xbeta_samps)

    #
    # Integrated importance sampling
    #



    # posterior samples

    data.frame(
      tau_hat = (mean(samples$tau)),
      # tau_map = bayestestR::map_estimate(samples$tau),
      elpd_psis_c =  loo::loo(logliks)$pointwise[, "elpd_loo"],
      elpd_waic_c = loo::waic(logliks)$pointwise[, "elpd_waic"],
      yhat_psis_c = psis_yhats$value,
      data_scale = scl,
      school_idx = 1:8
    ) %>% dplyr::bind_cols(
      iis_og, iis_vt, ghst
    ) %>%
      dplyr::mutate(
        i = 1:dplyr::n(),
        loop = 1:8
      )

  })
  # snow::stopCluster(cl)

  # pfit

}



pfit_radon_full <- function(info = radon_1, seed = 32897, ...) {
  xdat <- info$data %>%
    dplyr::mutate(floor = factor(floor), county = factor(county))

  purrr::map(1:length(info$models_mcv), function(mod_no) {
    message(sprintf("Model %s", mod_no))
    mod <- info$models_mcv[[mod_no]]
    rfit <- rstanarm::stan_glmer(mod, data = info$data, seed = seed,
                                 refresh = 0, ...)
    samples <- rstan::extract(rfit$stanfit)
    if ("beta" %in% names(samples)) {
      coeffs <- cbind(samples$alpha, samples$beta, samples$b)
    } else {
      coeffs <- cbind(samples$alpha, samples$b)
    }
    coeffs <- coeffs[, -ncol(coeffs)]
    X <- stats::model.matrix(info$models_axe[[mod_no]], data = xdat,
                             contrasts = contrasts_for_pooling(xdat))

    loops <- info$data$county
    j_idxes <- lapply(levels(xdat$county), function(lev) which(xdat$county == lev))

    Xbeta_samps <- coeffs %*% t(X)
    fixef <- !stringr::str_detect(colnames(X), "county")
    mu_samps <- coeffs[, fixef] %*% t(X[, fixef])

    y.pp <- sapply(1:nrow(X), function(j) {
      rnorm(n = nrow(Xbeta_samps), mean = Xbeta_samps[, j],
            sd = samples$aux)
    })
    #
    # IIS-C
    #
    message(sprintf("Running IIS-C"))
    iis_og <- iis_fm(loops,
                     mu_samps,
                     X_b = X[, !fixef], Y = info$data$log_radon,
                     sig_samps = samples$theta_L, tau_samps = samples$aux^2, Sigma = "diagonal")



    #
    #  IIS-M
    #
    message("Running IIS-M")
    iis_vt <- iis_im(loops, Xbeta_samps, X = X, Y = info$data$log_radon,
                     sig_samps = list("county" = samples$theta_L),
                     tau_samps = samples$aux^2, n_cores = 4,
                     Sigma = "diagonal")




    message("Ghosting")
    #
    #  Ghosting
    #
    ghst <- ghosting_c(
      seed, loops, mu_samps =  mu_samps,
      sig_samps = samples$theta_L, X_b = X[, !fixef],
      Y =  info$data$log_radon, tau_samps = samples$aux^2, Sigma = "diagonal"
    )

    # posterior samples
    data.frame(
      i = 1:nrow(X),
      tau_hat = (mean(samples$aux)),
      sighat =mean(sqrt(samples$theta_L)),
      tau_map = as.numeric(bayestestR::map_estimate(samples$aux)),
      sig2_map = as.numeric(bayestestR::map_estimate(samples$theta_L)),
      model = mod_no,
      loop = loops
    ) %>%
      dplyr::bind_cols(
        iis_og, iis_vt, ghst
      )

  })
}

pfit_radon_simul <- function(info = radon_2,
                             seed = 9871, n_cores = 5,
                             use_saved = T) { # arguments to pass to rstanarm::stan_g

  loop_over_radon_simul(
    info, seed, n_cores,
    export = list("use_saved", "mm_radon", "coeffs_radon", "radon_simul_samples"),
    expenv = environment(),
    FUN = function(mod_no, dat, n, perc, subset_idx, info, seed, use_saved) {

      sf <-  radon_simul_samples(perc, n, subset_idx, mod_no)
      message(sf)
      if (use_saved) {
        samples <- readr::read_rds(sf)
      } else {
        rfit <- rstanarm::stan_glmer(mod, data = dat, seed = seed, refresh = 0)
        samples <- rstan::extract(rfit$stanfit) #
        readr::write_rds(samples, sf )
      }

      # get fitted values
      X <- mm_radon(mod_no, dat, info$models_axe[[mod_no]])
      coeffs <- coeffs_radon(samples)
      test_mask <- dat$county == info$c_county
      fv <- apply(coeffs %*% t(X), 2, mean)

      # posterior samples
      data.frame(
        yhat_post = as.vector(fv[test_mask]),
        tau2hat_post = mean(samples$aux)^2,
        sig2hat_post = mean(samples$theta_L),
        n_tot = nrow(dat),
        tau2hat_map = as.vector(bayestestR::map_estimate(samples$aux)^2),
        sig2hat_map = as.vector(bayestestR::map_estimate(samples$theta_L))
      )



    }, use_saved = use_saved)


}
#
#
# pfit_fish <- function(info = fish) {
#     xdat <- info$data %>%
#         dplyr::mutate(floor = factor(floor), county = factor(county))
#
#     purrr::map(1:length(info$models_mcv), function(mod_no) {
#         mod <- info$models_mcv[[mod_no]]
#         rfit <- rstanarm::stan_glmer(mod, data = radon_1$data, seed = seed,
#                                      refresh = 0)
#         tau <- rstan::extract(rfit$stanfit, "sigma")$sigma
#         sname <- stringr::str_subset(names(rfit$stanfit), "Sigma\\[")
#         sigma2 <- rstan::extract(rfit$stanfit, sname)[[sname]]
#
#         X <- stats::model.matrix(info$models_axe[[mod_no]], data = xdat,
#                                  contrasts = contrasts_for_pooling(xdat))
#
#         list(tau_hat = mean(tau), sighat = list("county" = mean(sqrt(sigma2))),
#              model = mod_no)
# }


pfit_echo <- function(info = echo) {
  rfit <- rstanarm::stan_glm


  purrr::map(1:length(info$models_mcv), function(mod_no) {
    mod <- info$models_mcv[[mod_no]]
    rfit <- rstanarm::stan_glmer(mod, data = radon_1$data, seed = seed,
                                 refresh = 0)
    tau <- rstan::extract(rfit$stanfit, "sigma")$sigma
    sname <- stringr::str_subset(names(rfit$stanfit), "Sigma\\[")
    sigma2 <- rstan::extract(rfit$stanfit, sname)[[sname]]

    X <- stats::model.matrix(info$models_axe[[mod_no]], data = xdat,
                             contrasts = contrasts_for_pooling(xdat))

    list(tau_hat = mean(tau), sighat = list("county" = mean(sqrt(sigma2))),
         model = mod_no)
  })
}


#
# pfit_tadf <- function(info = tadf, seed = 23982, ...) {
#
#   df <- info$data
#
#   standat <- list(
#     K = length(unique(df$class_attr)),
#     N = nrow(df),
#     D = ncol(X),
#     y = df$class_attr,
#     x = X,
#     J1 = length(unique(df$instructor)),
#     J2 = length(unique(df$course)),
#     b1_ind = as.integer(df$instructor),
#     b2_ind = as.integer(df$course)
#
#   )
#   rfit <- rstan::stan("inst\\stan\\echo_polr.stan", data = standat,
#                       control = list("adapt_delta" = 0.95,
#                                      "max_treedepth" = 15),
#                       cores = 4, refresh = 0, seed = seed, ...)
#
#
# }
#


approx_bin_variance <- function(p, n) {
  p*(1-p)/(n*dnorm(qnorm(p))^2)
}
#
# pfit_tadf_bin <- function(info = tadf, seed = 2875) {
#   pfit <-  rstanarm::stan_glmer(low ~  ne + summer + class_size + (1|instructor) + (1|course),
#                                 family = "binomial", data = info$data,
#                                 refresh = 0, cores = 4, seed = seed)
#   set.seed(2875)
#
#   phats <- rstanarm::posterior_predict(pfit) %>%
#     apply(2, mean)
#   P <- approx_bin_variance(phats, 1)
#
#   samples <- rstan::extract(pfit$stanfit)
#   sigs <- samples$theta_L %>% apply(2, mean)
#
#   coeffs <- cbind(samples$alpha, samples$beta, samples$b)
#   coeffs <- coeffs[, -c(27, 53)] # remove "new" effs
#   Xbeta_samps <- coeffs %*% t(info$X)
#
#   X_effs <- info$X[, stringr::str_detect(colnames(info$X), "instructor")]
#   mu_samps <- Xbeta_samps - samples$b[, -c(1:27, 53)] %*% t(X_effs)
#
#   iis_og <- iis_fm(info$loops,
#                    mu_samps,
#                    X_b = X_effs, Y = qnorm(phats),
#                    sig_samps = samples$theta_L[,2], P = P, Sigma = "diagonal")
#
#
#
#   #
#   #  Integrated Importance sampling (Vehtari et al, 2016, JMLR)
#   #
#
#   iis_vt <- iis_im(info$loops, Xbeta_samps, X = info$X, Y = qnorm(phats),
#                    sig_samps =list("course" = samples$theta_L[,1],
#                                    "instructor" = samples$theta_L[,2]),
#                    P = P,
#                    Sigma = "diagonal")
#   #
#   ghst <- ghosting_c(
#     seed, info$loops, mu_samps =  mu_samps,
#     sig_samps = samples$theta_L[,2], X_b = X_effs,
#     Y =  qnorm(phats), P = P, Sigma = "diagonal"
#   )
#
#
#   list(pfit = pfit,  phats = phats, P = P, sig_course = sigs[1], sig_instr = sigs[2],
#        df = data.frame(
#          idx = 1:length(phats),
#          yhat_post = phats
#        ) %>% dplyr::bind_cols(iis_og, iis_vt, ghst)
#   )
# }
#

pfit_lol <- function(info = lol, seed = 18973, ...) {

  pfit <- rstanarm::stan_glmer(kills ~  position  + team +
                                 log_dpm + log_egpm +
                                 (1| player) + (1|champion),
                               family = "poisson", data = info$data,
                               seed = seed)
  set.seed(seed)
  phats <- rstanarm::posterior_predict(pfit) %>%
    apply(2, mean)
  P <- 1/phats

  Y <- log(info$data$kills + 1e-8)

  samples <- rstan::extract(pfit$stanfit)
  sigs <- samples$theta_L %>% apply(2, mean)

  coeffs <- cbind(samples$alpha, samples$beta, samples$b)
  coeffs <- coeffs[, -c(109, ncol(samples$b))] # remove "new" effs
  Xbeta_samps <- coeffs %*% t(info$X)

  X_effs <- info$X[, stringr::str_detect(colnames(info$X), "player")]
  mu_samps <- Xbeta_samps - samples$b[, -c(1:109, ncol(samples$b))] %*% t(X_effs)

  iis_og <- iis_fm(info$loops,  mu_samps, X_b = X_effs, Y = info$data$kills,
                   sig_samps = samples$theta_L[,2], P = P, Sigma = "car", glmm = T)



  #
  #  Integrated Importance sampling (Vehtari et al, 2016, JMLR)
  #

  iis_vt <- iis_im(info$loops, Xbeta_samps, X = info$X, Y = Y,
                   sig_samps =list("champion" = samples$theta_L[,1],
                                   "player" = samples$theta_L[,2]),
                   P = P,
                   Sigma = "diagonal")

  ghst <- ghosting_c(
    seed, info$loops, mu_samps =  mu_samps,
    sig_samps = samples$theta_L[,2], X_b = X_effs,
    Y =  info$data$kills, P = P, Sigma = "car"
  )


  list(pfit = pfit, phats = phats, P = P, sig_champion = sigs[1], sig_player = sigs[2],

       df = data.frame(
         idx = 1:length(Y),
         yhat_post = phats
       ) %>% dplyr::bind_cols(iis_og, iis_vt, ghst)
  )


}

car_inv <- function(D, W, alpha, tau) {
  tau^2*(D - alpha*W)
}

pfit_slc <- function(info = slc$data, seed = 29873, ...) {

  sdat <- list(n = info$n, p = info$p, X = info$X, y = info$y,
               log_offset = info$log_offset, W_n = info$W_n, W = info$W,
               D_sparse = rowSums(info$W))
  pfit <- rstan::sampling(stanmodels$sparsecar_slp,
                          data = sdat,
                          seed = seed, refresh = 0,
                          iter = 4e3, ...)
  set.seed(seed)
  samples <- rstan::extract(pfit)

  D <- diag(rowSums(info$W))
  Xbeta <- info$X %*% apply(samples$beta, 2, mean)

  phats <- Xbeta + apply(samples$phi, 2, mean)
  muhats <- exp(phats + info$log_offset)


  P <- as.vector(1/muhats)


  Y <- log(info$y + 1e-4)


  X <- cbind(info$X, diag(rep(1, length(Y))))
  coeffs <- cbind(samples$beta, samples$phi)
  Xbeta_samps <- coeffs %*% t(X)



  a <- Sys.time()
  spat_mu <- sapply(1:length(Y), function(i) {
    Nbhd_i <- which(info$W[i, ] != 0)

    samples$alpha*rowSums(samples$phi[, Nbhd_i, drop = F])/length(Nbhd_i)
  })



  mu_samps <- apply(samples$beta %*% t(info$X) + spat_mu, 1, function(samp) {
    samp + info$log_offset
  }) %>% t()
  b <- Sys.time()

  # iis_og <- iis_fm(
  #   info$loops, mu_samps, X_b = diag(rep(1, length(Y))),
  #   Y = Y, sig_samps = 1/samples$tau, P = P, Sigma = "car"
  # ) %>%
  #   dplyr::mutate_at(dplyr::vars(dplyr::starts_with("yhat_")), exp)


  iis_fm_mc <- iis_fm(
    info$loops, mu_samps, X_b = diag(rep(1, length(Y))),
    Y = info$y, sig_samps = 1/samples$tau, P = P, Sigma = "car", glmm = T
  )


  #
  #  Integrated Importance sampling (Vehtari et al, 2016, JMLR)
  #

  iis_vt <- iis_im(info$loops, Xbeta_samps, X = X, Y = phats,
                   sig_samps =list("alpha" = samples$alpha,
                                   "tau" = samples$tau^2),
                   P = P, Sigma = "car",
                   Sigma_objs = list(D = diag(rowSums(info$W)),
                                     W = info$W),
                   downsample = T) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::starts_with("yhat_")), exp)

  ghst <- ghosting_c(
    seed, info$loops, mu_samps =  mu_samps,
    sig_samps = 1/samples$tau^2, X_b = diag(rep(1, length(Y))),
    Y =  info$y, P = P, Sigma = "car", M = M
  ) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::starts_with("yhat_")), exp)



  list(
    P = P,
    tauhat_post = mean(samples$tau),
    alpha_post = mean(samples$alpha),
    muhat = Xbeta + apply(samples$phi, 2, mean),
    mutime = difftime(b, a, units = "secs"),
    df = data.frame(
      idx = 1:length(Y),
      yhat_post = muhats
    ) %>% dplyr::bind_cols(iis_fm_mc, iis_vt, ghst)

  )


}



stcarar_sigmainv <- function(rho_s, rho_t, W, tau2, P = 4 + N,
                             N = nrow(air$data$df), J = N/Ts, Ts=5) {
  sigmai_spat <- (rho_s*(diag(rowSums(W)) - W) +
                    diag(1 - rho_s, nrow = nrow(W)))/tau2
  sigmai_blck <- matrix(
    sapply(1:Ts, function(t) {
      sapply(1:ncol(sigmai_spat), function(j) {
        c(rep(0, J*(t - 1)), sigmai_spat[, j], rep(0, J*(Ts - t)))
      })
    }),
    ncol = N
  )
  Sigma_inv <- matrix(0, ncol = P, nrow = P)
  sigma_idxs <- (P - N + 1):P

  H <- matrix(0, N, N)
  H[(J + 1):N, 1:(N - J)] <- rho_t*diag(rep(1, N - J))
  I <- diag(rep(1, N))


  Sigma_inv[sigma_idxs, sigma_idxs] <- t(I - H) %*% sigmai_blck %*% (I - H)

  Sigma_inv

}

rMVNormP <- function(n, mu, Omega){
  p <- length(mu)
  Z <- matrix(rnorm(p*n), p, n)
  U <- chol(Omega) # By default R's chol fxn returns upper cholesky factor
  X <- backsolve(U, Z) # more efficient and stable than acctually inverting
  X <- sweep(X, 1, mu, FUN=`+`)
  return(X)
}

rMVNormC <- function(n, mu, Sigma){
  p <- length(mu)
  Z <- matrix(rnorm(p*n), p, n)
  L <- t(chol(Sigma)) # By default R's chol fxn returns upper cholesky factor
  X <- L%*%Z
  X <- sweep(X, 1, mu, FUN=`+`)
  return(X)
}


pfit_air <- function(info = air$data, seed = 69872, ...) {

  set.seed(seed)
  formula <- info$formula

  chain1 <- CARBayesST::ST.CARar(formula = formula, family = "poisson",
                                 data = info$df, W = info$W, burnin = 20000, n.sample = 220000,
                                 thin = 200)
  chain2 <- CARBayesST::ST.CARar(formula = formula, family = "poisson",
                                 data = info$df, W = info$W, burnin = 20000, n.sample = 220000,
                                 thin = 200)
  chain3 <- CARBayesST::ST.CARar(formula = formula, family = "poisson",
                                 data = info$df, W = info$W, burnin = 20000, n.sample = 220000,
                                 thin = 200)
  chain4 <- CARBayesST::ST.CARar(formula = formula, family = "poisson",
                                 data = info$df, W = info$W, burnin = 20000, n.sample = 220000,
                                 thin = 200)

  samples <- purrr::map(names(chain1$samples), function(param) {
    rbind(
      chain1$samples[[param]],
      chain2$samples[[param]],
      chain3$samples[[param]],
      chain4$samples[[param]]
    )

  }) %>% purrr::set_names(names(chain1$samples))

  # samples <- rstan::extract(pfit)


  D <- diag(rowSums(info$W))
  Xbeta <- info$X %*% apply(samples$beta, 2, mean)

  phats <- Xbeta + apply(samples$phi, 2, mean)
  muhats <- exp(phats + info$log_offset)



  rho <- apply(samples$rho, 2, mean)
  tau2 <- mean(samples$tau2)
  P <- as.vector(1/muhats)
  Y <- log(info$df$observed + 1e-5) - info$log_offset
  loops <- info$loops

  # axe_post <- X %*% V %*% t(X) %*% diag(1/P) %*% Y

  Xbeta_samps <- samples$beta %*% t(info$X)

  sample_idx <- sample(1:nrow(Xbeta_samps), 1000)
  cl <- snow::makeCluster(6)
  snow::clusterExport(
    cl,
    list("stcarar_sigmainv", "samples", "info", "air", "loops", "%>%", "Xbeta_samps"),
    envir = environment())
  start_iis <- Sys.time()
  {


    logliks <- do.call("rbind", snow::parLapply(cl, sample_idx, function(s) {

      a <- Sys.time()
      Sig_inv <-  stcarar_sigmainv(samples$rho[s, "rho.S"],
                                   samples$rho[s, "rho.T"],
                                   info$W, (samples$tau2[s]))
      Sigi <- Sig_inv[-c(1:4), -c(1:4)]
      Sig <- solve(Sigi)

      purrr::map_df(unique(loops), function(loop) {
        test_mask <- loop == loops


        draws <- condMVNorm::rcmvnorm(
          n = 201, mean = rep(0, nrow(Sig)), Sig,
          which(test_mask), which(!test_mask),
          samples$phi[s, !test_mask], check.sigma = F,
          method = c("chol")
        )



        loglik <- sapply(1:sum(test_mask), function(idx) {
          test_idx <- which(test_mask)[idx]

          log(mean(dpois(
            info$df$observed[test_idx],
            exp(draws[-201, idx]+ Xbeta_samps[s, test_idx] + info$log_offset[test_idx])
          )))


        })

        phi_test <- apply(draws, 2, mean)

        data.frame(
          i = which(test_mask),
          iter = s,
          yhat = exp(phi_test + Xbeta_samps[s, test_mask] + info$log_offset[test_mask]),
          loglik = sum(loglik),
          loop = loop,
          phi_draw = as.vector(draws[201])
        )

      })


    })) %>% dplyr::arrange(i, iter)

    snow::stopCluster(cl)

    # yhat psiis
    logliks_iis <- matrix(logliks$loglik, nrow = length(sample_idx), ncol = length(Y))
    mu_iis <- matrix(logliks$yhat, nrow = length(sample_idx), ncol = length(Y))
    r_eff_iis <- loo::relative_eff(exp(-logliks_iis), chain_id = rep(1:4, each = length(sample_idx)/4))
    psis_obj_iis <- loo::psis(-logliks_iis, r_eff = r_eff_iis, cores = 2 )

    yhat <- loo::E_loo(mu_iis, psis_obj_iis, type = 'mean',
                       log_ratios = -logliks_iis)

  }
  end_iis <- Sys.time()

  difftime(end_iis, start_iis, units = "hours")
  # yhat iis
  # draw posterior predictive samples
  yhat_iis <- logliks %>%
    dplyr::group_by(i) %>%
    dplyr::summarise(yhat_mean = mean(yhat * exp(-loglik - log(mean(exp(-loglik)))))) %>%
    dplyr::pull(yhat_mean)


    yhat_ghst <- (logliks %>% dplyr::group_by(i) %>%
      dplyr::summarise(yhat = mean(phi_draw)) %>% dplyr::pull(yhat)) +
      apply(Xbeta_samps, 2, mean) + info$log_offset



  df_means <- data.frame(
    loop = loops,
    yhat_psiis_fm = yhat$value,
    yhat_iis_fm = yhat_iis,
    pk_iis_fm = yhat$pareto_k,
    yhat_ghst_c = exp(yhat_ghst)

  )


  # dic


  Sig_inv <- stcarar_sigmainv(mean(samples$rho[, "rho.S"]),
                              mean(samples$rho[, "rho.T"]),
                              info$W, mean(samples$tau2))[-c(1:4), -c(1:4)]
  Sig <- solve(Sig_inv)
  phi <- apply(samples$phi, 2, mean)
  Xbeta <- apply(Xbeta_samps, 2, mean)
  dic <- -(logliks %>%
             dplyr::select(loop, loglik, iter) %>% unique() %>%
             dplyr::group_by(loop) %>%
             dplyr::summarise(
               e_loglik = mean(loglik)
             ) %>%
             dplyr::mutate(
               loglik_atmean = purrr::map_dbl(
                 loop, function(lp) {
                   test_mask <- lp == loops

                   draws <- condMVNorm::rcmvnorm(
                     n = 200, mean = rep(0, nrow(Sig)), Sig,
                     which(test_mask), which(!test_mask),
                     phi[!test_mask], check.sigma = F,
                     method = c("chol")
                   )

                   loglik <- sapply(1:sum(test_mask), function(idx) {
                     test_idx <- which(test_mask)[idx]

                     log(mean(dpois(
                       info$df$observed[test_idx],
                       exp(draws[, idx]+ Xbeta[test_idx] + info$log_offset[test_idx])
                     )))


                   })


                   sum(loglik)
                 }),
               p_d = -2*e_loglik + 2*loglik_atmean,
               dic = -2*loglik_atmean + 2*p_d
             ) %>%
             dplyr::pull(dic))/2


  # elpd psiis
  logliks_byg <- logliks %>% dplyr::select(loop, iter, loglik) %>%
    unique() %>% dplyr::pull(loglik) %>%
    matrix(nrow = length(sample_idx), ncol = length(unique(loops)))


  elpd_psiis_fm <-  loo::loo(
    logliks_byg,
    r_eff = loo::relative_eff(exp(logliks_byg), chain_id = rep(1:4, each = length(sample_idx)/4))
  )$pointwise[, "elpd_loo"]

  waic_fm <- loo::waic(logliks_byg)


  time <- as.double(difftime(end_iis, start_iis, units = "secs"))*6/length(unique(loops))
  df_elpds <- data.frame(
    loop = (unique(loops)),
    time_iis_fm = (time),
    elpd_psiis_fm = elpd_psiis_fm,
    elpd_waic_fm = waic_fm$pointwise[, "elpd_waic"],
    elpd_dic_fm = dic


  )

  iis_fm_mc <- dplyr::full_join(df_means, df_elpds, by = c())
  # waic

  #
  #  Integrated Importance sampling (Vehtari et al, 2016, JMLR)
  #

  # ghosting


  end <- Sys.time()
  time_yhat <- as.double(difftime(end, start, units = 'secs'))

  elpds <- purrr::map_dbl(1:ncol(yhats), function(n) {

    log(mean(dpois(Y[n], exp(yhats[, n]))))
  })
  end <- Sys.time()


  ghst <-  data.frame(
    time_ghst_c = (time_yhat)/J,
    elpd_ghst_c = elpds,
    yhat_ghst_c = apply(yhats, 2, mean),
    var_ghst = apply(yhats, 2, var)
  ) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::starts_with("yhat_")), exp)



  list(
    P = P,
    muhat = Xbeta + apply(samples$phi, 2, mean),
    rho_s = rho["rho.S"],
    rho_t = rho["rho.T"],
    tau2 = mean(samples$tau2),
    df = data.frame(
      idx = 1:length(Y),
      yhat_post = muhats,
      loop = loops
    ) %>% dplyr::full_join(dplyr::full_join(df_elpds, df_means %>%
                                              dplyr::mutate(idx = 1:length(Y))))

  )

}
