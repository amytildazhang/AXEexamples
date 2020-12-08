# Contains LCO method functions for all examples except SRD

# convenience function
match_to_x <- function(signames, X) {
    sapply(colnames(X), function(coefnm) {
        matches <- stringr::str_detect(coefnm, signames)
        if (any(matches)) {
            stringr::str_which(coefnm, signames)
        } else {
            NA
        }
    })
}

#
# iIS-A
#

# calculations done within each loop for iIS-A
iis_im.loop <- function(loop, loops, X, P, Sig_inv, Y, V = NULL) {
    test_mask = loops == loop
    X_t <- X[!test_mask, ]
    P_t <- P[!test_mask]
    X_test <- X[test_mask, , drop = F]
    P_test <- P[test_mask]

    n_test <- sum(test_mask)
    if (n_test < ncol(X) & !is.null(V)) {
        if (n_test == 1) {
            V_t <- V + V %*% t(X_test) %*% X_test %*% V/(P_test - as.vector(X_test %*% V %*% t(X_test)))
        } else {
            V_t <- V + V %*% t(X_test) %*% solve(diag(P_test) - X_test %*% V %*% t(X_test)) %*% X_test %*% V
        }
        V_t <- round(V_t, digits = 12)
    } else {
        V_t <- solve(t(X_t) %*% diag(1/P_t) %*% X_t + Sig_inv)
    }


    XVX <- X_test %*%
        V_t %*% t(X_t) %*% diag(1/P_t)

    axe_yhat <- XVX %*% Y[!test_mask]
    cov <- diag(P_test, nrow = length(P_test)) + (X_test %*% V_t %*% t(X_test))

    if (length(cov) > 1) {
      prec <- diag(1/P_test) - diag(1/P_test) %*% X_test %*% V %*%
        t(X_test) %*% diag(1/P_test)
    } else {
      prec <- 1/P_test - X_test %*% V %*% t(X_test)/P_test^2
    }


    dif <- Y[test_mask] - axe_yhat
    ll <- as.vector(-n_test /2 * log(2*pi) - log(det(cov))/2 - t(dif) %*% prec %*% dif/2)


    # # ll <- dnorm(Y[test_mask], mean = axe_yhat, sd = sd , log = T)
    # ll <- sum(dnorm(Y[test_mask], mean = axe_yhat, sd = sd , log = T))

    data.frame(loop = loop,
               i = which(test_mask),
               loglik = (ll),
               yhat = axe_yhat
    )
}
# function which takes samples, performs calculations a la iis_im.loop to get
# log ratios, then inputs them to loo::loo to obtain yhats, DIC, WAIC
iis_im <- function(loops, mu_samps, X, Y, sig_samps, tau_samps = NULL, P = NULL,
                  loop = NULL, n_cores = 0, Sigma = c("diagonal", "car"),
                  Sigma_objs = list(), downsample = T, family = "normal") {
    # all samples are provided as matrices with mcmc iter X datapoint
    # must have only one of tau or P null
    # all values are provided as variances

    # helpful mask to create Sigma^-1
    if (Sigma == "diagonal") {
      xmatch <- match_to_x(names(sig_samps), X)
      reff_mask <- !is.na(xmatch)
    }

    S <- length(sig_samps[[1]]) # number of MCMC iterations
    N <- nrow(X) # number of data points
    if (n_cores > 0) {
        cl <- snow::makeCluster(n_cores, outfile = '')
        snow::clusterExport(cl, list("iis_im.loop", "loops", "X", "P", "Y", "%>%"),
                            envir = environment())
    }


    if (downsample) {
      S <- sample(1:S, 1000) # to cut down on time...

    } else {
      S <- 1:nrow(mu_samps)
    }

    start_iis <- Sys.time()
    logliks_y_cp <- purrr::map_df(S, function(s) { # for each MC sample
        if (is.null(P)) {
            P <- rep(tau_samps[s], N)

        }
        # create Sigma^-1
        if (Sigma == "diagonal") {
          sigi_diag <- rep(0, ncol(X))
          sigs <- purrr::map_dbl(sig_samps, ~1/.[s])
          sigi_diag[reff_mask] <- sigs[xmatch[reff_mask]]
          Sig_inv <- diag(sigi_diag)
        } else if (Sigma == "car") {
          prec <- (sig_samps$tau[s])*(Sigma_objs$D - sig_samps$alpha[s]*Sigma_objs$W)
          fixed <- ncol(X) - ncol(prec)
          Sig_inv <- matrix(0, ncol(X), ncol(X))
          Sig_inv[-c(1:fixed), -c(1:fixed)] <-prec
        }

        # solve for posterior V
        V <- solve(t(X) %*% diag(1/P) %*% X + Sig_inv)
        if (is.null(loop)) {
            if (n_cores > 0) {

                df <- do.call(dplyr::bind_rows, snow::parLapply(cl, unique(loops), function(loop, V, Sig_inv) {
                    iis_im.loop(loop, loops, X, P, Sig_inv, Y, V)
                }, V, Sig_inv))
            } else {
                df <- purrr::map_df(unique(loops), function(loop) {

                    iis_im.loop(loop, loops, X, P, Sig_inv, Y, V)


                })

            }
        } else {

            df <- iis_im.loop(loop, loops, X, P, Sig_inv, Y, V)
        }

        df %>% dplyr::mutate(iter = s)

    }) %>% dplyr::arrange(i, iter)

    print(colnames(logliks_y_cp))

    if (n_cores > 0) {
        snow::stopCluster(cl)
    }
    end_iis <- Sys.time()

    logliks_iis <- matrix(logliks_y_cp$loglik, nrow = length(S))



    r_eff_iis <- loo::relative_eff(exp(-logliks_iis), chain_id = rep(1, length(S)))
    psis_obj_iis <- loo::psis(-logliks_iis, r_eff = r_eff_iis, cores = 2 )


    # yhats
    # draw posterior predictive samples
    # elpd <- loo::loo(logliks_iis,
    #                  r_eff = loo::relative_eff(exp(logliks_iis), chain_id = rep(1, length(S))))
    yhat_iis1 <- logliks_y_cp %>%
        dplyr::group_by(i) %>%
        dplyr::summarise(yhat_mean = mean(yhat * exp(-loglik -
                                                         log(mean(exp(-loglik)))))) %>%
        dplyr::pull(yhat_mean)

    ys_iis <- logliks_y_cp %>% dplyr::pull(yhat) %>%  matrix(nrow = length(S))

    yhat_iis2 <- loo::E_loo(ys_iis, psis_obj_iis, type = 'mean',
                            log_ratios = -logliks_iis)

    last_iis <- Sys.time()

    sighat <- loo::E_loo(do.call("cbind", purrr::map(1:ncol(logliks_iis), ~sig_samps[[1]][S])),
                         psis_obj_iis, type = 'mean',
                         log_ratios = -logliks_iis)
    logliks_byg <- logliks_y_cp %>% dplyr::select(loop, iter, loglik) %>%
        unique() %>% dplyr::pull(loglik) %>%
        matrix(nrow = nrow(mu_samps), ncol = length(unique(loops)))
    if (!is.null(loop)) {
        logliks_byg <- logliks_byg[, unique(loops) == loop, drop = F]
    }


    waic_im <- loo::waic(logliks_byg)

    sighats <- lapply(sig_samps, mean) %>% purrr::set_names(names(sig_samps))
    if (is.null(P)) {
      P <- rep(mean(tau_samps), length(Y))
    }
    if (Sigma == "diagonal") {
      sigi_diag <- rep(0, ncol(X))
      sigs <- purrr::map_dbl(sighats, ~1/.)
      sigi_diag[reff_mask] <- sigs[xmatch[reff_mask]]
      Sig_inv <- diag(sigi_diag)
    } else if (Sigma == "car") {
      prec <- (sighats$tau)^2*(Sigma_objs$D - sighats$alpha*Sigma_objs$W)
      fixed <- ncol(X) - ncol(prec)
      Sig_inv <- matrix(0, ncol(X), ncol(X))
      Sig_inv[-c(1:fixed), -c(1:fixed)] <- prec
    }
    V <- solve(t(X) %*% diag(1/P) %*% X + Sig_inv)

    elpd_dic <- -dic_im(loops, Sig_inv, X, Y, logliks_y_cp, P = P, V = V)/2


    if (!is.null(loop)) {
        loops <- loop
    }
    J <- length(unique(logliks_y_cp$loop))
    time_per_loop <- as.double(difftime(end_iis, start_iis, units = "secs"))*max(n_cores, 1)/J

    data.frame(
        loop = loops,
        time_iis_im = time_per_loop,
        elpd_iis_im = loo::loo(logliks_iis)$pointwise[, "elpd_loo"],
        yhat_psiis_im = yhat_iis2$value,
        yhat_iis_im = yhat_iis1,
        sighat_iis_im = sighat$value
    ) %>%
        dplyr::full_join(
            data.frame(
                loop = unique(loops),
                elpd_waic_im = waic_im$pointwise[, "elpd_waic"],
                elpd_dic_im = elpd_dic,
                impvar_iis_im = apply(logliks_byg, 2, sd)

            ), by = c()
        ) %>%
        dplyr::select(-loop)

}



dic_im.loglik_at_hats <- function(loop, loops, X, Y, P, siginvhat, V) {
    as.vector(unique(iis_im.loop(loop, loops, X, P, siginvhat, Y, V)$loglik))

}

dic_im <- function(loops, siginvhat, X, Y, logliks, P = NULL, tauhat = NULL, V) {
    # calculate DIC
    if (is.null(P)) {
        P <- rep(tauhat, length(Y))
    }

    logliks %>%
        dplyr::select(loop, loglik, iter) %>% unique() %>%
        dplyr::group_by(loop) %>%
        dplyr::summarise(
            e_loglik = mean(loglik)
        ) %>%
        dplyr::mutate(
            loglik_atmean = purrr::map_dbl(loop, ~dic_im.loglik_at_hats(., loops, X, Y, P, siginvhat, V)),
            p_d = -2*e_loglik + 2*loglik_atmean,
            dic = -2*loglik_atmean + 2*p_d
        ) %>%
        dplyr::pull(dic)
}



#
#  iIS-C
#

# calculations done within each loop for iIS-C
iis_fm.loglik <- function(Y_i, mu, P_i, sig, log = T,
                          Sigma = c("diagonal", "car"), Sigma_objs = list(), glmm = F) {

    n_i <- length(Y_i)
    dif <- Y_i - mu
    if (Sigma == "diagonal") {
      if (n_i == 1) {
        dnorm(Y_i, mean = mu, sd = sqrt(unique(P_i) + sig), log = log)
      } else {
        ones <- rep(1, n_i)

        - n_i /2 * log(2*pi) -  sum(P_i)/2 - log(sig)/2 - log((1/sig + sum(1/P_i)))/2 -
          t(dif) %*% solve(diag(P_i) + sig^2 * ones %*% t(ones)) %*% dif/2

      }
    } else if (Sigma == "car" & glmm) {
      sum(sapply(1:length(Y_i), function(idx) {
        si_samps <- rnorm(200, mu[idx], sd = sqrt(sig))
        log(mean(dpois((Y_i[idx]), exp(si_samps))))


      }))
    } else if (Sigma == "car") {
      sum(dnorm(Y_i, mean = mu, sd = sqrt(P_i + sig), log = T))
    } else if (Sigma == "st_car") {

    }
}
# function which takes samples, performs calculations a la iis_fm.loop to get
# log ratios, then inputs them to loo::loo to obtain yhats, DIC, WAIC
iis_fm <- function(loops, mu_samps,  X_b, Y, sig_samps, tau_samps = NULL, P = NULL,
                   Sigma = c("diagonal", "car"), glmm = F) {
    # X is without the b_i
    # b_samps is an MCMC x data point matrix, with the b_i for each data point within
    # all values are provided as variances
    # either tau_samps or P must be null

    N <- length(Y)
    start_iis <- Sys.time()

    logliks <- purrr::map_df(1:nrow(mu_samps), function(s) {
        if (is.null(P)) {
            P_t <- rep(tau_samps[s], N)

        } else {
            P_t <- P
        }
        purrr::map_df(unique(loops), function(loop) {
          test_mask <- loop == loops

          loglik <- iis_fm.loglik(
                Y[test_mask], mu_samps[s, test_mask],
                (P_t[test_mask]), sig_samps[s],
                Sigma = Sigma,
                glmm = glmm)

          if (glmm) {
            yhat <- exp(mu_samps[s, test_mask])
          } else {
            yhat <- mu_samps[s, test_mask]
          }

            data.frame(
                i = which(test_mask),
                iter = s,
                yhat = yhat,
                loglik = as.vector(loglik),
                loop = loop
            )

        })


    }) %>% dplyr::arrange(i, iter)


    # yhat psiis
    logliks_iis <- matrix(logliks$loglik, nrow = nrow(mu_samps), ncol = ncol(mu_samps))
    mu_iis <- matrix(logliks$yhat, nrow = nrow(mu_samps), ncol = ncol(mu_samps))
    r_eff_iis <- loo::relative_eff(exp(-logliks_iis), chain_id = rep(1:4, each = nrow(mu_samps)/4))
    psis_obj_iis <- loo::psis(-logliks_iis, r_eff = r_eff_iis, cores = 2 )

    yhat <- loo::E_loo(mu_iis, psis_obj_iis, type = 'mean',
                       log_ratios = -logliks_iis)
    end_iis <- Sys.time()

    # yhat iis
    # draw posterior predictive samples
    yhat_iis <- logliks %>%
        dplyr::group_by(i) %>%
        dplyr::summarise(yhat_mean = mean(yhat * exp(-loglik - log(mean(exp(-loglik)))))) %>%
        dplyr::pull(yhat_mean)




    # sighat psiis
    sighat <- loo::E_loo(do.call("cbind", purrr::map(1:ncol(logliks_iis), ~sig_samps)),
                         psis_obj_iis, type = 'mean',
                         log_ratios = -logliks_iis)
    if (!is.null(tau_samps)) {
      tauhat <- loo::E_loo(
        do.call("cbind", purrr::map(1:ncol(logliks_iis), ~tau_samps)),
        psis_obj_iis, type = 'mean',
        log_ratios = -logliks_iis
      )
    } else {
      tauhat <- list(value = NA)
    }


    df_means <- data.frame(
      loop = loops,
      yhat_psiis_fm = yhat$value,
      yhat_iis_fm = yhat_iis,
      pk_iis_fm = sighat$pareto_k,
      sig2hat_iis_fm = sighat$value,
      tau2hat_iis_fm = tauhat$value
    )


    # dic

    dic <- -dic_fm(
      loops, apply(mu_samps, 2, mean), mean(sig_samps),
      Y, logliks, P = P, tauhat = mean(tau_samps),
      Sigma = Sigma, glmm = glmm
    )/2



    # elpd psiis
    logliks_byg <- logliks %>% dplyr::select(loop, iter, loglik) %>%
      unique() %>% dplyr::pull(loglik) %>%
      matrix(nrow = nrow(mu_samps), ncol = length(unique(loops)))


    elpd_psiis_fm <-  loo::loo(
      logliks_byg,
      r_eff = loo::relative_eff(exp(logliks_byg), chain_id = rep(1:4, each = nrow(mu_samps)/4))
    )$pointwise[, "elpd_loo"]

    waic_fm <- loo::waic(logliks_byg)


    time <- as.double(difftime(end_iis, start_iis, units = "secs"))/length(unique(loops))
    df_elpds <- data.frame(
      loop = (unique(loops)),
      time_iis_fm = (time),
      elpd_psiis_fm = elpd_psiis_fm,
      elpd_waic_fm = waic_fm$pointwise[, "elpd_waic"],
      elpd_dic_fm = dic,
      impvar_iis_fm = apply(logliks_byg, 2, sd)


    )
    # waic


    dplyr::full_join(
      df_means,
      df_elpds, by = c()
    ) %>% dplyr::select(-loop)

}



dic_fm.loglik_at_hats <- function(loop, loops, Y, muhat, P, sighat,
                                  Sigma = c("car", "diagonal"), Sigma_objs = list(),
                                  glmm = F) {
        test_mask <- loop == loops
        as.vector(iis_fm.loglik(
          Y[test_mask], muhat[test_mask],
          (P[test_mask]), sighat,  Sigma = Sigma, glmm))

}

dic_fm <- function(loops, muhat, sighat,Y, logliks, P = NULL, tauhat = NULL,
                   Sigma = c("diagonal", "car"), glmm = F) {
    # calculate DIC
    if (is.null(P)) {
        P <- rep(tauhat, length(Y))
    }

   logliks %>%
        dplyr::select(loop, loglik, iter) %>% unique() %>%
        dplyr::group_by(loop) %>%
        dplyr::summarise(
            e_loglik = mean(loglik)
        ) %>%
        dplyr::mutate(
            loglik_atmean = purrr::map_dbl(
              loop,
              ~dic_fm.loglik_at_hats(., loops, Y, muhat, P, sighat, Sigma = Sigma, glmm = glmm)),
            p_d = -2*e_loglik + 2*loglik_atmean,
            dic = -2*loglik_atmean + 2*p_d
        ) %>%
       dplyr::pull(dic)
  }


#
#  GHOST
#

ghosting_c <- function(seed, loops, mu_samps, sig_samps,
                       X_b, Y, tau_samps = NULL, P = NULL,
                       Sigma = c("diagonal", "car"), M = NULL) {
    set.seed(seed)
    start <- Sys.time()

    J <- length(unique(loops))

      b_samps <- do.call("rbind", purrr::map(1:length(sig_samps), function(s) {
        rnorm(J, sd = sqrt(sig_samps[s]))
      }) )

    # S x N + S x J x N
    yhats <- mu_samps + b_samps %*% t(X_b)



    end <- Sys.time()
    time_yhat <- as.double(difftime(end, start, units = 'secs'))

    elpds <- purrr::map_dbl(1:nrow(X_b), function(n) {
        if (is.null(P)) {
            sd = sqrt(tau_samps)
        } else {
            sd = sqrt(P[n])
        }

      if (Sigma == "diagonal") {
        log(mean(dnorm(Y[n], yhats[, n], sd = sd)))

      } else if (Sigma == "car") {
        log(mean(dpois(Y[n], exp(yhats[, n]))))
      }
    })
    end <- Sys.time()


    data.frame(
        time_ghst_c = (time_yhat)/J,
        elpd_ghst_c = elpds,
        yhat_ghst_c = apply(yhats, 2, mean),
        var_ghst = apply(yhats, 2, var)
    )

}



#
# Radon subsets -- load saved posterior samples and then obtain approximations
#
lco_radon_simul <- function(info = radon_2, seed = 8971, n_cores = 5,
                            method = c("iis_fm", "ghost", "iis_im")) {
    loop_over_radon_simul(
        info, seed, n_cores,
        export = list("mm_radon", "coeffs_radon", "radon_simul_samples",
                      "iis_fm", "ghosting_c", "iis_im", "iis_fm.loglik",
                      "iis_im.loop", "match_to_x", "dic_fm", "dic_im",
                      "dic_fm.loglik_at_hats", "dic_im.loglik_at_hats"),
        expenv = environment(),
        FUN = function(mod_no, dat, n, perc, subset_idx, info, seed, method) {

            sf <-  radon_simul_samples(perc, n, subset_idx, mod_no)
            samples <- readr::read_rds(sf)
            coeffs <- coeffs_radon(samples)

            X <- mm_radon(mod_no, dat, info$models_axe[[mod_no]])

            test_mask <- dat$county == info$c_county
            cv_loop <- rep(1, sum(test_mask))


            X_test <- X[test_mask, , drop = F]
            Y_test <- dat$log_radon[test_mask]
            Xbeta_samps <- coeffs %*% t(X_test)

            if (method %in% c("iis_fm", "ghost")) {

                X_effs <- X_test[, stringr::str_detect(colnames(X), "county")]
                mu_samps <- Xbeta_samps - samples$b[, -c(n + 1)] %*% t(X_effs)

                X_b <- X_test[, sprintf("county%s", info$c_county), drop = F]

                if (method== "iis_fm") {
                    iis_fm(cv_loop,
                           mu_samps,
                           X_b = X_b, Y = Y_test,
                           sig_samps = samples$theta_L, tau_samps = samples$aux^2)
                } else {
                    ghosting_c(
                        seed, cv_loop, mu_samps =  mu_samps,
                        sig_samps = samples$theta_L, X_b = X_b,
                        Y =  info$data$log_radon[test_mask], tau_samps = samples$aux^2
                    )
                }
            } else if (method == "iis_im") {

                iis_im(loops = dat$county %in% info$c_county,
                       Xbeta_samps, X = X, Y = dat$log_radon,
                       sig_samps = list("county" = samples$theta_L),
                       tau_samps = samples$aux^2, loop = TRUE,
                       Sigma = "diagonal")

            }
        }, method = method)

}




