##### helper functions for MCV


.collate_folds <- function(objs) {
    obj_names <- names(objs[[1]])

    loops <- purrr::map_lgl(objs, ~length(.$idx) > 1)
    if (!all(loops))  {
        i <- 1
    }  else {
        i <- which(loops)[1]
    }

    purrr::map(obj_names, function(obnm) {
        ex <- objs[[i]][[obnm]]

        if (checkmate::test_atomic_vector(ex)) {
            op <- "c"
        } else {
            op <- "rbind"
        }

        bind <- do.call(op, purrr::map(objs, ~.[[obnm]]))
        names(bind) <- NULL

        bind
    }) %>% purrr::set_names(obj_names)


}

###### functions within each fold

crossval_eight <-  function(loop, data, loops, X, seed, ...) {
    a <- Sys.time()
    test_mask <- loop == loops
    cv_dat <- data

    cv_dat$train_idx <- which(!test_mask)
    cv_dat$test_idx <- as.array(which(test_mask))
    cv_dat$H = sum(test_mask)

    fit_cv <- rstan::stan("inst/stan/eightschools_sim.stan", data = cv_dat,
                              refresh = 0, seed = seed, ...)

    cv_betas <- t(do.call("cbind",
                          rstan::extract(fit_cv, pars = c("mu", "theta"))))
    cv_samps <- X %*% cv_betas
    cv_taus <- rstan::extract(fit_cv, pars = "tau")$tau

    cv_elppd <- -log(mean(dnorm(data$y[test_mask],
                                mean = cv_samps[test_mask, ],
                                sd = data$sigma[test_mask])))
    P <- data$sigma^2
    X_t <- X[!test_mask, ]
    P_t <- P[!test_mask]
    logliks_y <- lapply(1:length(cv_taus), function(s) {
        sigi_diag <- rep(0, ncol(X))
        sigi_diag[2:9] <- 1/cv_taus[s]^2
        Sig_inv <- diag(sigi_diag)


        V <- solve(t(X_t) %*% diag(1/P_t) %*% X_t + Sig_inv)
        axe_yhat <- X[test_mask, , drop = F] %*%
            V %*% t(X_t) %*% diag(1/P_t) %*%
            data$y[!test_mask]

        sd <- X[test_mask, , drop = F] %*%
            V %*% t(X[test_mask, , drop = F])
        list(
            loglik = dnorm(data$y[test_mask], mean = axe_yhat, sd = sqrt(as.vector(sd)), log = T),
            yhat = axe_yhat
        )


    })
    logliks <- sapply(logliks_y, function(l) l$loglik)
    yhat <- sapply(logliks_y, function(l) l$yhat)


    cv_lco_elppd <- -log(mean(exp(logliks)))

b <- Sys.time()


    list(idx = which(test_mask), cv_yhat = mean(cv_samps[test_mask, ]),
         sighat_mcv = mean(cv_taus), time_cv = difftime(b, a, units = "secs"),
         cv_elppd = cv_elppd, cv_lco = cv_lco_elppd, dif = mean(yhat) - mean(cv_samps[test_mask, ]))
}


crossval_radon <- function(loop, data, loops, X, formula, seed, ...) {
    a <- Sys.time()
    test_mask <- loop == loops
    rfit <- rstanarm::stan_glmer(
        formula,
        data = dplyr::filter(data, !test_mask),
        seed = seed, refresh = 0, ...
    )
    yhat <- rstanarm::posterior_predict(rfit, newdata = dplyr::filter(data, test_mask),
                                        draws = 500) %>% apply(2, mean)
    b <- Sys.time()
    stopifnot(length(yhat) == sum(test_mask))
    samples <- rstan::extract(rfit$stanfit)

    tau <- mean(samples$aux)
    sig2 <- mean(samples$theta_L)

    Y_t <- data$log_radon[test_mask]
    if ("beta" %in% names(samples)) {
        coeff_samps <- cbind(samples$alpha, samples$beta, samples$b)
    } else {
        coeff_samps <- cbind(samples$alpha,  samples$b)
    }

    # since test county is not provided, it is present in samples as last column in samples$b
    test_col <- which(apply(X, 2, function(col) all(col == as.numeric(test_mask))))
    ind <- X[, test_col]
    X <- cbind(X[, -test_col], ind)

    mean <- coeff_samps %*% t(X)
    elpd <- sapply(1:sum(test_mask), function(i) {
        log(mean(dnorm(Y_t[i], mean = mean[, i], sd = samples$aux)))
    })

    data.frame(
        idx = which(test_mask),
        yhat_cv = yhat,
        tauhat_cv = tau,
        time_cv = as.double(difftime(b, a, units = "secs")),
        sighat_cv = sqrt(sig2),
        elpd_cv = elpd
    )

}

crossval_radon_simul <- function(mod_no, dat, n, perc, subset_idx, info, seed, ...) {
    mod <- info$models_mcv[[mod_no]]

    test_mask <- dat$county == info$c_county
    train_mask <- !test_mask

    # fit to just training data
    start_time <- Sys.time()
    cv_fit <- rstanarm::stan_glmer(
        mod,
        data = dat[train_mask, ],
        seed = seed,
        refresh = 0, ...
    )
    yhat_cv <-  rstanarm::posterior_predict(
        cv_fit,
        newdata = dplyr::filter(dat, test_mask),
        draws = 500
    ) %>%
        apply(2, mean)

    end_time <- Sys.time()


    cv_samps <- rstan::extract(cv_fit$stanfit)

    # log lik
    elpd <- rstanarm::log_lik(cv_fit, newdata = dplyr::filter(dat, test_mask)) %>%
        exp %>% apply(2, mean) %>% log()


    X <- mm_radon(mod_no, dat, info$models_axe[[mod_no]], cv = T)

    county_col <- which(apply(X, 2, function(col) {
        all(col == as.numeric(test_mask))
    }))
    X <- cbind(X[, -county_col], X[county_col])
    colnames(X) <- c(colnames(X)[-ncol(X)], sprintf("county%s", info$c_county))
    coeffs <- coeffs_radon(cv_samps, cv = T)

    tau_samps <- cv_samps$aux
    lpd <- do.call("cbind", purrr::map(1:nrow(coeffs), function(mit) {
        dnorm(
            dat$log_radon[test_mask],
            mean = X[test_mask, ] %*% t(coeffs[mit, , drop = F]),
            sd = tau_samps[mit]
        )
    })) %>% apply(1, mean) %>% log()


    # cond_lp

    sigma2 <- mean(cv_samps$theta_L)
    sig_inv <- rep(0, ncol(X))
    sig_inv[stringr::str_detect(colnames(X), "county")] <- 1/sigma2
     cond_lpd <- axe_ll(diag(sig_inv), mean(tau_samps), X, test_mask, dat$log_radon)

    data.frame(
        y = dat$log_radon[test_mask],
        rst_elpd = elpd,
        elpd_cv = lpd,
        elpd_cv_cond = cond_lpd,

        cv_yhat = yhat_cv,
        tau2_cv =  mean(tau_samps)^2,
        sig2_cv = sigma2
    ) %>%
        dplyr::mutate(
            time_cv = as.double(difftime(end_time, start_time, units = "secs"))
        )



}




crossval_lol <- function(loop, data, loops, seed, ...) {
    test_mask <- loop == loops
    a <- Sys.time()
    cv_fit <- rstanarm::stan_glmer(
        kills ~ position  + team + (1 | player) + (1|champion) +
            log_dpm + log_egpm,
        family = "poisson", data = dplyr::filter(data, !test_mask), refresh = 0,
        seed = seed, refresh = 0, ...
    )

    yhats <- rstanarm::posterior_predict(
        cv_fit,
        newdata = dplyr::filter(data, test_mask)
    ) %>%
    apply(2, mean)
    b <- Sys.time()


    data.frame(
        time_cv = as.double(difftime(b, a, units = "secs")),
        loop = loop,
        idx = which(test_mask),
        yhat_cv = yhats
    )

}


crossval_slc <- function(loop, info, loops, seed, ...) {
    message(loop)
    test_mask <- loop == loops
    a <- Sys.time()

    D <- diag(rowSums(info$W))
    D_rr = D[!test_mask, !test_mask]
    W_rr = info$W[!test_mask, !test_mask]
    W_ir <- info$W[test_mask, !test_mask]

    X_cv <- rbind(info$X[!test_mask, ], info$X[test_mask, , drop = F])
    W_cv <- cbind(
        rbind(W_rr, W_ir),
        c((W_ir), 0)
    )
    loff_cv <- c(info$log_offset[!test_mask], info$log_offset[test_mask])
    sdat_cv <- list(n = info$n,
                 p = info$p, X = X_cv, y_train = info$y[!test_mask],
                 log_offset = loff_cv, W = W_cv, test_idx = info$n)
    cv_fit <- rstan::stan("inst/stan/carslp_cv.stan",
                        data = sdat_cv,
                        seed = seed, refresh = 0, iter = 10e3, thin = 5
    )
    b <- Sys.time()


    samples <- rstan::extract(cv_fit)


    set.seed(seed)


    yhats <- samples$y_test

    logliks <- sapply(yhats, function(y) {
        dpois(info$y[test_mask], yhats)
    }) %>% mean() %>% log()


    data.frame(
        time_cv = as.double(difftime(b, a, units = "secs")),
        loop = loop,
        idx = which(test_mask),
        yhat_cv = mean(yhats),
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

    chain1 <- CARBayesST::ST.CARar(formula = formula, family = "poisson",
                                   data = cvdat, W = info$W, burnin = 20000, n.sample = 220000,
                                   thin = 100)
    chain2 <- CARBayesST::ST.CARar(formula = formula, family = "poisson",
                                   data = cvdat, W = info$W, burnin = 20000, n.sample = 220000,
                                   thin = 100)
    chain3 <- CARBayesST::ST.CARar(formula = formula, family = "poisson",
                                   data = cvdat, W = info$W, burnin = 20000, n.sample = 220000,
                                   thin = 100)

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
        elpd_cv = logliks
    )

}


#####

# "..." captures parameters passed to underlying STAN model.
mcv_eight <- function(info = eight, n_cores = 7,
                      data_scale = seq(0.1, 4, by = 0.1),
                      seed = 238972, ...) {


    X <- model.matrix(~school, data = info$df,
                      contrasts.arg = list(school = contrasts(info$df$school, contrasts = F)))
    loops <- 1:8
    purrr::map(data_scale, function(scl) {
        scld_dat <- scale_dat(info$stan_list, scl)

        if (n_cores > 1) {
            cl <- snow::makeCluster(n_cores, outfile = '')
            snow::clusterExport(cl, list("%>%")) # pkg: change to envir = environment()

            yhats <- snow::parLapply(cl,loops, crossval_eight,
                                     scld_dat, loops, X, seed)

            snow::stopCluster(cl)
        } else {
            yhats <- lapply(loops, crossval_eight,
                            scld_dat, loops, X, seed, ...)
        }
        objs <- .collate_folds(yhats)
        c(objs, data_scale = scl)

    })
}

mcv_radon_full <- function(info = radon_1, models = radon_1$models_mcv,
                           seed = 32897, n_cores = 7, ...) {
    xdat <- info$data %>%
        dplyr::mutate(
            floor = factor(floor), county = factor(county)
        )
    loops <- as.character(info$data$county)
    y <- info$data$log_radon

    purrr::map(1:length(models), function(mod_no) {
        print(mod_no)
        X <- model.matrix(info$models_axe[[mod_no]], data = xdat,
                          contrasts = contrasts_for_pooling(xdat))

        if (n_cores > 1) {
            cl <- snow::makeCluster(n_cores, outfile = '')
            snow::clusterExport(cl, list("X", "%>%"),
                                envir = environment())

            objs <- snow::parLapply(cl, unique(loops), crossval_radon,
                                    info$data, loops,
                                    X, models[[mod_no]], seed, ...)

            snow::stopCluster(cl)
        } else {
            objs <- lapply(unique(loops), crossval_radon,
                           info$data, loops,
                           X, models[[mod_no]], seed)
        }
        objs <- dplyr::bind_rows(objs)

        objs %>%
            dplyr::arrange(idx) %>%
            dplyr::mutate(
                loop = as.numeric(factor(loops))
            )
    })


}


mcv_radon_simul <- function(info = radon_2,
                            seed = 9871, n_cores = 5, ...) {

    loop_over_radon_simul(
        info = info, seed = seed,n_cores = n_cores,
        export = list("mm_radon", "coeffs_radon", "axe_ll", "axe", "axe_v"),
        expenv = environment(),
        FUN = crossval_radon_simul,

        ...)
}


######




mcv_lol <- function(info = lol,
                    seed = 23982, n_cores = 7, ...) {
    loops <- info$loops

    if (n_cores > 1) {
        cl <- snow::makeCluster(n_cores, outfile = '')
        snow::clusterExport(cl, list("%>%"),
                            envir = environment())

        objs <- snow::parLapply(cl, unique(loops), crossval_lol,
                                info$data, loops,  seed, ...)

        snow::stopCluster(cl)
    } else {
        objs <- lapply(unique(loops), crossval_lol,
                       info$data, loops, seed, ...)
    }

    objs

}


mcv_slc <- function(info = slc, seed = 29873, n_cores = 3, ...) {
    loops <- info$data$loops

    if (n_cores > 1) {
        cl <- snow::makeCluster(n_cores, outfile = '')
        snow::clusterExport(cl, list("%>%", "stanmodels"),
                            envir = environment())

        objs <- snow::parLapply(cl, unique(loops), crossval_slc,
                                info$data, loops,  seed, ...)

        snow::stopCluster(cl)
    } else {
        objs <- lapply(unique(loops), crossval_slc,
                       info$data, loops, seed, ...)
    }

    dplyr::bind_rows(objs)
}



mcv_air <- function(info = air, seed = 69872, n_cores = 3, ...) {
    loops <- info$data$loops

    if (n_cores > 1) {
        cl <- snow::makeCluster(n_cores, outfile = '')
        snow::clusterExport(cl, list("%>%"),
                            envir = environment())

        objs <- snow::parLapply(cl, unique(loops), crossval_air,
                                info$data, loops,  seed, info$data$formula, ...)

        snow::stopCluster(cl)
    } else {
        objs <- lapply(unique(loops), crossval_air,
                       info$data, loops, seed, info$data$formula, ...)
    }

    dplyr::bind_rows(objs)
}




