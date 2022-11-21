# Run and save posterior samples for all models/examples

for (file in list.files("data")) {
    load(sprintf("data/%s", file))
}
pfit_eight <- function(info = eight, seed = 238972, ...) {
    mod <- get_cmdstan_obj("eight")
    savefolder <- "pretrained/eight_schools"
    purrr::map((info$data_scale), function(scl) {
        message(scl)
        scld_dat <- scale_dat(info$stan_list, scl)
        # fit to full data
        pfit <- mod$sample(
            data = scld_dat,
            seed = seed, refresh = 0,
            parallel_chains = 4,
            ...)
        # fnames <- pfit$save_output_files(dir = savefolder,
        #                                  basename = scl)
        #

        a <- pfit$draws(format = "draws_matrix")
        map <- unclass(a[which.max(a[, "lp__"]), ])[1,]
        map_df <- data.frame(
            param = names(map),
            map = map
        ) %>%
            filter(param %in% c("mu", "tau") | str_detect(param, "theta"))

        write.csv(map_df, file = sprintf("%s/%.1f_map.csv", savefolder, scl),
                  row.names = FALSE)

        pfit
    })
}



stanglmer_samps <- function(rfit, exclude = '') {
    samples <- rstan::extract(rfit$stanfit)

    map_idx <- which.max(samples$lp__)
    if ("beta" %in% names(samples)) {
        coeffs <- cbind(samples$alpha, samples$beta, samples$b)
    } else {
        coeffs <- cbind(samples$alpha, samples$b)
    }

    sampnames <- names(rfit$stanfit)
    exclude_pars <- c("mean_PPD", "log-posterior")
    exclude_random <- (stringr::str_detect(sampnames, "_NEW_")) #|
        # stringr::str_detect(sampnames, exclude)) &
        # !stringr::str_detect(sampnames, 'Sigma')

    parnames <- setdiff(sampnames[!exclude_random], exclude_pars)
    samps <- cbind(coeffs[, !exclude_random[1:ncol(coeffs)]],
                   samples$theta_L, samples$aux)
    colnames(samps) <-parnames
    samps[map_idx, ]
}

save_prior_json <- function(list, savefolder, prefix) {
    jsd <- rjson::toJSON(list)
    fileConn <- file(sprintf("%s/%sprior.json", savefolder, prefix))
    writeLines(jsd, fileConn)
    close(fileConn)

}

save_map <- function(map, savefile) {
    map_df <- data.frame(
        param = names(map),
        map = map
    )
    write.csv(map_df, savefile,
              row.names = FALSE)

}
save_stanglmer_info <- function(savefolder, prefix, rfit, exclude = '') {
    # save prior info
    save_prior_json(rfit$prior.info, savefolder, prefix )
    map <- (stanglmer_samps(rfit, exclude))
    save_map(map, sprintf("%s/%smap.csv", savefolder, prefix))
}


pfit_r1 <- function(info = radon_1, seed = 32897, ...) {
    savefolder <- "pretrained/r1"
    xdat <- info$data %>%
        dplyr::mutate(floor = factor(floor), county = factor(county))
    write.csv(xdat, sprintf("%s/data.csv", savefolder))

    purrr::map(1:length(info$models_mcv), function(mod_no) {
        # message(sprintf("Model %s", mod_no))
        fname <- sprintf("%s/model%s", savefolder, mod_no)
        rfit <- rstanarm::stan_glmer(info$models_mcv[[mod_no]],
                                     data = info$data, seed = seed,
                                     refresh = 0, sample_file = fname)
        save_stanglmer_info(savefolder, sprintf("model%s_", mod_no), rfit,
                            exclude = 'county')
        # save X
        X <- model.matrix(info$models_axe[[mod_no]],
                          data = xdat,
                          contrasts = contrasts_for_pooling(xdat)
        )
        write.csv(X, file = sprintf("%s/model%s_X.csv", savefolder, mod_no), row.names = FALSE)

    })
}


pfit_r2 <- function(info = radon_2,
                    seed = 9871, n_cores = 5) { # arguments to pass to rstanarm::stan_g

    loop_over_radon_simul(
        info, seed, n_cores,
        export = list("mm_radon", "coeffs_radon", "radon_simul_samples", "radon2_prefix",
                      "save_stanglmer_info", "save_prior_json", "stanglmer_samps"),
        expenv = environment(),
        FUN = function(mod_no, dat, n, perc, subset_idx, info, seed) {
            sf <- radon_simul_samples(perc, n, subset_idx, mod_no)
            message(sf)
            mod <- info$models_mcv[[mod_no]]
            a <- Sys.time()
            rfit <- rstanarm::stan_glmer(mod, data = dat, seed = seed,
                                         refresh = 0, iter = 1e4, thin = 5)
            b <- Sys.time()
            prefix <- sprintf("%s_", radon2_prefix(perc, n, subset_idx, mod_no))
            savefolder <- "pretrained/r2"
            save_stanglmer_info(savefolder, prefix, rfit, exclude = 'county')

            # save X
            dat <- dplyr::mutate(dat, floor = factor(floor), county = factor(county))
            X <- model.matrix(info$models_axe[[mod_no]],
                              data = dat,
                              contrasts = contrasts_for_pooling(dat)
            )
            write.csv(X, file = sprintf("%s/%sX.csv", savefolder, prefix), row.names = FALSE)
            write.csv(dat, file = sprintf("%s/%sdata.csv", savefolder, prefix), row.names = FALSE)

            data.frame(
                file = sf,
                time_post = difftime(b, a, units = "secs")
            )
        })
}

pfit_lol <- function(info = lol, seed = 18973, ...) {
    savefolder <- "pretrained/lol"
    pfit <- rstanarm::stan_glmer(kills ~ position + team +
                                     log_dpm + log_egpm +
                                     (1 | player) + (1 | champion),
                                 family = "poisson", data = info$data,
                                 seed = seed,
                                 sample_file = sprintf("%s/lol", savefolder)
    )

    # save X

    write.csv(1/pfit$fitted.values, sprintf("%s/phi.csv", savefolder, row.names = FALSE))
    write.csv(info$X, sprintf("%s/X.csv", savefolder), row.names = FALSE)
    # save data frame
    write.csv(info$data, sprintf("%s/data.csv", savefolder), row.names = FALSE)
    save_stanglmer_info(savefolder, "", pfit, exclude = 'player')
    # save samples
}

pfit_slc <- function(info = slc$data, seed = 29873, ...) {
    savefolder <- "pretrained/slc"
    cmdmod <- get_cmdstan_obj("slc")
    sdat <- list(
        n = info$n, p = info$p, X = info$X, y = info$y,
        log_offset = info$log_offset, W_n = info$W_n, W = info$W,
        D_sparse = rowSums(info$W)
    )

    pfit <- cmdmod$sample(data = sdat, refresh = 0, seed = seed,
                          iter_warmup = 2e3, iter_sampling = 2e3,
                          ...)

    slc_prior <- list(
        'prior' = NULL,
        'prior_intercept' = list(
            location = rep(0.0, info$p),
            scale = rep(1.0, info$p),
            adjusted_scale = NULL,
            dist = 'normal'
        ),
        'prior_tau' = list(
            dist = 'gamma',
            shape = 2.0,
            scale = 0.5
        ),
        'prior_alpha' = list(
            dist = 'uniform'
        )
    )

    save_prior_json(slc_prior, savefolder, '')

    samps <- pfit$draws(format = "draws_matrix")
    map_idx <- which.max(samps[, "lp__"])
    map <- samps[map_idx, -1]
    save_map(savefolder, "", unclass(map)[1,])

    X <- cbind(info$X, diag(rep(1, nrow(info$X))))
    write.csv(X, sprintf("%s/X.csv", savefolder), row.names = FALSE)

    samps <- unclass(samps)
    betas <- samps[, str_detect(colnames(samps), "beta|phi")] %>% apply( 2, mean)
    phats <- exp(X %*% betas + info$log_offset)
    names(phats) <- NULL
    write.csv(1/phats, sprintf("%s/phi.csv", savefolder, row.names = FALSE))

    write.csv(info$W, sprintf("%s/W.csv", savefolder), row.names = FALSE)

    write.csv(data.frame(y = info$y, log_offset = info$log_offset, loop = info$loops ),
              sprintf("%s/data.csv", savefolder), row.names = FALSE)
}


sampnames <- function(name, matrix) {
    if (name == "rho") {
        labels <- colnames(matrix) %>% str_replace(".", "_")
    } else if (ncol(matrix) > 1) {
        labels <- sprintf("%s%s", name, 1:ncol(matrix))
    } else {
        labels <- name
    }
    colnames(matrix) <- labels
    matrix
}
pfit_air <- function(info = air$data, seed = 69872, ...) {
    set.seed(seed)
    savefolder <- "pretrained/air"
    formula <- info$formula
    samples <- lapply(1:4, function(chain) {
        fit <- CARBayesST::ST.CARar(
            formula = formula, family = "poisson",
            data = info$df, W = info$W, burnin = 20000, n.sample = 220000,
            thin = 200, AR = 1
        )
        samps <- sampnames("beta", fit$samples$beta)
        for (par in c("phi", "rho", "tau2")) {
            parsamps <- sampnames(par, fit$samples[[par]])
            samps <- cbind(samps, parsamps)
        }
        samps
    })

    samps <- do.call("rbind", samples)
    betas <- str_detect(colnames(samps), "beta")
    phis <- str_detect(colnames(samps), "phi")

    Xbeta <- info$X %*% t(samps[, betas])

    phats_all <- Xbeta + t(samps[, phis])
    lp <- apply(phats, 2, function(phs) {
        mus = exp(phs + info$log_offset)
        sum(dpois(info$df$observed, mus, log = TRUE))

    })

    phats <- exp(apply(phats_all, 1, mean) + info$log_offset)
    write.csv(1/phats, sprintf("%s/phi.csv", savefolder, row.names = FALSE))

    map_idx <- which.max(lp)
    map <- samps[map_idx, ]
    save_map(map, sprintf("%s/map.csv", savefolder))


    write.csv(info$df, sprintf("%s/data.csv", savefolder), row.names = FALSE)
    write.csv(info$W, sprintf("%s/W.csv", savefolder), row.names = FALSE)



        write.csv(info$X, sprintf("%s/X.csv", savefolder), row.names = FALSE)

    air_prior <- list(
        'prior' = NULL,
        'prior_intercept' = list(
            location = rep(0.0, ncol(info$X)),
            scale = rep(1000.0, ncol(info$X)),
            adjusted_scale = NULL,
            dist = 'normal'
        ),
        'prior_tau2' = list(
            dist = 'invgamma',
            shape = 1.0,
            scale = 0.01
        ),
        'prior_rho_s' = list(
            dist = 'uniform'
        ),
        'prior_rho_t' = list(
            dist = 'uniform'
        )
    )
    save_prior_json(air_prior, savefolder, '')

}

#
a <- Sys.time()
eight$pfit_files <- pfit_eight()
b <- Sys.time()

difftime(b, a, units = "secs")

a <- Sys.time()
radon_1$pfit_files <- pfit_r1()
b <- Sys.time()
difftime(b, a, units = "secs")

a <- Sys.time()
radon_2$pfit_files <- pfit_r2()
b <- Sys.time()
 difftime(b, a, units = "secs")
system.time(pfit_r1())
system.time(pfit_r2())
system.time(pfit_lol())
system.time(pfit_slc())
system.time(pfit_air())
