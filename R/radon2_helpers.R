
radon2_prefix <- function(perc, n, iter, mod_no) {
  sprintf("perc%s-n%s-iter%s-mod%s",
          perc * 10, n, iter, mod_no)
}
radon_simul_samples <- function(perc, n, iter, mod_no, cv = FALSE) {
  topfolder <- file.path("pretrained", "r2")
  if (cv) {
    topfolder <- file.path(topfolder, "mcv")
  }
  file.path(topfolder, radon2_prefix(perc, n, iter, mod_no))
}

loop_inner <- function(n, info, seed, FUN, ...) {
  purrr::map_df(info$data_perc, function(perc) {
    schem <- info$schema[[sprintf("perc%s", perc)]][[sprintf("cluster%s", n)]]

    purrr::map_df((1:nrow(schem)), function(subset_idx) {
      dat <- info$data |>
        dplyr::filter(county %in% c(info$c_county, schem[subset_idx, ])) |>
        dplyr::mutate(county = factor(county), floor = factor(floor))



      purrr::map_df(seq_along(info$models_mcv), function(mod_no) {
        FUN(mod_no, dat, n, perc, subset_idx, info, seed, ...) |>
          dplyr::mutate(
            model = mod_no,
            idx = 1:dplyr::n(),
            iter = subset_idx,
            n_clusters = n,
            perc = perc
          )
      })
    })
  })
}


loop_over_radon_simul <- function(info = radon_2, seed = 9871, n_cores = 5,
                                  export = NULL, expenv = environment(), FUN,
                                  ...) {
  if (n_cores > 1) {
    cl <- snow::makeCluster(n_cores, outfile = "")
    snow::clusterExport(cl, c(
      "info", "seed", "|>", "contrasts_for_pooling",
      export
    ),
    envir = expenv
    )
    snow::clusterCall(cl, function(x) {
      library(tidyverse)
      for (file in list.files("R")) {
        source(file.path("R", file))
      }

    })
  }

  set.seed(seed)


  if (n_cores > 1) {
    tryCatch(
      {
        do.call(
          dplyr::bind_rows,
          snow::parLapply(cl, info$n_clusters, loop_inner,
                          info, seed, FUN, ...)
        )
      },
      finally = snow::stopCluster(cl)
    )
  } else {
    purrr::map_df(info$n_clusters, loop_inner,
                  info = info, seed = seed, FUN = FUN, ...
    )
  }
}


samples_radon <- function(perc, n, subset_idx, mod_no) {
  sf <- radon_simul_samples(perc, n, subset_idx, mod_no)
  readr::read_rds(sf)
}

mm_radon <- function(mod_no, dat, model, cv = FALSE, cv_fold = "OLMSTED") {
  X <- stats::model.matrix(model,
                           data = dat,
                           contrasts = contrasts_for_pooling(dat)
  )

  if (mod_no != 1 && !cv) {
    X[, 1] <- 1
  }
  X
}

coeffs_radon <- function(samples, cv = FALSE) {
  if ("beta" %in% names(samples)) {
    coeffs <- cbind(samples$alpha, samples$beta, samples$b)
  } else {
    coeffs <- cbind(samples$alpha, samples$b)
  }
  if (!cv) {
    coeffs <- coeffs[, -ncol(coeffs)]
  }

  coeffs
}


#
# Radon subsets -- load saved posterior samples and then obtain approximations
#
lco_radon_simul <- function(info = radon_2, seed = 8971, n_cores = 5,
                            method = c("iis", "ghost", "vt")) {
  loop_over_radon_simul(
    info, seed, n_cores,
    export = list(
      "mm_radon", "coeffs_radon", "radon_simul_samples",
      "iis", "ghosting_c", "iis.loglik", "match_to_x"
    ),
    expenv = environment(),
    FUN = function(mod_no, dat, n, perc, subset_idx, info, seed, method) {
      sf <- radon_simul_samples(perc, n, subset_idx, mod_no)
      samples <- readr::read_rds(sf)
      coeffs <- coeffs_radon(samples)

      X <- mm_radon(mod_no, dat, info$models_axe[[mod_no]])

      test_mask <- dat$county == info$c_county
      cv_loop <- rep(1, sum(test_mask))


      X_test <- X[test_mask, , drop = FALSE]
      Y_test <- dat$log_radon[test_mask]
      Xbeta_samps <- coeffs %*% t(X_test)

      if (method %in% c("iis", "ghost")) {
        X_effs <- X_test[, stringr::str_detect(colnames(X), "county")]
        mu_samps <- Xbeta_samps - samples$b[, -c(n + 1)] %*% t(X_effs)

        X_b <- X_test[, sprintf("county%s", info$c_county), drop = FALSE]

        if (method == "iis") {
          iis(cv_loop,
              mu_samps,
              X_b = X_b, Y = Y_test,
              sig_samps = samples$theta_L, tau_samps = samples$aux^2
          )
        } else if (method == "ghost") {
          ghosting_c(
            seed, cv_loop,
            mu_samps = mu_samps,
            sig_samps = samples$theta_L, X_b = X_b,
            Y = info$data$log_radon[test_mask], tau_samps = samples$aux^2
          )
        }
      } else if (method == "vt") {
        vehtari(loops = dat$county %in% info$c_county,
                Xbeta_samps, X = X, Y = dat$log_radon,
                sig_samps = list("county" = samples$theta_L),
                tau_samps = samples$aux^2, loop = TRUE,
                Sigma = "diagonal", family = "normal",
                glmm = FALSE, beta_samps = coeffs
        )

      }
    }, method = method
  )
}
