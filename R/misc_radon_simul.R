

radon_simul_samples <- function(perc, n, iter, mod_no) {
    sprintf("data-raw/radon2_samples/perc%s_n%s_i%s_model%s.rds", perc, n, iter, mod_no)
}

loop_inner <- function(n, info, seed, FUN, ...) {
    purrr::map_df(info$data_perc, function(perc) {

        schem <- info$schema[[sprintf("perc%s", perc)]][[sprintf("cluster%s", n)]]

        purrr::map_df(1:nrow(schem), function(subset_idx) {
            dat <- dplyr::filter(info$data, county %in% c(info$c_county, schem[subset_idx, ])) %>%
                dplyr::mutate(county = factor(county), floor = factor(floor))



            purrr::map_df(1:length(info$models_mcv), function(mod_no) {
                FUN(mod_no, dat, n, perc, subset_idx, info, seed, ...) %>%
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


loop_over_radon_simul <- function(info = radon_2,
                                  seed = 9871, n_cores = 5,
                                  export = NULL, expenv = environment(), FUN, ...) {

    if (n_cores > 1) {
        cl <- snow::makeCluster(n_cores, outfile = '')
        snow::clusterExport(cl, c("info", "seed", "%>%", "contrasts_for_pooling",
                                  export),
                            envir = expenv)
    }

    set.seed(seed)


    if (n_cores > 1) {
        tryCatch({
            do.call(
                dplyr::bind_rows,
                snow::parLapply(cl, (info$n_clusters), loop_inner, info, seed, FUN, ...)
            )
        }, finally = snow::stopCluster(cl))

    } else {
        purrr::map_df(info$n_clusters, loop_inner,
                      info = info, seed = seed, FUN = FUN, ...)
    }




}


samples_radon <- function(perc, n, subset_idx, mod_no) {
    sf <-  radon_simul_samples(perc, n, subset_idx, mod_no)
    readr::read_rds(sf)

}

mm_radon <- function(mod_no, dat, model, cv = F, cv_fold = "OLMSTED") {
    X <- stats::model.matrix(model, data = dat,
                             contrasts = contrasts_for_pooling(dat))

    if (mod_no != 1 & !cv) {
        X[, 1] <- 1
    }
    X

}

coeffs_radon <- function(samples, cv = F) {
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

