# for eight schools data
scale_dat <- function(dat, scl) {
    dat$y <- dat$y * scl
    dat
}


contrasts_for_pooling <- function(df, categorical_names = NULL) {
    if (is.null(categorical_names)) {
        categorical_mask <- purrr::map_lgl(df, function(x) {
            is.character(x) | is.factor(x)
        })
        categorical_names <- colnames(df)[categorical_mask]
    }
    purrr::map(categorical_names, function(x) {
        stats::contrasts(df[[x]], contrasts = F)
    }) |> purrr::set_names(categorical_names)
}


# convenience function for CAR model precision matrix
car_inv <- function(D, W, alpha, tau) {
    tau^2 * (D - alpha * W)
}


# convenience function for getting Sigma^-1 (for AXE) with spatio-temporal CAR
# as defined in paper
stcarar_sigmainv <- function(rho_s, rho_t, W, tau2, P = 4 + N,
                             N = nrow(air$data$df), J = N / Ts, Ts = 5) {
    sigmai_spat <- (rho_s * (diag(rowSums(W)) - W)) +
        diag(1 - rho_s, nrow = nrow(W)) / tau2
    sigmai_blck <- matrix(
        sapply(1:Ts, function(t) {
            sapply(1:ncol(sigmai_spat), function(j) {
                c(rep(0, J * (t - 1)), sigmai_spat[, j], rep(0, J * (Ts - t)))
            })
        }),
        ncol = N
    )
    Sigma_inv <- matrix(0, ncol = P, nrow = P)
    sigma_idxs <- (P - N + 1):P

    H <- matrix(0, N, N)
    H[(J + 1):N, 1:(N - J)] <- rho_t * diag(rep(1, N - J))
    I <- diag(rep(1, N))


    Sigma_inv[sigma_idxs, sigma_idxs] <- t(I - H) %*% sigmai_blck %*% (I - H)

    Sigma_inv
}


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

sigma_diagonal <- function(X, sigs, extra) {
    xmatch <- match_to_x(names(sigs), X)
    reff_mask <- !is.na(xmatch)

    sigi_diag <- rep(0, ncol(X))
    sigs <- purrr::map_dbl(sigs, ~ 1 / .)
    sigi_diag[reff_mask] <- sigs[xmatch[reff_mask]]
    diag(sigi_diag)
}


sigma_car <- function(X, sigs, extra) {
    prec <- (sigs$tau) * (extra$D - sigs$alpha * extra$W)
    fixed <- ncol(X) - ncol(prec)
    Sig_inv <- matrix(0, ncol(X), ncol(X))
    Sig_inv[-c(1:fixed), -c(1:fixed)] <- prec

    Sig_inv
}
