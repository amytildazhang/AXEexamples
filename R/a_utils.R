# variance of normal approximation of binomial with probit link
approx_bin_variance <- function(p, n) {
  p * (1 - p) / (n * dnorm(qnorm(p))^2)
}


# for eight schools data
scale_dat <- function(dat, scl) {
  dat$y <- dat$y * scl
  dat
}


contrasts_for_pooling <- function(df, categorical_names = NULL) {
  if (is.null(categorical_names)) {
    categorical_mask <- purrr::map_lgl(df, ~ is.character(.) |
      is.factor(.))
    categorical_names <- colnames(df)[categorical_mask]
  }
  purrr::map(categorical_names, ~ stats::contrasts(df[[.]],
    contrasts = F
  )) %>% purrr::set_names(categorical_names)
}


# convenience function for CAR model precision matrix
car_inv <- function(D, W, alpha, tau) {
  tau^2 * (D - alpha * W)
}


# convenience function for getting Sigma^-1 (for AXE) with spatio-temporal CAR
# as defined in paper
stcarar_sigmainv <- function(rho_s, rho_t, W, tau2, P = 4 + N,
                             N = nrow(air$data$df), J = N / Ts, Ts = 5) {
  sigmai_spat <- (rho_s * (diag(rowSums(W)) - W) +
    diag(1 - rho_s, nrow = nrow(W))) / tau2
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


#
# cond_lpd_radon <- function(model, sig, tau, dat, test_mask, schem, cv = F) {
#     X <- stats::model.matrix(info$models_axe[[mod_no]], data = dat,
#                              contrasts = contrasts_for_pooling(dat))
#     county_col <- which(apply(X, 2, function(col) {
#         all(col == as.numeric(test_mask))
#     }))
#     if (cv) {
#         X <- X[, -county_col]
#     }
#     train_mask <- !test_mask
#
#
#     # cond_lp
#     X_t <- X[train_mask, ]
#     n_0 <- ncol(X_t) - ncol(schem)
#     siginv <- diag(c(rep(0, n_0), rep(1/sigma2, ncol(schem))))
#     V_t <- solve((t(X_t) %*% X_t)/tau^2 + siginv)
#     axe_yhat <- (X[test_mask, , drop = F] %*% V_t %*% t(X_t) %*% dat$log_radon[train_mask])/tau^2
#     Sigma <- X[test_mask, ] %*% V_t %*% t(X[test_mask, ])
#     # assume independence rather than drawing from multivariate normal;
#     # diagonals are not > off-diagonals because of repeated measures
#
#     dnorm(
#         dat$log_radon[test_mask],
#         mean = axe_yhat,
#         sd = sqrt(diag(Sigma)), log = T
#     )
#
# }



#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`
