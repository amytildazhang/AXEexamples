
#
#  GHOST
#

ghosting_c <- function(seed, loops, mu_samps, sig_samps,
                       X_b, Y, tau_samps = NULL, P = NULL,
                       Sigma = c("diagonal", "car"), M = NULL, glmm = FALSE,
                       log_offset = 0) {
  set.seed(seed)
  start <- Sys.time()

  J <- length(unique(loops))

  b_samps <- do.call(
    "rbind",
    purrr::map(seq_along(sig_samps), function(s) {
      rnorm(J, sd = sqrt(sig_samps[s]))
    }))

  # S x N + S x J x N
  yhats <- mu_samps + b_samps %*% t(X_b)

  if (glmm) {
    yhats <- exp(yhats + log_offset)
  }
  end <- Sys.time()
  time_yhat <- as.double(difftime(end, start, units = "secs"))
  end <- Sys.time()


  data.frame(
    time_ghst_c = (time_yhat) / J,
    yhat_ghst_c = apply(yhats, 2, mean),
    var_ghst = apply(yhats, 2, var)
  )
}
