hess_pois <- function(X, beta, y, mu = 0, ...) {
  mult <- -exp(as.vector(X %*% beta) + mu)
  if (length(mult) > 1) {
    t(X) %*% diag(mult) %*% X
  } else {
    t(X) %*% X * mult
  }
}

gr_pois <- function(X, beta, y, mu = 0) {
  # mu is, for example, a log-offset term
  dif <- y - exp(X %*% beta + mu)
  apply(diag(dif) * X, 2, sum)
}

ll_pois <- function(X, beta, y, mu = 0, ...) {
  dpois(y, exp(X %*% beta + mu), log = TRUE)
}

hess_norm <- function(beta, Sigma_inv) {
  # treating Sigma_inv as constant
  -2 * Sigma_inv
}
gr_norm <- function(beta, Sigma_inv, ...) {
  if (length(beta) > 1) {
    # assumes beta^tSigma^-1beta kernel, i.e. prior means are 0
    - 2 * Sigma_inv %*% beta

  } else {
    -2 * Sigma_inv * beta
  }
}

ll_norm <- function(beta, mu, Sigma_inv, ...) {
  if (length(Sigma_inv) > 1) {
    dif <- beta - mu
    det <- determinant(Sigma_inv, logarithm = TRUE)$modulus
    beta <- as.matrix(beta, ncol = 1)
    as.vector(det - length(beta) / 2 * log(2 * pi) -
                t(dif) %*% Sigma_inv %*% dif / 2)

  } else {
    dnorm(beta, mean = mu, sd = sqrt(1 / Sigma_inv), log = TRUE)
  }

}


.lpl_approx_pois <- function(beta0, Sigma_inv, X, y, mu) {
  hess <- hess_pois(X, beta0, y, mu) + hess_norm(beta0, Sigma_inv)
  ll_norm(beta0, mu = beta0, Sigma_inv = - hess)
}

lpl_approx_pois <- function(y, X, beta, Sigma_inv, mu, tau_inv = NULL,
                            dist = "poisson") {
  model_fun <- function(pars, y, Sigma_inv, X, mu, ...) {
    sum(ll_norm(pars, 0, Sigma_inv)) + sum(ll_pois(X, pars, y, mu))
  }

  grad_fun <- function(pars, y, Sigma_inv, X, mu, ...) {
    gr_norm(pars, Sigma_inv) + gr_pois(X, pars, y, mu)
  }

  inits <- beta
  mode <- optim(
    inits, model_fun, control = list(fnscale = -1),
    gr = grad_fun,
    y = y, Sigma_inv = Sigma_inv, X = X, mu = mu, tau_inv = tau_inv,
    method = ifelse(length(beta > 1), "BFGS", "Brent")
  )

  list(
    ll =  .lpl_approx_pois(beta0 = mode$par, X = X,
                           Sigma_inv = Sigma_inv, y = y, mu = mu),
    beta = mode$par
  )

}


.lpl_approx_norm <- function(beta0, Sigma_inv, X, y, mu,
                             tau_mat = NULL, tau_inv = NULL) {
  if (is.null(tau_mat)) {
    hess <- -2 * tau_inv * t(X) %*% X + hess_norm(beta0, Sigma_inv)

  } else {
    hess <- -2 * t(X) %*% tau_mat %*% X + hess_norm(beta0, Sigma_inv)
  }
  ll_norm(beta0, mu = beta0, Sigma_inv = - hess)
}


lpl_approx_norm <- function(y, X, beta, Sigma_inv, mu,
                            tau_inv = NULL, dist = "normal") {

  if (length(tau_inv) > 1) {
    tau_mat <- diag(tau_inv, nrow = nrow(X), ncol = nrow(X))
    tau_inv <- tau_mat
  } else {
    tau_mat <- NULL
  }
  model_fun <- function(pars, y, Sigma_inv, X, tau_inv, ...) {
    eta <- X %*% pars
    sum(ll_norm(pars, 0, Sigma_inv)) +
      sum(ll_norm(beta = y, mu = eta, Sigma_inv = tau_inv))
  }

  grad_fun <- function(pars, y, Sigma_inv, X,
                       tau_mat = NULL, tau_inv = NULL, ...) {
    if (is.null(tau_mat)) {
      as.vector(gr_norm(pars, Sigma_inv) +
                  2 * tau_inv * t(X) %*% (y - X %*% pars))

    } else {
      as.vector(gr_norm(pars, Sigma_inv) +
                  2 * t(X) %*% tau_mat %*% (y - X %*% pars))
    }
  }


  inits <- beta
  mode <- optim(
    inits, model_fun, control = list(fnscale = -1),
    gr = grad_fun, tau_mat = tau_mat,
    y = y, Sigma_inv = Sigma_inv, X = X, mu = mu, tau_inv = tau_inv,
    method = ifelse(length(beta > 1), "BFGS", "Brent")
  )

  list(
    ll =  .lpl_approx_norm(beta0 = mode$par, X = X, Sigma_inv = Sigma_inv,
                           y = y, mu = mu, tau_mat = tau_mat),
    beta = mode$par
  )

}

lpl_approx.iis_im <- function(y, X, beta, Sigma_inv, mu, dist = "poisson") {
  # laplace approximation to get f(beta | y_train, Sigma, tau)

  hess_loglik <- hess_pois(X = X_test, beta = beta_mode, y = Y_test, mu = mu)
  V_inv <- t(X) %*% diag(1 / P) %*% X + diag(-1 / diag(hess_loglik))
  ll_norm(beta, mode = beta_mode, Sigma_inv = V_inv)

}
