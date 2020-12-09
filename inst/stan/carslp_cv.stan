// all matrices and data provided must be ordered so that test_idx = n

data {
  int<lower = 1> n;
  int<lower = 1> p;
  matrix[n, p] X;
  int<lower = 0> y_train[n - 1];
  vector[n] log_offset;
  matrix<lower = 0, upper = 1>[n, n] W;

  int<lower=1, upper=n> test_idx;
}
transformed data{
  vector[n] zeros;
  matrix<lower = 0>[n, n] D;
  {
    vector[n] W_rowsums;
    for (i in 1:n) {
      W_rowsums[i] = sum(W[i, ]);
    }
    D = diag_matrix(W_rowsums);
  }
  zeros = rep_vector(0, n);
}
parameters {
  vector[p] beta;
  vector[n] phi;
  real<lower = 0> tau;
  real<lower = 0, upper = 1> alpha;
}

transformed parameters {
  vector[n] mu;
  vector[n - 1] mu_train;

  mu  = X * beta + phi + log_offset;
 mu_train = mu[1:(n - 1)];
}
model {
  phi ~ multi_normal_prec(zeros, tau * (D - alpha * W));
  beta ~ normal(0, 1);
  tau ~ gamma(2, 2);
  y_train ~ poisson_log(mu_train);
}

generated quantities {
  int y_test = poisson_log_rng(mu[n]);
}
