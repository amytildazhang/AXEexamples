// saved as 8schools.stan
data {
    int<lower=1> J;         // number of schools
    int<lower=1> N;         // total sample size
    int<lower=0, upper=N> H; // number held-out


    real y[N];              // estimated treatment effects
    real<lower=0> sigma[J]; // standard error of effect estimates

    int<lower=0, upper=J> sigma_idx[N];

    int<lower=1, upper=N> train_idx[N-H];
    int<lower=1, upper=N> test_idx[H];
}

transformed data {
    int<lower=1, upper=J> school_train[N-H] = sigma_idx[train_idx];
    int<lower=1, upper=J> school_test[H] = sigma_idx[test_idx];

}

parameters {
    real mu;                // population treatment effect
    real<lower=0> tau;      // standard deviation in treatment effects
    vector[J] eta;          // unscaled deviation from mu by school
}
transformed parameters {
    vector[J] theta = tau * eta;        // school treatment effects
}

model {
    eta ~ normal(0, 1);
    y[train_idx] ~ normal(mu + theta[school_train], sigma[school_train]);
}

generated quantities {
    vector[N] log_lik;

    for (n in 1:N) {
        log_lik[n] = normal_lpdf(y[n] | mu + theta[sigma_idx[n]], sigma[sigma_idx[n]]);
    }

}
