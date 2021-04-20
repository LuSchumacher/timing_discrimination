data {
  int<lower=1>                   T;         // number of trials 
  int<lower=1>                   S;         // number of subjects
  int<lower=1>                   sub[T];    // subject number
  int                            resp[T];   // response
  vector<lower=0>[T]             rt;        // response time
  vector[T]                      d1;        // stimulus 1
  vector[T]                      d2;        // stimulus 2
}


parameters {
  //hyper parameter
  real<lower=0>                  mu_a;       // mean threshold
  real<lower=0>                  sigma_a;    // SD threshold
  real<lower=0>                  mu_ndt;     // mean non-decision time
  real<lower=0>                  sigma_ndt;  // SD non-decision time
  
  real<lower=0, upper=1>         mu_z0;      // mean starting point
  real<lower=0, upper=1>         sigma_z0;   // SD starting point
  real                           mu_bz;      // mean starting point beta
  real<lower=0>                  sigma_bz;   // SD starting point beta
   
  real                           mu_v0;      // mean drift intercept
  real<lower=0>                  sigma_v0;   // SD drift intercept
  real                           mu_b1v;     // mean drift b1
  real<lower=0>                  sigma_b1v;  // SD drift b1
  
  real<lower=0, upper=1>         mu_g;       // mean g parameter
  real<lower=0>                  sigma_g;    // SD g parameter
  
  //individual parameter  
  real<lower=0>                  a[S];       // individual threshold
  real<lower=0>                  ndt[S];     // individual non-decision time
  real<lower=0, upper=1>         z0[S];      // individual starting point
  real                           bz[S];      // individual beta starting point
  real                           v0[S];      // individual drift intercept
  real                           b1v[S];     // individual drift beta 1
  real<lower=0, upper=1>         g[S];       // individual g
}

transformed parameters {
  vector[T] delta;
  vector[T] beta;
  vector[T] X_1;
  vector[T] X_2;
  real      I;
  
  // Internal representation updating
  I = d1[1];
  X_1[1] = d1[1];
  X_2[1] = d2[1];

  for (i in 2:T) {
    if (sub[i] != sub[i-1]) {
      I = d1[i];
      X_1[i] = d1[i];
      X_2[i] = d2[i];
    } else {
      I = g[sub[i]]*I + (1-g[sub[i]])*d1[i];
      X_1[i] = I;
      I = g[sub[i]]*I + (1-g[sub[i]])*d2[i];
      X_2[i] = I;
    }
  }
  
  // Compute  trial-by-trial delta
  for (i in 1:T){
    delta[i] = v0[sub[i]] + b1v[sub[i]] * ((X_1[i] - X_2[i])/(X_1[i] + X_2[i]));
  }
  
    // Compute trial-by-trial beta
  for (i in 1:T){
    beta[i] = z0[sub[i]] + bz[sub[i]] * X_1[i];
  }
}

model {
  // hyper priors
  mu_a      ~ gamma(2, 2);
  sigma_a   ~ gamma(1, 5);
  mu_ndt    ~ gamma(3, 15);
  sigma_ndt ~ gamma(1, 5);

  mu_z0     ~ beta(5, 5);
  sigma_z0  ~ gamma(1, 5);
  mu_bz     ~ normal(0, 5);
  sigma_bz  ~ gamma(1, 5);

  mu_v0     ~ normal(0, 5);
  sigma_v0  ~ gamma(1, 5);
  mu_b1v    ~ normal(0, 5);
  sigma_b1v ~ gamma(1, 5);
  
  mu_g      ~ beta(1, 1);
  sigma_g   ~ gamma(1, 5);

  // individual priors
  for (i in 1:S) {
    a[i]     ~ normal(mu_a, sigma_a) T[0,];
    ndt[i]   ~ normal(mu_ndt, sigma_ndt) T[0,];
    z0[i]    ~ normal(mu_z0, sigma_z0);
    bz[i]   ~ normal(mu_bz, sigma_bz);
    v0[i]    ~ normal(mu_v0, sigma_v0);
    b1v[i]   ~ normal(mu_b1v, sigma_b1v);
    g[i]     ~ normal(mu_g, sigma_g)T[0,1];
  }

  for (i in 1:T) {
    if (resp[i] == 1) {
      rt[i] ~ wiener(a[sub[i]], ndt[sub[i]], beta[i], delta[i]);
    } else {
        rt[i] ~ wiener(a[sub[i]], ndt[sub[i]], 1-beta[i], -delta[i]);
    }
  }
  
}

generated quantities {
  vector[T] log_lik;
  for (i in 1:T) {
    if(resp[i]==1) {
      log_lik[i] = wiener_lpdf(rt[i] | a[sub[i]], ndt[sub[i]], beta[i], delta[i]);
    } else {
      log_lik[i] = wiener_lpdf(rt[i] | a[sub[i]], ndt[sub[i]], 1-beta[i], -delta[i]);
    }
  }
}
