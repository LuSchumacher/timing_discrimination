data {
  int<lower=1>                   T;         // number of trials 
  int<lower=1>                   S;         // number of subjects
  int<lower=1>                   sub[T];    // subject number
  int<lower=1, upper=2>          resp[T];   // response
  real                           min_rt;    // smallest rt in the data
  vector<lower=0>[T]             rt;        // response time
  //vector[T]                      d1;        // stimulus 1
  vector[T]                      d2;        // stimulus 2
}

parameters {
  //hyper parameter
  real<lower=0>                  mu_a;       // mean threshold
  real<lower=0>                  sigma_a;    // SD threshold
  real<lower=0>                  mu_ndt;     // mean non-decision time
  real<lower=0>                  sigma_ndt;  // SD non-decision time
  
  real<lower=0, upper=1>         mu_z0;      // mean starting point intercept
  real<lower=0, upper=1>         sigma_z0;   // SD starting point intercept
  
  real                           mu_v0;      // mean drift intercept
  real<lower=0>                  sigma_v0;   // SD drift intercept
  real                           mu_b2v;     // mean drift beta 2
  real<lower=0>                  sigma_b2v;  // mean drift beta 2
  
  //individual parameter  
  real<lower=0>                  a[S];       // individual threshold
  real<lower=0>                  ndt[S];     // individual non-decision time
  real<lower=0, upper=1>         z0[S];      // individual starting point intercept
  real                           v0[S];      // individual drift intercept
  real                           b2v[S];     // individual drift beta 2
  
}

transformed parameters {
  vector[T] delta;
  
  // Compute  trial-by-trial parameters
  for (i in 1:T){
    delta[i] = v0[sub[i]] + b2v[sub[i]] * d2[i];
  }

}

model {
  
  // hyper priors
  mu_a      ~ gamma(4, 3);
  sigma_a   ~ gamma(1.5, 3);
  mu_ndt    ~ normal(3, 0.5) T[0,];
  sigma_ndt ~ gamma(1.5, 1);

  mu_z0     ~ beta(5, 5);
  sigma_z0  ~ gamma(1, 5);

  mu_v0     ~ normal(0, 2);
  sigma_v0  ~ gamma(1, 2);
  mu_b2v    ~ normal(0, 2);
  sigma_b2v ~ gamma(1, 2);
  

  // individual priors
  for (i in 1:S) {
    a[i]     ~ normal(mu_a, sigma_a) T[0,];
    ndt[i]   ~ normal(mu_ndt, sigma_ndt) T[0,];
    z0[i]    ~ normal(mu_z0, sigma_z0);
    v0[i]    ~ normal(mu_v0, sigma_v0);
    b2v[i]   ~ normal(mu_b2v, sigma_b2v);
  }

  for (i in 1:T) {
    if (resp[i] == 1) {
      rt[i] ~ wiener(a[sub[i]], ndt[sub[i]], z0[sub[i]], delta[i]);
    } else {
        rt[i] ~ wiener(a[sub[i]], ndt[sub[i]], 1-z0[sub[i]], -delta[i]);
    }
  }
   
}

generated quantities {
  vector[T] log_lik;
  for (i in 1:T) {
    if(resp[i]==1) {
      log_lik[i] = wiener_lpdf(rt[i] | a[sub[i]], ndt[sub[i]], z0[sub[i]], delta[i]);
    } else {
      log_lik[i] = wiener_lpdf(rt[i] | a[sub[i]], ndt[sub[i]], 1-z0[sub[i]], -delta[i]);
    }
  }
}
