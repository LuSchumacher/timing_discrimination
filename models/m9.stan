data {
  int                       T;         // number of trials 
  int                       S;         // number of subjects
  int                       sub[T];    // subject number
  int                       resp[T];   // response
  vector[T]                 rt;        // response time
  vector[T]                 d1;        // stimulus 1
  vector[T]                 d2;        // stimulus 2
}

parameters {
  // hyper parameter
  real                      mu_a;      // mean threshold
  real<lower=0>             sigma_a;   // sd threshold
  real                      mu_ndt;    // mean non-decision time
  real<lower=0>             sigma_ndt; // sd non-decision time

  real                      mu_z0;     // mean starting point intercept
  real<lower=0>             sigma_z0;  // sd starting point intercept
  real                      mu_bz;     // mean starting point beta
  real<lower=0>             sigma_bz;  // sd starting point beta
  
  real                      mu_v0;     // mean drift intercept
  real<lower=0>             sigma_v0;  // sd drift intercept
  real                      mu_b1v;    // mean drift beta 1
  real<lower=0>             sigma_b1v; // sd drift beta 1
  real                      mu_b2v;    // mean drift beta 2
  real<lower=0>             sigma_b2v; // sd drift beta 2
  
  real                      mu_g;      // mean IRM updating param
  real<lower=0>             sigma_g;   // sd IRM updating param
  
  real                      z_a[S];   
  real                      z_ndt[S];  
  real                      z_z0[S];   
  real                      z_bz[S];  
  real                      z_v0[S];   
  real                      z_b1v[S];  
  real                      z_b2v[S];  
  real                      z_g[S];  
}

transformed parameters {
  vector[T]                 delta;    // trial-by-trial drift
  vector[T]                 beta;     // trial-by-trial threshold
  vector[T]                 I;        // internal reference stimulus 1
  
  // individual parameter   
  real<lower=0>             a[S];     // individual threshold
  real<lower=0>             ndt[S];   // individual non-decision time
  real<lower=0, upper=1>    z0[S];    // individual starting point intercept
  real                      bz[S];    // individual starting point intercept
  real                      v0[S];    // individual drift intercept
  real                      b1v[S];   // individual drift beta 1
  real                      b2v[S];   // individual drift beta 2
  real<lower=0, upper=1>    g[S];     // IRM updating param

  // non-centered parametrization
  for (s in 1:S){
    a[s]                  = exp(mu_a + z_a[s] * sigma_a);
    ndt[s]                = exp(mu_ndt + z_ndt[s] * sigma_ndt);
    z0[s]                 = Phi(mu_z0 + z_z0[s] * sigma_z0);
    bz[s]                 = mu_bz + z_bz[s] * sigma_bz;
    v0[s]                 = mu_v0 + z_v0[s] * sigma_v0;
    b1v[s]                = mu_b1v + z_b1v[s] * sigma_b1v;
    b2v[s]                = mu_b2v + z_b2v[s] * sigma_b2v;
    g[s]                  = Phi(mu_g + z_g[s] * sigma_g);
    }
  
  // internal representation updating
  I[1] = d1[1];

  for (i in 2:T) {
    if (sub[i] != sub[i-1]) {
      I[i] = d1[i];
    } else {
      I[i] = g[sub[i]]*I[i-1] + (1-g[sub[i]])*d1[i];
    }
  }
  
  // compute  trial-by-trial delta and beta
  for (i in 1:T){
    delta[i] = v0[sub[i]] + b1v[sub[i]] * I[i] + b2v[sub[i]] * d2[i];
    beta[i] = z0[sub[i]] + bz[sub[i]] * I[i];
  }
}

model {
  // priors
  mu_a       ~ normal(0,3);
  sigma_a    ~ normal(0,3);
  mu_ndt     ~ normal(-1,1);
  sigma_ndt  ~ normal(0,3);

  mu_z0      ~ normal(0,1);
  sigma_z0   ~ normal(0,1);
  mu_bz      ~ normal(0,3);
  sigma_bz   ~ normal(0,2);

  mu_v0      ~ normal(0,3);
  sigma_v0   ~ normal(0,2);
  mu_b1v     ~ normal(0,3);
  sigma_b1v  ~ normal(0,2);
  mu_b2v     ~ normal(0,3);
  sigma_b2v  ~ normal(0,2);
  
  mu_g       ~ normal(0,1);
  sigma_g    ~ normal(0,1);
  
  z_a        ~ normal(0,1);
  z_ndt      ~ normal(0,1); 
  z_z0       ~ normal(0,1);
  z_bz       ~ normal(0,1);
  z_v0       ~ normal(0,1);
  z_b1v      ~ normal(0,1);
  z_b2v      ~ normal(0,1);

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

