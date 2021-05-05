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
  //hyper parameter
  real                      mu_a;      // mean threshold
  real<lower=0>             sd_a;      // SD threshold

  real                      mu_ndt;    // mean non-decision time
  real<lower=0>             sd_ndt;    // SD non-decision time
  // real                      mu_ndt_var;// mean non-decision time trial-by-trial variability
  // real<lower=0>             sd_ndt_var;// SD non-decision time trial-by-trial variability
  
  real                      mu_z0;     // mean starting point intercept
  real<lower=0>             sd_z0;     // SD starting point intercept
  
  real                      mu_v0;     // mean drift intercept
  real<lower=0>             sd_v0;     // SD drift intercept
  real                      mu_b1v;    // mean drift beta 1
  real<lower=0>             sd_b1v;    // mean drift beta 1
  real                      mu_b2v;    // mean drift beta 2
  real<lower=0>             sd_b2v;    // mean drift beta 2
  
  real                      z_a[S];   
  real                      z_ndt[S];  
  real                      z_z0[S];   
  real                      z_v0[S];   
  real                      z_b1v[S];  
  real                      z_b2v[S];  
  
  // real                      mu_ndt_var;// mean non-decision time trial-by-trial variability
  // real<lower=0>             sd_ndt_var;// SD non-decision time trial-by-trial variability
  // real                      z_ndt_var;
  // real                      z_eps;
  real<lower=0>             ndt_var;
  real                      eps[T];

  
}

transformed parameters {
  real                      delta[T];  // trial-by-trial drift
  real<lower=0>             ndt_t[T];  // trial-by-trial  ndt

  //individual parameter   
  real<lower=0>             a[S];      // individual threshold
  real<lower=0>             ndt[S];    // individual non-decision time
  // real<lower=0>             ndt_var[S];// individual non-decision time trial-by-trial variability
  real<lower=0>             z0[S];     // individual starting point intercept
  real                      v0[S];     // individual drift intercept
  real                      b1v[S];    // individual drift beta 1
  real                      b2v[S];    // individual drift beta 2
  
  // real<lower=0>             ndt_var;
  // real                      eps[T];

  for (s in 1:S){
    a[s]                  = exp(mu_a + z_a[s] * sd_a);
    ndt[s]                = exp(mu_ndt + z_ndt[s] * sd_ndt);
    // ndt_var[s]            = exp(mu_ndt_var + z_ndt_var[s] * sd_ndt_var);
    z0[s]                 = Phi(mu_z0 + z_z0[s] * sd_z0);
    v0[s]                 = mu_v0 + z_v0[s] * sd_v0;
    b1v[s]                = mu_b1v + z_b1v[s] * sd_b1v;
    b2v[s]                = mu_b2v + z_b2v[s] * sd_b2v;
    }

  // ndt_var = exp(mu_ndt_var + z_ndt_var * sd_ndt_var);
  // eps = (0 + z_eps * ndt_var);

  // Compute  trial-by-trial delta
  for (i in 1:T){
    ndt_t[i] = ndt[sub[i]] + eps[i];
    delta[i] = v0[sub[i]] + b1v[sub[i]] * d1[i] + b2v[sub[i]] * d2[i];
  }

}

model {
  // priors
  mu_a     ~ normal(0,2);
  sd_a     ~ normal(0,2);
  
  mu_ndt   ~ normal(-1,2);
  sd_ndt   ~ normal(0,2);
  // mu_ndt_var   ~ normal(-2,1);
  // sd_ndt_var   ~ normal(0,2);
  ndt_var  ~ gamma(1,1);
  eps      ~ normal(0,ndt_var);

  mu_z0    ~ normal(0,1);
  sd_z0    ~ normal(0,1);

  mu_v0    ~ normal(0,3);
  sd_v0    ~ normal(0,2);
  mu_b1v   ~ normal(0,3);
  sd_b1v   ~ normal(0,2);
  mu_b2v   ~ normal(0,3);
  sd_b2v   ~ normal(0,2);

  z_a      ~ normal(0,1);
  z_ndt    ~ normal(0,1); 
  // z_ndt_var~ normal(0,1); 
  z_z0     ~ normal(0,1);
  z_v0     ~ normal(0,1);
  z_b1v    ~ normal(0,1);
  z_b2v    ~ normal(0,1);


  for (i in 1:T) {
    // ndt_t[i] ~ uniform(ndt[sub[i]]-ndt_var[sub[i]],ndt[sub[i]]+ndt_var[sub[i]]);
    if (resp[i] == 1) {
      rt[i] ~ wiener(a[sub[i]], ndt_t[i], z0[sub[i]], delta[i]);
    } else {
      rt[i] ~ wiener(a[sub[i]], ndt_t[i], 1-z0[sub[i]], -delta[i]);
    }
  }
}

generated quantities {
  vector[T] log_lik;
  for (i in 1:T) {
    if(resp[i]==1) {
      log_lik[i] = wiener_lpdf(rt[i] | a[sub[i]], ndt_t[i], z0[sub[i]], delta[i]);
    } else {
      log_lik[i] = wiener_lpdf(rt[i] | a[sub[i]], ndt_t[i], 1-z0[sub[i]], -delta[i]);
    }
  }
}
