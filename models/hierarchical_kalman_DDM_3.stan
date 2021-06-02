data {
  int<lower=1>                   T;         // number of trials 
  int<lower=1>                   S;         // number of subjects
  int<lower=1>                   sub[T];    // subject number
  int<lower=1, upper=2>          resp[T];   // response
  vector<lower=0>[T]             rt;        // response time
  vector[T]                      d1;        // stimulus 1
  vector[T]                      d2;        // stimulus 2
}

parameters {
  //hyper parameter
  real                      mu_a;      // mean threshold
  real<lower=0>             sigma_a;      // sigma threshold
  real                      mu_ndt;    // mean non-decision time
  real<lower=0>             sigma_ndt;    // sigma non-decision time
  real                      mu_z0;     // mean starting point intercept
  real<lower=0>             sigma_z0;     // sigma starting point intercept
  // real                      mu_bz;     // mean starting point beta
  // real<lower=0>             sigma_bz;     // sigma starting point beta
  real                      mu_v0;     // mean drift intercept
  real<lower=0>             sigma_v0;     // sigma drift intercept
  real                      mu_b1v;    // mean drift beta 1
  real<lower=0>             sigma_b1v;    // mean drift beta 1
  // real                      mu_initVar;
  // real<lower=0>             sigma_initVar;
  real                      mu_noise;
  real<lower=0>             sigma_noise;
  real                      mu_q;
  real<lower=0>             sigma_q;
  
  
  //individual parameter  
  real<lower=0>                  a[S];       // individual threshold
  real<lower=0>                  ndt[S];     // individual non-decision time
  real<lower=0, upper=1>         z0[S];      // individual starting point intercept
  real                           v0[S];      // individual drift intercept
  real                           b1v[S];     // individual drift beta 1
  real                           q[S];       // individual g
  real<lower=0>                  noise[S];   // individual perpect noise

}

transformed parameters {
  // drift and bias
  vector[T]               delta;
  // noisy percept
  // vector[T]               xd1;
  // vector[T]               xd2;
  // final stimulus estimate
  real                    I;
  vector[T]               I_1;
  vector[T]               I_2;
  // kalman filter parameters
  real<lower=0, upper=1>  gain;
  real<lower=0>           priorVar;


  // noisy perception
  // xd1 = d1 + percept_noise_1[1];
  // xd2 = d2 + percept_noise_2[1];
  
  // intial value
  priorVar = 0;
  I = 0;
  // first trial
  // update kalman gain
  gain = (priorVar + q[1]) / (priorVar + q[1] + noise[1]);
  // update posterior variance
  priorVar = gain * noise[1];
  // calculate percept
  I = (1-gain) * I + gain * d1[1];
  I_1[1] = I;
  
  gain = (priorVar + q[1]) / (priorVar + q[1] + noise[1]);
  priorVar = gain * noise[1];
  I = (1-gain) * I + gain * d2[1];
  I_2[1] = I;
  delta[1] = v0[1] + b1v[1] * (I_1[1] - I_2[1]);
  
  
  // iterate over all trials
  for (i in 2:T) {
    
    if (sub[i] != sub[i-1]) {
      priorVar = 0;
      I = 0;
      }
    
    // stimulus 1
    gain = (priorVar + q[sub[i]]) / (priorVar + q[sub[i]] + noise[sub[i]]);
    priorVar = gain * noise[sub[i]]; 
    I = (1-gain) * I + gain * d1[i];
    I_1[i] = I;
    // stimulus 2
    gain = (priorVar + q[sub[i]]) / (priorVar + q[sub[i]] + noise[sub[i]]);
    priorVar = gain * noise[sub[i]];
    I = (1-gain) * I + gain * d2[i];
    I_2[i] = I;

    // trial-by-trial drift rate
    delta[i] = v0[sub[i]] + b1v[sub[i]] * (I_1[i] - I_2[i]);
  }
}

model {
  // priors
  mu_a        ~ gamma(2, 2);
  sigma_a     ~ gamma(1, 5);
  mu_ndt      ~ gamma(3, 15);
  sigma_ndt   ~ gamma(1, 5);

  mu_z0       ~ beta(5, 5);
  sigma_z0    ~ gamma(1, 5);
  mu_v0       ~ normal(0, 5);
  sigma_v0    ~ gamma(1, 5);
  mu_b1v      ~ normal(0, 5);
  sigma_b1v   ~ gamma(1, 5);
  mu_q        ~ normal(0,10);
  sigma_q     ~ gamma(1, 5);
  mu_noise    ~ gamma(2, 1);
  sigma_noise ~ gamma(1, 5);

  // individual priors
  for (i in 1:S) {
    a[i]     ~ normal(mu_a, sigma_a) T[0,];
    ndt[i]   ~ normal(mu_ndt, sigma_ndt) T[0,];
    z0[i]    ~ normal(mu_z0, sigma_z0);
    v0[i]    ~ normal(mu_v0, sigma_v0);
    b1v[i]   ~ normal(mu_b1v, sigma_b1v);
    q[i]     ~ normal(mu_q, sigma_q);
    noise[i] ~ normal(mu_noise, sigma_noise) T[0,];
  }


  for (i in 1:T) {
    // percept_noise_1[i] ~ normal(0,2);
    // percept_noise_2[i] ~ normal(0,2);
    if(resp[i]==1) {
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
