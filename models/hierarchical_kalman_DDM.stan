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
  real<lower=0>             sd_a;      // SD threshold
  real                      mu_ndt;    // mean non-decision time
  real<lower=0>             sd_ndt;    // SD non-decision time
  real                      mu_z0;     // mean starting point intercept
  real<lower=0>             sd_z0;     // SD starting point intercept
  // real                      mu_bz;     // mean starting point beta
  // real<lower=0>             sd_bz;     // SD starting point beta
  real                      mu_v0;     // mean drift intercept
  real<lower=0>             sd_v0;     // SD drift intercept
  real                      mu_b1v;    // mean drift beta 1
  real<lower=0>             sd_b1v;    // mean drift beta 1
  // real                      mu_initVar;
  // real<lower=0>             sd_initVar;
  real                      mu_noise;
  real<lower=0>             sd_noise;
  real                      mu_q;
  real<lower=0>             sd_q;
  
  real                      z_a[S];   
  real                      z_ndt[S];  
  real                      z_z0[S];   
  // real                      z_bz[S];  
  real                      z_v0[S];   
  real                      z_b1v[S];  
  real                      z_noise[S]; 
  real                      z_q[S];
  // real                      z_initVar[S];
  
  real<lower=0>           percept_noise_1[T];
  real<lower=0>           percept_noise_2[T];

}

transformed parameters {
  // drift and bias
  vector[T]               delta;
  // noisy percept
  vector[T]               xd1;
  vector[T]               xd2;
  // final stimulus estimate
  real                    I;
  vector[T]               I_1;
  vector[T]               I_2;
  // kalman filter parameters
  real<lower=0, upper=1>  gain;
  real<lower=0>           initVar;
  real<lower=0>           priorVar;

  //individual parameter   
  real<lower=0>             a[S];      // individual threshold
  real<lower=0>             ndt[S];    // individual non-decision time
  real<lower=0, upper=1>    z0[S];     // individual starting point intercept
  // real                      bz[S];     // individual starting point intercept
  real                      v0[S];     // individual drift intercept
  real                      b1v[S];    // individual drift beta 1
  real                      q[S];
  real                      noise[S];

  for (s in 1:S){
    a[s]                  = exp(mu_a + z_a[s] * sd_a);
    ndt[s]                = exp(mu_ndt + z_ndt[s] * sd_ndt);
    z0[s]                 = Phi(mu_z0 + z_z0[s] * sd_z0);
    // bz[s]                 = mu_bz + z_bz[s] * sd_bz;
    v0[s]                 = mu_v0 + z_v0[s] * sd_v0;
    b1v[s]                = mu_b1v + z_b1v[s] * sd_b1v;
    q[s]                  = exp(mu_q + z_q[s] * sd_q);
    noise[s]              = exp(mu_noise + z_noise[s] * sd_noise);
    }

  // noisy perception
  xd1 = d1 + percept_noise_1[1];
  xd2 = d2 + percept_noise_2[1];
  
  // intial value
  initVar = 0;
  I = 0;
  // first trial
  // update kalman gain
  gain = (initVar + q[1]) / (initVar + q[1] + noise[1]);
  // update posterior variance
  priorVar = gain * noise[1];
  // calculate percept
  I = (1-gain) * I + gain * xd1[1];
  I_1[1] = I;
  
  gain = (initVar + q[1]) / (initVar + q[1] + noise[1]);
  priorVar = gain * noise[1];
  I = (1-gain) * I + gain * xd2[1];
  I_2[1] = I;
  delta[1] = v0[1] + b1v[1] * (I_1[1] - I_2[1]);
  
  // iterate over all trials
  for (i in 2:T) {
    if (sub[i] != sub[i-1]) {
      initVar = 0;
      I = 0;
      // stimulus 1
      gain = (initVar + q[sub[i]]) / (initVar + q[sub[i]] + noise[sub[i]]);
      priorVar = gain * noise[sub[i]];
      I = (1-gain) * I + gain * xd1[i];
      I_1[i] = I;
      // stimulus 2
      gain = (initVar + q[sub[i]]) / (initVar + q[sub[i]] + noise[sub[i]]);
      priorVar = gain * noise[sub[i]];
      I = (1-gain) * I + gain * xd2[i];
      I_2[i] = I;
    } else {
      // stimulus 1
      gain = (priorVar + q[sub[i]]) / (priorVar + q[sub[i]] + noise[sub[i]]);
      priorVar = gain * noise[sub[i]]; 
      I = (1-gain) * I + gain * xd1[i];
      I_1[i] = I;
      // stimulus 2
      gain = (priorVar + q[sub[i]]) / (priorVar + q[sub[i]] + noise[sub[i]]);
      priorVar = gain * 2;
      I = (1-gain) * I + gain * xd2[i];
      I_2[i] = I;
    }
    // trial-by-trial drift rate
    delta[i] = v0[sub[i]] + b1v[sub[i]] * (I_1[i] - I_2[i]);
  }
}

model {
  // priors
  mu_a     ~ normal(0,3);
  sd_a     ~ normal(0,2);
  mu_ndt   ~ normal(-1,1);
  sd_ndt   ~ normal(0,2);
  mu_z0    ~ normal(0,1);
  sd_z0    ~ normal(0,2);
  mu_v0    ~ normal(0,3);
  sd_v0    ~ normal(0,2);
  mu_b1v   ~ normal(0,3);
  sd_b1v   ~ normal(0,2);
  mu_q     ~ normal(2,1);
  sd_q     ~ normal(0,2);
  // mu_noise ~ normal(0,3);
  // sd_noise ~ normal(0,2);

  z_a      ~ normal(0,1);
  z_ndt    ~ normal(0,1); 
  z_z0     ~ normal(0,1);
  z_v0     ~ normal(0,1);
  z_b1v    ~ normal(0,1);
  z_q      ~ normal(0,1);
  // z_noise  ~ normal(0,1);


  for (i in 1:T) {
    percept_noise_1[i] ~ normal(0,2);
    percept_noise_2[i] ~ normal(0,2);
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
