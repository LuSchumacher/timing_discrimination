data {
  int<lower=1>            T;         // number of trials 
  int<lower=1>            sub[T];    // subject number
  int<lower=1, upper=2>   resp[T];   // response
  vector<lower=0>[T]      rt;        // response time
  vector[T]               d1;        // stimulus 1
  vector[T]               d2;        // stimulus 2
}

parameters {
  real<lower=0, upper=5>  a;
  real<lower=0, upper=1>  ndt;
  real<lower=0, upper=1>  z0;
  real                    v0;
  real                    b1v;
  
  real<lower=0>           percept_noise_1[T];
  real<lower=0>           percept_noise_2[T];
  real<lower=0>           kalman_prior_var_start;
  real<lower=0>           noise_sd;
  real                    kalman_q;
}

transformed parameters {
  // drift and bias
  vector[T]               delta;
  // noisy percept
  vector[T]               xd1;
  vector[T]               xd2;
  // final stimulus estimate
  vector[T]               I_1;
  vector[T]               I_2;
  // kalman filter parameters
  vector[T]               kalman_k;
  vector[T]               kalman_prior_var;

  // starting values
  xd1[1] = d1[1] + percept_noise_1[1]; // df[1] + normal(0,SD) --> SD is a parameter that should be estimated
  xd2[1] = d2[1] + percept_noise_2[1];
  kalman_prior_var[1] = kalman_prior_var_start;

  
  for (i in 2:T) {
      xd1[i] = d1[i] + percept_noise_1[i];
      xd2[i] = d2[i] + percept_noise_2[i];
    if (sub[i] != sub[i-1]) {
      kalman_prior_var[i] = kalman_prior_var_start;
    } else {
      kalman_k[i] = (kalman_prior_var[i-1] + kalman_q) / (kalman_prior_var[i-1] + kalman_q + noise_sd);
      kalman_prior_var[i] = kalman_k[i] * noise_sd;
      I_1[i] = (1-kalman_k[i]) * I_1[i-1] + kalman_k[i] * xd1[i];
      kalman_k[i] = (kalman_prior_var[i-1] + kalman_q) / (kalman_prior_var[i-1] + kalman_q + noise_sd);
      kalman_prior_var[i] = kalman_k[i] * noise_sd;
      I_2[i] = (1-kalman_k[i]) * I_2[i-1] + kalman_k[i] * xd1[i];
    }
  }
  
  delta = v0 + b1v * (I_1 - I_2);
}

model {
  a ~ gamma(1,1);
  z0 ~ beta(10,10);
  v0 ~ normal(0,5);
  b1v ~ normal(0,1);
  ndt ~ gamma(1,1);
  
  noise_sd ~ gamma(1,0.2);
  //kalman_prior_var_start ~ gamma(1,0.2);
  kalman_q ~ uniform(-1,1);
  
  for (i in 1:T) {
    percept_noise_1[i] ~ normal(0,noise_sd);
    percept_noise_2[i] ~ normal(0,noise_sd);
    if(resp[i]==1) {
      rt[i] ~ wiener(a, ndt, z0, delta[i]);
    } else {
      rt[i] ~ wiener(a, ndt, z0, -delta[i]);
    }
  }
}

generated quantities {
  vector[T] log_lik;
  for (i in 1:T) {
    if(resp[i]==1) {
      log_lik[i] = wiener_lpdf(rt[i] | a, ndt, z0, delta[i]);
    } else {
      log_lik[i] = wiener_lpdf(rt[i] | a, ndt, z0, -delta[i]);
    }
  }
}
