// Bacon as a Stan model
// 25.03.2017 Andrew Dolman
data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> nu;
  vector[N] depth;
  vector[N] obs_age;
  vector[N] obs_err;
  vector[K] c_depth_bottom;
  vector[K] c_depth_top;
  int which_c[N];
  real<lower = 0> delta_c;

  // hyperparameters
  real<lower = 0> acc_mean;
  real<lower = 0> acc_var;

  real<lower = 0> mem_mean;
  real<lower = 0> mem_strength;
}
parameters {
  real<lower = 0, upper = 1> R;
  vector<lower = 0>[K] alpha;
  real<lower = 0> age0;
}
transformed parameters{
  // transformed hyperparameters
  real<lower=0> acc_alpha;
  real<lower=0> acc_beta;

  real<lower=0> mem_alpha;
  real<lower=0> mem_beta;

  real<lower = 0, upper = 1> w;

  vector[K] x;
  vector[K+1] c_ages;
  vector[N] Mod_age;

  acc_alpha = acc_mean^2 / acc_var;
  acc_beta = acc_mean / acc_var;

  mem_alpha = mem_strength * mem_mean;
  mem_beta = mem_strength * (1-mem_mean);

  w = R^(delta_c/1);

  x[1] = alpha[1];

  //tail(x, N - 1) = w*head(x, N - 1) + (1-w)*alpha;
  for(i in 2:K){
    x[i] = w*x[i-1] + (1-w)*alpha[i];
    }


  // Get the Mod_ages
  c_ages[1] = 0;
  c_ages[2:(K+1)] = cumulative_sum(alpha * delta_c);
  Mod_age = age0 + c_ages[which_c] + alpha[which_c] .* (depth - c_depth_top[which_c]);
}
model {
  age0 ~ normal(obs_age[1], 100);
  R ~ beta(mem_alpha, mem_beta);
  alpha ~ gamma(acc_alpha, acc_beta);
  //obs_age ~ normal(Mod_age, obs_err);
  obs_age ~ student_t(nu, Mod_age, obs_err);
  }
