// Hamstr with arbitrary breaks at each level
// 25.03.20203 Andrew Dolman
data{
  // age control points
  int<lower=0> N;
  vector[N] depth;
  vector[N] obs_age;
  vector[N] obs_err;
  real min_age;

  // resolution of age-depth model
  int<lower=1> n_lvls; // number of hierarchical levels, not including overall mean
  int<lower=0> K_fine; // number of highest resolution sections
  int<lower=0> K_tot; // total no of gamma parameters

  array[K_tot-1] int parent1;
  vector[K_tot-1] wts1; // weights of parents

  array[K_tot-1] int parent2;
  vector[K_tot-1] wts2; // weights of parents

  // modelled depths
  vector[K_fine] c_depth_bottom;
  vector[K_fine] c_depth_top;
  real<lower = 0> delta_c; // width of each highest resolution section


  // hyperparameters for the gamma innovations

  // prior for the oversall mean accumulation rate
  real<lower = 0> acc_mean_prior;

  // shape of the gamma distributions
  real<lower = 0> acc_shape;
  real<lower = 0> multi_parent_adj;

  // scale the shape parameter to control the total variance of the alpha
  // innovations for the number of hierarchical levels
  int<lower=0, upper=1> scale_shape;

  // hyperparameters for prior distribution on memory strength (the AR1 coefficient)
  real<lower = 0> mem_mean;
  real<lower = 0> mem_strength;

  // observation error model parameters

  int<lower=0> nu; // degrees of freedom of t error distribution
  array[N] int which_c;

  int<lower=0, upper=1> scale_R; // scale the AR1 coefficient or not
  int<lower=0, upper=1> inflate_errors; // use error inflation model or not


  // hyperparameters for the error inflation model
  real<lower = 0> infl_shape_shape;
  real<lower = 0> infl_shape_mean;

  real<lower = 0> infl_sigma_sd;


  // Additional data for bioturbation model

  // use bioturbation model or not
  int<lower=0, upper=1> model_bioturbation;

  int<lower=0> I;
  array[I, N] int smooth_i;

  real<lower = 0> L_prior_mean;
  real<lower = 0> L_prior_shape;

  vector[N*model_bioturbation] n_ind;

  // Additional data for modelling displacement
  int<lower=0, upper=1> model_displacement;
  real<lower = 0> D_prior_mean;
  real<lower = 0> D_prior_shape;

  int<lower=0, upper=1> smooth_s;

  // Model hiatuses

  int<lower=0, upper=1> model_hiatus;
  real H_top;
  real H_bottom;

}
transformed data{

  real min_depth = min(c_depth_top);
  real max_depth = max(c_depth_bottom);
  real data_age_range = max(obs_age) - min(obs_age);

  // inverse scale of the prior on L
  real L_rate;

  int<lower = 0, upper = 1> sample_L;

  // transform mean and strength of memory beta distribution to alpha and beta
  // as used to parameterise beta dist in stan function
  real<lower=0> mem_alpha = mem_strength * mem_mean;
  real<lower=0> mem_beta = mem_strength * (1-mem_mean);

  // position of the first highest resolution innovation (alpha)
  int<lower = 1> first_K_fine = K_tot - K_fine+1;

  // scale shape
  real<lower = 1> acc_shape_adj;
  if (scale_shape == 1){
    acc_shape_adj = acc_shape * n_lvls / multi_parent_adj;
  } else{
    acc_shape_adj = acc_shape;
  }

  L_rate = L_prior_shape / L_prior_mean;

  // bioturbation depth can be fixed rather than sampled
  if (L_prior_shape == 0) {
    sample_L = 0;
  } else {
    sample_L = 1;
  }

}
parameters {
  // AR1 coeffiecient at 1 depth unit
  real<lower = 0, upper = 1> R;

  // the hierarchical gamma innovations in one long vector that will be indexed
  vector<lower = 0>[K_tot] alpha;

  // the age at the first modelled depth
  real<lower = min_age> age0;

  // the measurement error inflation factors
  // these have length 0 if inflate_errors == 0 meaning that the parameters are
  // in scope, so the model runs, but are zero length so nothing is sampled
  array[inflate_errors] real<lower = 0> infl_mean;
  array[inflate_errors] real<lower = 0> infl_shape_1;
  vector<lower = 0>[inflate_errors ? N : 0] infl;

  array[model_bioturbation * sample_L] real<lower = 0> L;

  vector<lower = 0>[model_bioturbation ? N : 0] bt_error;

  array[model_displacement] real<lower = 0> D;

  array[model_hiatus] real<lower = H_top, upper = H_bottom> H_depth;
  array[model_hiatus] real<lower = 0, upper = data_age_range> H_length;
}

transformed parameters{

  // the AR1 coefficient scaled for the thickness of the modelled sediment sections
  real<lower = 0, upper = 1> w;

  // the highest resolution AR1 correlated innovations
  vector[K_fine] x;

  // the modelled ages
  vector[K_fine+1] c_ages;

  // the modelled ages interpolated to the positions of the data
  vector[N] Mod_age;

  // the inflated observation errors
  array[inflate_errors] real<lower = 0> infl_shape;

  vector[(inflate_errors == 1 || model_displacement == 1) ? N : 0] obs_err_infl;

  // latent bioturbation corrected age
  vector[model_bioturbation ? N : 0] bt_age;
  vector[(model_bioturbation == 1 || model_displacement == 1) ? N : 0] smooth_x;

  // age heterogeneity due to bioturbation at locations of observed ages
  vector[model_bioturbation ? N : 0] age_het;

  vector<lower = 0>[model_displacement ? N : 0] disp_yrs;

  if (scale_R == 1){
    w = R^(delta_c);
  } else {
    w = R;
  }

  // only the "fine" alphas
  // the first innovation
  x[1] = alpha[first_K_fine];

  // the remaining innovations with the AR1 parameter applied
  for(i in 2:K_fine){
    x[i] = w*x[i-1] + (1-w)*alpha[i + first_K_fine -1];
  }


  // the cumulative sum of the highest resolution innovations
  c_ages[1] = age0;
  c_ages[2:(K_fine+1)] = age0 + cumulative_sum(x * delta_c);


  // hiatus vector
  if (model_hiatus == 1){
    for (i in 2:(K_fine+1)){
      if (H_depth[1] < c_depth_top[i-1]) c_ages[i] = c_ages[i] + H_length[1];
      }
      }


  // age model interpolated to the positions of the observations
  Mod_age = c_ages[which_c] + x[which_c] .* (depth - c_depth_top[which_c]);



  if (model_bioturbation == 1 || model_displacement == 1){
    if (smooth_s == 1){
      for (n in 1:N){
        smooth_x[n] =  mean(x[smooth_i[,n]]);
      }
    } else {
      for (n in 1:N){
        smooth_x[n] =  x[which_c[n]];
      }
    }
  }


  if (model_bioturbation == 1){

    if (sample_L == 1){
      age_het = L[1] * smooth_x;
      // the modelled (shifted) gamma distributed bioturbation error
      // subtract the age_het to centre the error around the obs_age
      bt_age = obs_age + bt_error - age_het;

    }
    if (sample_L == 0){
      age_het = L_prior_mean * smooth_x;
      // the modelled (shifted) gamma distributed bioturbation error
      // subtract the age_het to centre the error around the obs_age
      bt_age = obs_age + bt_error - age_het;


    }
  }

  if (inflate_errors == 1 && model_displacement == 1){
    disp_yrs = D[1] * smooth_x;
    obs_err_infl = sqrt((obs_err .* obs_err) + (infl .* infl) + (disp_yrs .* disp_yrs));
    infl_shape[1] = infl_shape_1[1] + 1;
  } else if (inflate_errors == 1 && model_displacement == 0){
    for (n in 1:N)
    obs_err_infl[n] = sqrt((obs_err[n])^2 + (infl[n])^2);
    infl_shape[1] = infl_shape_1[1] + 1;
  } else if (inflate_errors == 0 && model_displacement == 1){
    disp_yrs = D[1] * smooth_x;
    obs_err_infl = sqrt((obs_err .* obs_err) + (disp_yrs .* disp_yrs));
  }


}
model {

  // the overall mean accumulation rate
  // weak half normal prior
  alpha[1] ~ normal(0, 10*acc_mean_prior);

  // the gamma distributed innovations

  {
  vector[K_tot-1] parent_mean;
  parent_mean = alpha[parent1] .* wts1 + alpha[parent2] .* wts2;

  // prior parameterised by use set shape and the value of it's parent section
  alpha[2:K_tot] ~ gamma(acc_shape_adj, acc_shape_adj ./ parent_mean);
  }

  // the memory parameters
  R ~ beta(mem_alpha, mem_beta);

  // the observation error inflation model

  if (inflate_errors){
    infl_shape ~ gamma(infl_shape_shape, infl_shape_shape / infl_shape_mean);
    infl_mean ~ normal(0, infl_sigma_sd);
    infl ~ gamma(infl_shape[1], infl_shape[1] / infl_mean[1]);
  }

  // bioturbation error model

  // parameters that are zero length do not get sampled
  L ~ gamma(L_prior_shape, L_rate);
  D ~ gamma(D_prior_shape, D_prior_shape / D_prior_mean);

  // additional error in ages due to age-heterogeneity
  bt_error ~ gamma(n_ind, n_ind ./ age_het);

  if (model_bioturbation == 1){
    // bioturbation depth

    if (inflate_errors == 1 || model_displacement == 1){
      // the Likelihood of the data given the model
      bt_age ~ student_t(nu, Mod_age, obs_err_infl);
    } else {
      bt_age ~ student_t(nu, Mod_age, obs_err);
    }

  } else {
    if (inflate_errors == 1 || model_displacement == 1){
      // the Likelihood of the data given the model
      obs_age ~ student_t(nu, Mod_age, obs_err_infl);
    } else {
      obs_age ~ student_t(nu, Mod_age, obs_err);
    }
  }
}
