// Measurement model in the appendix
functions {
  row_vector center(row_vector a) {
    return a - mean(a);
  }
  
  vector center_col(vector a) {
    return a - mean(a);
  }
}

data {
  int<lower = 1> N;
  int<lower = 1> M;
  int<lower = 2> K;
  int<lower = 1> D;
  int<lower = 1, upper = K> y[N, M];
  int<lower = 1, upper = D> d[M];
  vector[M] s;
}

parameters {
  matrix[K + 1, M] z_m;
  matrix[K, N] z_n;
  matrix[D, N] z_theta;
  cholesky_factor_corr[K + 1] L_omega_m;
  cholesky_factor_corr[K] L_omega_n;
  cholesky_factor_corr[D] L_omega_theta;
  vector<lower = 0>[K + 1] tau_m;
  vector<lower = 0>[K] tau_n;
}

transformed parameters {
  matrix[N, D] theta;
  vector[M] a;
  matrix[M, K] c;
  matrix[N, K] eta;
  vector[N] log_lik;
  
  {
    matrix[M, K + 1] x_m = (diag_pre_multiply(tau_m, L_omega_m) * z_m)';
    matrix[N, K] x_n = (diag_pre_multiply(tau_n, L_omega_n) * z_n)';
    
    theta = (L_omega_theta * z_theta)';
    
    for (m in 1:M)
      c[m] = center(x_m[m, 1:K]);
    a = exp(x_m[, K + 1]) .* s;
    
    for (n in 1:N)
      eta[n, ] = center(x_n[n, ]);
    for (k in 1:K)
      eta[, k] = center_col(eta[, k]);
    
    for (n in 1:N) {
      log_lik[n] = 0;
      for (m in 1:M) {
        vector[K] p;
        for (k in 1:K)
          p[k] = (k - (K + 1) * 0.5) * a[m] * theta[n, d[m]] + c[m, k] + eta[n, k];
        log_lik[n] += categorical_logit_lpmf(y[n, m] | p);
      }
    }
  }
}

model {
  to_vector(z_m) ~ std_normal();
  to_vector(z_n) ~ std_normal();
  to_vector(z_theta) ~ std_normal();
  
  tau_m ~ cauchy(0, 1);
  tau_n ~ cauchy(0, 1);
  L_omega_m ~ lkj_corr_cholesky(2);
  L_omega_n ~ lkj_corr_cholesky(2);
  L_omega_theta ~ lkj_corr_cholesky(2);
  
  target += sum(log_lik);
}

generated quantities {
  matrix[K + 1, K + 1] omega_m = multiply_lower_tri_self_transpose(L_omega_m);
  matrix[K, K] omega_n = multiply_lower_tri_self_transpose(L_omega_n);
  matrix[D, D] omega_theta = multiply_lower_tri_self_transpose(L_omega_theta);
  matrix[K + 1, K + 1] sigma_m = quad_form_diag(omega_m, tau_m);
  matrix[K, K] sigma_n = quad_form_diag(omega_n, tau_n);
  vector[M] v_theta[N];
  vector[M] v_c[N];
  vector[M] v_eta[N];
  
  for (n in 1:N)
    for (m in 1:M) {
      int k = y[n, m];
      v_theta[n, m] = (k - (K + 1) * 0.5) * a[m] * theta[n, d[m]];
      v_c[n, m] = c[m, k];
      v_eta[n, m] = eta[n, k];
    }
}
