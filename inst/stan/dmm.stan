// =============================================================================
// Dirichlet-Multinomial Mixture Model for eDNA Metabarcoding
// eDNAstructure package — Stan model
// =============================================================================
//
// GENERATIVE MODEL:
//
//   Community compositions:
//     pi_k  ~ Dirichlet(conc * 1_S)          for k = 1..K
//
//   Community membership (covariate-driven softmax regression):
//     Gamma_i  = softmax(beta_0k + beta_1k*cov1_i + ... + beta_Pk*covP_i)
//     Community K is the reference category (linear predictor = 0)
//
//   Observed counts (z marginalized out):
//     x_i ~ DirichletMultinomial(N_i, alpha * pi_k)   summed over k
//
//   Overdispersion:
//     alpha ~ gamma(alpha_shape, alpha_rate)
//       alpha >> 1 : near-multinomial (low overdispersion)
//       alpha ~  1 : high overdispersion (typical for eDNA)
//
// PARAMETERS TRACKED:
//   pi[K, S]              — posterior mean community compositions
//   beta[K-1, P+1]        — softmax regression coefficients (covariate effects)
//   alpha                 — global overdispersion scalar
//   community_probs[N, K] — posterior sample-level community membership probs
//   z_hat[N]              — MAP community assignment per sample
//   log_lik[N]            — per-sample log-likelihood (for LOO-CV)
//
// =============================================================================

data {
  int<lower=1> N;               // Number of samples
  int<lower=1> S;               // Number of taxa (or haplotypes)
  int<lower=2> K;               // Number of latent communities
  int<lower=0> P;               // Number of covariates (0 = intercept-only)

  array[N, S] int<lower=0> X;  // Count matrix [N x S]
  matrix[N, P] covariates;      // Covariate matrix [N x P], should be Z-scored

  real<lower=0> conc;           // Dirichlet concentration for pi priors
                                // conc < 1 = sparse communities (recommended for eDNA)
                                // conc = 1 = flat/uninformative
                                // conc > 1 = even compositions

  real<lower=0> alpha_shape;    // Gamma prior shape for overdispersion alpha
  real<lower=0> alpha_rate;     // Gamma prior rate  for overdispersion alpha
                                // Mean = shape/rate. Default: shape=5, rate=2 => mean=2.5
}

transformed data {
  // Per-sample read totals (N_i in the model)
  array[N] int<lower=0> Ni;
  for (i in 1:N) Ni[i] = sum(X[i]);

  // Precompute log multinomial coefficients for efficiency
  vector[N] log_mc;
  for (i in 1:N) {
    real lmc = lgamma(Ni[i] + 1.0);
    for (s in 1:S) lmc -= lgamma(X[i, s] + 1.0);
    log_mc[i] = lmc;
  }

  // Design matrix: intercept column prepended to covariates
  // When P=0 (intercept-only model), covariates is N x 0 matrix
  matrix[N, P + 1] X_design;
  if (P > 0) {
    X_design = append_col(rep_vector(1.0, N), covariates);
  } else {
    X_design = rep_matrix(1.0, N, 1);
  }
}

parameters {
  // K community compositions (K simplices over S taxa)
  array[K] simplex[S] pi;

  // Softmax regression coefficients [K-1 x (P+1)]
  // Community K is the reference: its linear predictor is fixed at 0
  matrix[K - 1, P + 1] beta;

  // Global overdispersion: single scalar shared across all communities
  real<lower=0> alpha;
}

transformed parameters {
  // Log-scale mixing weights [N x K] (community K = reference = 0)
  matrix[N, K] log_mixing_weights;
  {
    if (K > 1) {
      matrix[N, K - 1] eta = X_design * beta';
      log_mixing_weights = append_col(eta, rep_vector(0.0, N));
    } else {
      log_mixing_weights = rep_matrix(0.0, N, 1);
    }
  }
}

model {
  // ── Priors ────────────────────────────────────────────────────────────────
  for (k in 1:K)
    pi[k] ~ dirichlet(rep_vector(conc, S));

  to_vector(beta) ~ normal(0, 1.0);

  alpha ~ gamma(alpha_shape, alpha_rate);

  // ── Marginalized likelihood ───────────────────────────────────────────────
  for (i in 1:N) {
    vector[K] log_weights = log_softmax(log_mixing_weights[i]');
    vector[K] lp;
    for (k in 1:K) {
      vector[S] alpha_k = alpha * pi[k];
      real dm = log_mc[i] + lgamma(alpha) - lgamma(Ni[i] + alpha);
      for (s in 1:S)
        dm += lgamma(X[i, s] + alpha_k[s]) - lgamma(alpha_k[s]);
      lp[k] = log_weights[k] + dm;
    }
    target += log_sum_exp(lp);
  }
}

generated quantities {
  // Posterior community membership probabilities per sample
  matrix[N, K] community_probs;
  // MAP community assignment per sample
  array[N] int<lower=1, upper=K> z_hat;
  // Per-sample log-likelihood for LOO-CV
  vector[N] log_lik;

  for (i in 1:N) {
    vector[K] log_weights = log_softmax(log_mixing_weights[i]');
    vector[K] lp;
    for (k in 1:K) {
      vector[S] alpha_k = alpha * pi[k];
      real dm = log_mc[i] + lgamma(alpha) - lgamma(Ni[i] + alpha);
      for (s in 1:S)
        dm += lgamma(X[i, s] + alpha_k[s]) - lgamma(alpha_k[s]);
      lp[k] = log_weights[k] + dm;
    }
    log_lik[i]         = log_sum_exp(lp);
    vector[K] post     = softmax(lp);
    community_probs[i] = post';
    z_hat[i]           = categorical_rng(post);
  }
}
