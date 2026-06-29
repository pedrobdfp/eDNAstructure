# eDNAstructure

**Dirichlet-Multinomial Mixture Models for eDNA Metabarcoding Community Structure**

`eDNAstructure` is an R package for fitting Bayesian Dirichlet-Multinomial Mixture (DMM) models to environmental DNA (eDNA) read count data from metabarcoding surveys. Given a sample Ã— taxon count matrix and optional environmental covariates, the model identifies latent ecological communities, estimates their taxonomic compositions, and quantifies how environmental gradients drive community membership â€” all within a fully Bayesian framework with principled uncertainty quantification.

---

## Installation

### Step 1 â€” Install a C++ toolchain

Stan compiles models to C++ and requires a toolchain on your machine:

- **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- **macOS**: Run `xcode-select --install` in Terminal
- **Linux**: Install `build-essential` (Ubuntu/Debian) or equivalent

### Step 2 â€” Install rstan

```r
install.packages("rstan")
```

Verify it works before proceeding:

```r
library(rstan)
example(stan_model, package = "rstan", run.dontrun = TRUE)
```

If you see sampling output without errors, Stan is ready. Full guide: <https://mc-stan.org/rstan/articles/rstan.html>

### Step 3 â€” Install eDNAstructure

```r
install.packages("remotes")
remotes::install_github("pedrobdfp/eDNA_structure", upgrade = "never")
```

During installation, a large amount of black text will appear â€” this is the Stan model compiling to C++. It only happens once. Every subsequent call to `eDNA_dmm()` goes straight to sampling with no compilation output.

### Dependencies

Installed automatically:

| Package | Purpose |
|---------|---------|
| `rstan` (â‰¥ 2.21) | Bayesian inference via Stan |
| `ggplot2` (â‰¥ 3.4) | All visualizations |
| `dplyr`, `tidyr` | Data manipulation |
| `vegan` (â‰¥ 2.6) | NMDS ordination |
| `posterior` (â‰¥ 1.4) | MCMC diagnostics (ESS, Rhat) |
| `loo` (â‰¥ 2.6) | Leave-one-out cross-validation |
| `scales` | Axis formatting |

---

## Quick start

```r
library(eDNAstructure)
library(dplyr)    # for pipe and data manipulation
library(ggplot2)  # for plot customization

# Option A â€” use the built-in example dataset
data <- get_example_data()

# Option B â€” simulate your own dataset with known ground truth
data <- simulate_eDNA_survey(
  n_communities         = 4,
  n_species             = 40,
  samples_per_community = 5,
  seed                  = 2026
)

# Inspect raw species composition before fitting
plot_true_compositions(
  data$counts,
  metadata  = data$covariates,
  facet_var = "TrueCommunity"
)

# Fit the model with a given number of communities (K)
fit <- eDNA_dmm(
  counts     = data$counts,
  covariates = data$covariates[, c("Depth", "Distance_shore")],
  K          = 4
)

print(fit)
summary(fit)

# Or select K using LOO cross-validation
loo_result <- eDNA_loo(data$counts, data$covariates[, c("Depth", "Distance_shore")],
                       K_range = 2:5)
loo_result$plot


# The loo_result object stores all fitted models â€” no need to refit
# Extract the K=4 model directly:
fit <- loo_result$fits[["K4"]]

# Or using the K value as a number:
K_best <- 4
fit <- loo_result$fits[[paste0("K", K_best)]]

# Confirm what you have:
print(fit)

##You can also plot the results!

# Structure plot â€” one bar per sample, colored by community membership probability
eDNA_dmm_structure(fit, metadata = data$covariates,
                   facet_var = "TrueCommunity", sort_var = "Depth")

# NMDS ordination colored by community assignment
eDNA_dmm_nmds(fit)$plot

# Prior vs posterior distributions for covariate effects
eDNA_dmm_beta(fit)$plot
```

> **For a complete walkthrough** â€” including step-by-step simulation, data formatting, K selection, all visualization options, parameter recovery, and troubleshooting â€” see the **[full tutorial vignette](vignettes/tutorial.Rmd)**. It is designed to be read start to finish and assumes no prior familiarity with Bayesian mixture models.

---

## Input data format

### Count matrix

The primary input to `eDNA_dmm()` is a **sample Ã— taxon** matrix of non-negative integer read counts:

- **Rows** = samples (stations, replicates, individuals, etc.)
- **Columns** = taxa or ASVs â€” taxonomic annotation is not required
- **Values** = raw integer read counts (do not normalize)

```r
data$counts[1:3, 1:5]
#          Sp_1  Sp_2  Sp_3  Sp_4  Sp_5
# STN_001   412   310   121    73     0
# STN_002   389   275    98    61    14
# STN_003    52    41   487   312   208
```

If your data are in long format, convert them first:

```r
library(tidyr)
count_matrix <- long_df |>
  pivot_wider(names_from = taxon, values_from = reads, values_fill = 0) |>
  tibble::column_to_rownames("sample_id") |>
  as.matrix()
```

### Covariate data frame

A **sample Ã— covariate** data frame in the same row order as the count matrix. Covariates are Z-score standardized internally by default.

```r
head(data$covariates)
#   sample_id  TrueCommunity  Depth  Distance_shore
# 1   STN_001              1     82             198
# 2   STN_002              1     79             204
# 3   STN_003              2     11             197
```

---

## Functions

### `eDNA_dmm()` â€” Fit the DMM

The core function. Fits a Dirichlet-Multinomial Mixture model via Stan and returns an `edna_dmm_fit` object.

```r
fit <- eDNA_dmm(
  counts           = my_counts,   # sample Ã— taxon integer count matrix
  covariates       = my_covs,     # sample Ã— covariate data frame, or NULL
  K                = 4,           # number of latent communities to fit
  scale_covariates = TRUE,        # Z-score standardize covariates (strongly recommended)
  chains           = 1,           # number of MCMC chains (see note on label switching below)
  iter             = 4000,        # total iterations per chain (including warmup)
  warmup           = 2000,        # warmup iterations to discard
  adapt_delta      = 0.95,        # HMC target acceptance rate; increase to 0.99 if divergences
  max_treedepth    = 12,          # increase to 14â€“15 if "max treedepth exceeded" warnings
  seed             = 13,          # random seed for reproducibility
  conc             = 0.5,         # Dirichlet prior concentration: < 1 = sparse communities
  alpha_shape      = 5,           # Gamma prior shape for overdispersion parameter alpha
  alpha_rate       = 2,           # Gamma prior rate  (prior mean = shape/rate = 2.5)
  verbose          = TRUE         # print sampling progress
)
```

The returned `edna_dmm_fit` object contains:

| Element | Description |
|---------|-------------|
| `sample_info` | Data frame: posterior membership probabilities and MAP assignment per sample |
| `pi_mean` | Matrix [K Ã— S]: posterior mean community compositions |
| `beta_summary` | Data frame: covariate coefficient summaries with ESS and reliability |
| `alpha_mean` | Scalar: posterior mean overdispersion |
| `stan_fit` | Raw `rstan::stanfit` object for advanced diagnostics |

> **On single chains:** Mixture models suffer from label switching across chains â€” "Community 1" in chain A may map to "Community 2" in chain B, making multi-chain Rhat diagnostics meaningless. A single long chain sidesteps this. Use within-chain ESS (reported by `summary()`) as your convergence criterion.

---

### `eDNA_dmm_structure()` â€” Structure bar plot

Produces a STRUCTURE-style plot: one vertical bar per sample, divided into colored segments by posterior community membership probability.

```r
p <- eDNA_dmm_structure(
  fit,
  metadata         = my_metadata,    # data frame with additional sample variables
  sample_id_col    = "sample_id",    # column in metadata matching sample IDs in fit
  facet_var        = "TrueCommunity", # facet panels by this variable (e.g. site, year, depth)
  sort_var         = "Depth",        # sort samples within each panel by this variable
  community_colors = NULL,           # named hex vector (e.g. c("Community 1" = "#E63946"))
                                     # or NULL for automatic HCL palette
  bar_width        = 0.9,            # bar width (0â€“1); 1 = no gaps between bars
  x_text           = FALSE,          # show sample ID labels on x-axis?
  base_size        = 11,             # base font size in points
  title            = NULL,           # plot title; NULL = auto-generated
  subtitle         = NULL,           # plot subtitle; NULL = auto-generated
  legend_position  = "bottom"        # "bottom", "right", "left", "top", or "none"
)
```

Returns a `ggplot2` object â€” save with `ggsave()` or extend with additional ggplot2 layers.

---

### `eDNA_dmm_nmds()` â€” NMDS ordination

Runs NMDS on community dissimilarities and plots samples colored by their MAP community assignment. Point **size** reflects assignment certainty: larger points are more confidently assigned to a single community.

```r
result <- eDNA_dmm_nmds(
  fit,
  k                = 2,          # NMDS dimensions (2 or 3); increase if stress > 0.2
  nmds_axes        = c(1, 2),    # which two axes to display; e.g. c(1,3) for axes 1 and 3
  distance         = "bray",     # dissimilarity metric passed to vegan::vegdist()
  use_edna_index   = TRUE,       # apply eDNA index transform before computing distances
  trymax           = 100,        # maximum random NMDS starts (more = less risk of local optima)
  seed             = 42,         # random seed for NMDS
  show_ellipse     = TRUE,       # draw 95% confidence ellipse per community?
  ellipse_type     = "t",        # ellipse type: "t" (robust) or "norm" (normal-based)
  community_colors = NULL,       # named hex vector or NULL for automatic palette
  size_range       = c(1.5, 5),  # point size range: c(min, max) mapped to
                                 # 50% certainty (smallest) â†’ 100% certainty (largest)
  alpha            = 0.85,       # point transparency (0 = invisible, 1 = opaque)
  base_size        = 13,
  title            = NULL,
  subtitle         = NULL,
  legend_position  = "right"
)

result$plot   # ggplot2 object
result$nmds   # vegan::metaMDS object (access stress value, species scores, etc.)
```

---

### `eDNA_dmm_beta()` â€” Covariate effects

Overlays the prior and posterior distributions for each softmax regression coefficient. A posterior pulled away from the prior is evidence that the covariate genuinely predicts community membership.

```r
result <- eDNA_dmm_beta(
  fit,
  layout             = "joint",    # "joint": communities overlaid per covariate panel
                                   # "separate": one row per community, one column per covariate
  covariates_to_plot = NULL,       # character vector of covariate names to include, or NULL for all
  show_intercept     = FALSE,      # include the intercept term?
  beta_prior_sd      = 1.0,        # prior SD â€” must match the Stan model (default: Normal(0,1))
  n_prior_samples    = 4000,       # prior draws for the density curve (more = smoother)
  community_colors   = NULL,       # named hex vector or NULL
  prior_color        = "grey60",   # fill color for the prior density
  prior_alpha        = 0.35,       # prior density transparency
  posterior_alpha    = 0.55,       # posterior density transparency
  show_annotations   = NULL,       # NULL = auto (shown for K=2 only); TRUE or FALSE to override
  base_size          = 13,
  title              = NULL,
  subtitle           = NULL
)

result$plot    # ggplot2 object
result$table   # data frame: mean, 90% CI, P(direction), ESS, reliability per coefficient
```

---

### `eDNA_loo()` â€” K selection via LOO cross-validation

Fits models across a range of K values and compares them using Leave-One-Out cross-validation. Returns an elbow plot and a comparison table to guide K selection.

```r
loo_result <- eDNA_loo(
  counts           = my_counts,
  covariates       = my_covs,
  K_range          = 2:5,        # integer vector of K values to evaluate
  scale_covariates = TRUE,
  chains           = 1,
  iter             = 4000,
  warmup           = 2000,
  adapt_delta      = 0.95,
  seed             = 13,
  conc             = 0.5,
  alpha_shape      = 5,
  alpha_rate       = 2,
  verbose          = TRUE
)

loo_result$plot        # LOO-ELPD elbow plot (higher = better; look for the elbow)
loo_result$loo_table   # data frame: K, LOO-ELPD, SE
loo_result$loo_compare # loo::loo_compare() output
loo_result$fits        # named list of edna_dmm_fit objects, one per K
```

---

### `plot_true_compositions()` â€” Raw species composition

Visualizes the observed species frequencies per sample as stacked bars â€” the same layout as `eDNA_dmm_structure()`, allowing direct before/after comparison. Most useful before fitting to inspect the raw community signal, and with simulated data where true community labels are known.

```r
p <- plot_true_compositions(
  counts,
  metadata        = my_metadata,    # data frame for faceting and sorting
  sample_id_col   = "sample_id",
  facet_var       = "TrueCommunity", # facet by known or hypothesized grouping
  sort_var        = "Depth",
  top_n           = 20,             # show top N taxa individually; rest collapsed to "Other"
  bar_width       = 0.95,
  base_size       = 11,
  title           = NULL,
  subtitle        = NULL,
  legend_position = "none"          # default none â€” too many taxa for a useful legend
)
```

Returns a `ggplot2` object.

---
### eDNA_dmm_compositions() â€” Posterior community compositions
Visualizes the posterior mean taxonomic composition of each latent community as stacked bars â€” the model's estimate of what each community "looks like" in species space. Colors match those used in eDNA_dmm_structure() and plot_true_compositions() for direct comparison.
```r
rp <- eDNA_dmm_compositions(
  fit,
  top_n           = 20,      # show top N taxa; rest collapsed to "Other"
  base_size       = 13,
  title           = NULL,
  subtitle        = NULL,
  legend_position = "right", # "right", "bottom", "left", "top", or "none"
  bar_width       = 0.7      # bar width (0â€“1)
)
```

Returns a ggplot2 object. The x-axis labels show community numbers (1, 2, 3â€¦). Pair with eDNA_dmm_structure() to connect community identities to sample assignments.

---

### `get_example_data()` â€” Built-in example dataset

Returns the built-in simulated dataset: 20 samples Ã— 40 taxa across 4 communities separated by depth and distance from shore. Generated by `simulate_eDNA_survey()` with known ground truth, so fitted parameters can be compared to the true values.

```r
data <- get_example_data()
# data$counts                â€” integer matrix [20 Ã— 40]
# data$covariates            â€” data frame: sample_id, TrueCommunity, Depth, Distance_shore
# data$community_compositions â€” true composition matrix [4 Ã— 40]
# data$metab_df              â€” raw simulated metabarcoding reads
# data$sample_metadata       â€” full simulation metadata
```

---

### Simulation pipeline

`eDNAstructure` includes a mechanistic simulation pipeline for generating eDNA datasets with known community structure. Use it for method validation, power analysis, or teaching. The full tutorial (`vignettes/tutorial.Rmd`) walks through the pipeline in detail.

```r
# Full pipeline in one call
sim <- simulate_eDNA_survey(
  n_communities             = 4,           # number of communities K
  n_species                 = 40,          # number of species S
  samples_per_community     = 5,           # sampling stations per community
  community_covariate_means = NULL,        # K Ã— P matrix of covariate means per community
                                           # NULL = default 2-covariate depth Ã— shore design
  covariate_sds             = NULL,        # K Ã— P SDs (NULL = 5 for all)
  mean_read_depth           = 10000,       # mean reads per sample
  bio_reps                  = 1,           # biological replicates per station
  seq_reps                  = 1,           # sequencing technical replicates per bio rep
  spillover                 = 0.15,        # fraction of species frequency leaking between communities
  shedding_error            = 0.3,         # lognormal SD on per-organism eDNA shedding
  decay_rate                = 0.1,         # exponential distance decay of eDNA signal
  seed                      = 42
)
# sim$counts      â€” ready for eDNA_dmm()
# sim$covariates  â€” ready for eDNA_dmm()
```

Or run each step individually for full control:

```r
community_mat   <- generate_community_compositions(
  n_communities     = 4,
  n_species         = 40,
  n_dominant_range  = c(2, 3),     # 2â€“3 high-frequency dominant species per community
  n_unique_low_freq = 4,           # species present only in one community
  community_groups  = c(1,2,1,2),  # group structure for shared species
  n_group_shared    = 4,           # species shared within each group
  spillover         = 0.15,        # cross-community frequency leakage
  seed              = 1
)

contrib_obj <- generate_contributors(
  community_compositions = community_mat,
  samples_per_community  = 5,
  n_contributors_range   = c(20, 100),   # organisms per sample
  n_distribution         = "Negative Binomial",  # overdispersed organism counts
  size_range             = c(1, 10),     # body size range (affects eDNA shedding)
  distance_range         = c(0, 100),    # distance from sampler (affects eDNA decay)
  seed                   = 1
)

eDNA_obj <- generate_eDNA(
  contributors_list = contrib_obj$contributors_list,
  shedding_rate     = 1000,   # baseline molecules shed per unit size
  beta              = 0.75,   # allometric exponent (metabolic scaling)
  decay_rate        = 0.1,    # exponential distance decay
  shedding_error    = 0.3,    # lognormal noise on shedding (SD on log scale)
  bottle_volume     = 0.1,    # fraction of local eDNA pool captured per bottle
  bio_reps          = 1
)

metab_df <- simulate_metabarcoding(
  eDNA_obj,
  mean_read_depth = 10000,   # mean total reads per sample
  read_depth_sd   = 0.2,     # lognormal SD for read depth variation
  error_sd        = 0.05,    # Gaussian noise added to species frequencies
  rep             = 1        # sequencing technical replicates per bottle
)

sample_metadata <- generate_sample_covariates(
  contributors_list    = contrib_obj$contributors_list,
  community_covariates = matrix(           # community-specific covariate means
    c(80, 200, 10, 200, 80, 20, 10, 20),
    nrow = 4, byrow = TRUE,
    dimnames = list(NULL, c("Depth", "Distance_shore"))
  )
)
```

| Function | Purpose |
|----------|---------|
| `simulate_eDNA_survey()` | Full pipeline in one call |
| `generate_community_compositions()` | K community frequency vectors over S species |
| `generate_contributors()` | Organisms shedding eDNA per sample |
| `generate_eDNA()` | Shedding, exponential decay, bottle sub-sampling |
| `simulate_metabarcoding()` | Amplification bias and multinomial read counts |
| `generate_sample_covariates()` | Environmental metadata drawn from community-specific distributions |

---

## The model

For sample *i*, the DMM marginalizes over a latent community assignment *z*_i:

1. **Compositions**: Ï€_k ~ Dirichlet(conc Â· **1**_S) for k = 1â€¦K
2. **Membership**: P(*z*_i = k) = softmax(Î²_0k + Î²_1k Â· x_1i + â€¦ + Î²_Pk Â· x_Pi), community K = reference
3. **Counts**: **x**_i | *z*_i = k ~ DirichletMultinomial(N_i, Î± Â· Ï€_k)

The global overdispersion Î± absorbs both technical (PCR, sequencing) and ecological compositional variance. Marginalizing over *z*_i makes inference exact.

---

## Frequently asked questions

**Why only one chain?**
Label switching: "Community 1" in chain A may map to "Community 2" in chain B. Multi-chain Rhat values are pathological even when each chain converges perfectly. One long chain avoids this. Check within-chain ESS instead (printed by `summary()`).

**I have divergent transitions. What do I do?**
Increase `adapt_delta` toward `0.99`. If they persist, try lower K or verify your count matrix has no all-zero rows or columns.

**Can I use raw ASVs instead of taxonomy-collapsed counts?**
Yes. The model treats each column as a compositional unit and does not use taxonomy. ASVs give finer resolution; taxa collapse dimensionality and often converge faster.

**How do I include year as a covariate?**
Pass it as a numeric column. But if you have only a few discrete years, the linearity assumption may be too strong â€” consider fitting without year and testing it post-hoc via multinomial regression on the posterior assignments.

**The first run takes forever â€” is something wrong?**
No. Stan compiles the model to C++ on the first call after installation (1â€“2 minutes). All subsequent calls skip compilation. This is normal behavior for any rstan-based package.

---

## Citation

If you use this package in published research, please cite:

> BrandÃ£o-Dias et al. (year). Multinomial mixture models from environmental DNA reveal
> depth stability and dynamic surface turnover of marine vertebrate communities. *Under review.*

> BrandÃ£o-Dias et al. (year). eDNAstructure: Dirichlet-Multinomial Mixture Models for
> eDNA Metabarcoding Community Structure. R package version 0.1.0.
> https://github.com/pedrobdfp/eDNA_structure

Please also cite Stan:

> Carpenter B. et al. (2017). Stan: A probabilistic programming language.
> *Journal of Statistical Software*, 76(1).

---

## License

MIT © Pedro Brandão-Dias. See [LICENSE](LICENSE) for details.
