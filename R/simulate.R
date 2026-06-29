# =============================================================================
# Simulation functions for eDNA metabarcoding community structure
# =============================================================================
# These functions implement a mechanistic simulation pipeline:
#   1. generate_community_compositions() — K community frequency vectors
#   2. generate_contributors()           — organisms shedding eDNA per sample
#   3. generate_eDNA()                   — eDNA shedding, decay, sub-sampling
#   4. simulate_metabarcoding()          — amplification bias + read counts
#   5. generate_sample_covariates()      — environmental metadata
#   6. simulate_eDNA_survey()            — full pipeline in one call
# =============================================================================


#' Generate community composition matrices for eDNA simulation
#'
#' @description
#' Builds K community compositions over S species. Each community has a set of
#' dominant species (high frequency, mostly exclusive to that community),
#' unique low-frequency species (present only in that community), optional
#' group-shared species (shared within a set of related communities), and a
#' shared background pool. Spillover parameters control how much a species'
#' frequency leaks into non-native communities.
#'
#' @param n_communities A positive integer: number of communities (K). Default `4`.
#' @param n_species A positive integer: total number of species (S). Default `50`.
#' @param n_shared An integer or NULL: number of species in the shared background
#'   pool. NULL (default) uses all remaining species after dominants, uniques,
#'   and group-shared are assigned.
#' @param spillover A number in [0, 1]: mean fraction of a species' total
#'   frequency that leaks into non-native communities. Default `0.15`.
#'   `0` = perfectly exclusive communities; `1` = no community structure.
#' @param spillover_concentration A positive number: concentration of the Beta
#'   prior governing per-species spillover variance. Higher values force all
#'   species to have spillover close to the mean. Default `10`.
#' @param community_groups An integer vector of length K assigning communities
#'   to groups (e.g., `c(1, 2, 1, 2)` means communities 1 and 3 share a group).
#'   Group-shared species appear in all communities of the same group but not
#'   others. Pass `NULL` (default) for no group structure.
#' @param n_group_shared A non-negative integer: number of species shared within
#'   each group (not across groups). Default `0`.
#' @param group_shared_freq A length-2 numeric vector: min and max frequency for
#'   group-shared species. Default `c(0.02, 0.08)`.
#' @param group_shared_presence A number in [0, 1]: probability that a
#'   group-shared species is present in a given community of its group.
#'   Default `0.9`.
#' @param n_dominant_range A length-2 integer vector: min and max number of
#'   dominant species per community. Default `c(2, 3)`.
#' @param n_unique_low_freq A non-negative integer: number of unique low-frequency
#'   species per community. Default `3`.
#' @param dominant_freq A length-2 numeric vector: min and max frequency for
#'   dominant species. Default `c(0.15, 0.30)`.
#' @param low_unique_freq A length-2 numeric vector: min and max frequency for
#'   unique low-frequency species. Default `c(0.005, 0.04)`.
#' @param shared_freq A length-2 numeric vector: min and max frequency for
#'   shared background species. Default `c(0.001, 0.05)`.
#' @param shared_presence A number in [0, 1]: probability that a shared
#'   background species appears in a given community. Default `0.4`.
#' @param seed An integer random seed for reproducibility, or `NULL`. Default `NULL`.
#'
#' @return A numeric matrix [K × S] of community compositions (rows sum to 1).
#'   Row names are `"Community_1"` through `"Community_K"`. Column names are
#'   `"Sp_1"` through `"Sp_S"`. The matrix carries attributes documenting which
#'   species are dominants, uniques, group-shared, shared, and absent per community.
#'
#' @seealso [generate_contributors()], [simulate_eDNA_survey()]
#'
#' @examples
#' cm <- generate_community_compositions(
#'   n_communities = 3,
#'   n_species     = 30,
#'   seed          = 42
#' )
#' dim(cm)       # 3 x 30
#' rowSums(cm)   # all 1
#'
#' @export
generate_community_compositions <- function(
    n_communities           = 4,
    n_species               = 50,
    n_shared                = NULL,
    spillover               = 0.15,
    spillover_concentration = 10,
    community_groups        = NULL,
    n_group_shared          = 0,
    group_shared_freq       = c(0.02, 0.08),
    group_shared_presence   = 0.9,
    n_dominant_range        = c(2, 3),
    n_unique_low_freq       = 3,
    dominant_freq           = c(0.15, 0.30),
    low_unique_freq         = c(0.005, 0.04),
    shared_freq             = c(0.001, 0.05),
    shared_presence         = 0.4,
    seed                    = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  if (spillover < 0 || spillover > 1)
    rlang::abort("`spillover` must be between 0 and 1.")
  if (spillover_concentration <= 0)
    rlang::abort("`spillover_concentration` must be > 0.")

  K <- n_communities
  S <- n_species

  if (is.null(community_groups)) {
    community_groups <- seq_len(K)
    n_group_shared   <- 0
  } else {
    if (length(community_groups) != K)
      rlang::abort("`community_groups` must have length equal to `n_communities`.")
  }
  group_levels <- unique(community_groups)
  n_groups     <- length(group_levels)

  community_mat <- matrix(0, K, S,
    dimnames = list(paste0("Community_", seq_len(K)),
                    paste0("Sp_", seq_len(S))))
  idx <- 1L

  # Dominant species
  community_dominants <- vector("list", K)
  for (k in seq_len(K)) {
    n_dom <- if (n_dominant_range[1] == n_dominant_range[2]) {
      n_dominant_range[1]
    } else {
      sample(n_dominant_range[1]:n_dominant_range[2], 1)
    }
    community_dominants[[k]] <- idx:(idx + n_dom - 1L)
    idx <- idx + n_dom
  }

  # Unique low-frequency species
  community_uniques <- vector("list", K)
  for (k in seq_len(K)) {
    community_uniques[[k]] <- idx:(idx + n_unique_low_freq - 1L)
    idx <- idx + n_unique_low_freq
  }

  # Group-shared species
  group_shared <- vector("list", n_groups)
  names(group_shared) <- as.character(group_levels)
  if (n_group_shared > 0) {
    for (g_idx in seq_along(group_levels)) {
      group_shared[[g_idx]] <- idx:(idx + n_group_shared - 1L)
      idx <- idx + n_group_shared
    }
  }

  if (idx > S + 1L)
    rlang::abort(paste0(
      "Not enough species (n_species = ", S, "). ",
      "Reduce n_dominant_range, n_unique_low_freq, or n_group_shared, ",
      "or increase n_species."
    ))

  # Shared background pool
  n_remaining <- S - (idx - 1L)
  if (is.null(n_shared)) {
    n_shared_use <- n_remaining
  } else {
    if (n_shared < 0)
      rlang::abort("`n_shared` must be >= 0.")
    if (n_shared > n_remaining)
      rlang::abort(paste0("`n_shared` (", n_shared, ") exceeds remaining species (",
                          n_remaining, ")."))
    n_shared_use <- n_shared
  }
  shared_species <- if (n_shared_use > 0) idx:(idx + n_shared_use - 1L) else integer(0)
  absent_species <- if (n_shared_use < n_remaining) (idx + n_shared_use):S else integer(0)

  # Fill frequencies
  for (k in seq_len(K)) {
    for (s in community_dominants[[k]])
      community_mat[k, s] <- stats::runif(1, dominant_freq[1], dominant_freq[2])
    for (s in community_uniques[[k]])
      community_mat[k, s] <- stats::runif(1, low_unique_freq[1], low_unique_freq[2])
    if (n_group_shared > 0) {
      g_idx <- which(group_levels == community_groups[k])
      for (s in group_shared[[g_idx]]) {
        if (stats::runif(1) < group_shared_presence)
          community_mat[k, s] <- stats::runif(1, group_shared_freq[1], group_shared_freq[2])
      }
    }
    for (s in shared_species) {
      if (stats::runif(1) < shared_presence)
        community_mat[k, s] <- stats::runif(1, shared_freq[1], shared_freq[2])
    }
  }

  # Per-species exclusivity and spillover
  unique_set       <- unlist(community_uniques)
  group_shared_set <- unlist(group_shared)
  exclusivity      <- numeric(S)
  exclusivity[unique_set]       <- 1
  exclusivity[group_shared_set] <- 1
  spillable <- setdiff(seq_len(S), c(unique_set, group_shared_set))
  mean_excl <- 1 - spillover
  a <- mean_excl       * spillover_concentration
  b <- (1 - mean_excl) * spillover_concentration
  exclusivity[spillable] <- stats::rbeta(length(spillable), a, b)

  for (s in seq_len(S)) {
    if (s %in% unique_set || s %in% group_shared_set || s %in% absent_species) next
    col_s <- community_mat[, s]
    total <- sum(col_s)
    if (total == 0) next
    e_s          <- exclusivity[s]
    native_part  <- col_s * e_s
    spill_total  <- total * (1 - e_s)
    rand_weights <- stats::rgamma(K, shape = 2)
    rand_weights <- rand_weights / sum(rand_weights)
    community_mat[, s] <- native_part + spill_total * rand_weights
  }

  community_mat <- community_mat / rowSums(community_mat)

  attr(community_mat, "dominants")        <- community_dominants
  attr(community_mat, "uniques")          <- community_uniques
  attr(community_mat, "group_shared")     <- group_shared
  attr(community_mat, "shared")           <- shared_species
  attr(community_mat, "absent")           <- absent_species
  attr(community_mat, "exclusivity")      <- exclusivity
  attr(community_mat, "spillover")        <- spillover
  attr(community_mat, "community_groups") <- community_groups
  community_mat
}


#' Generate eDNA contributors for each simulated sample
#'
#' @description
#' For each community × sample replicate combination, draws a set of individual
#' organisms (contributors) whose species are sampled from the community's
#' frequency vector. Each contributor is assigned a body size and a distance
#' from the sampler. Species flagged as grouped (schooling) share distances
#' within groups; ungrouped (solitary) species get independent distances.
#'
#' @param community_compositions A [K × S] community composition matrix from
#'   [generate_community_compositions()].
#' @param samples_per_community A positive integer: number of independent samples
#'   (stations) drawn from each community. Default `3`.
#' @param n_contributors_range A length-2 integer vector: min and max number of
#'   individual organisms contributing eDNA per sample. Default `c(20, 200)`.
#' @param grouped_species_vec A logical vector of length S: `TRUE` if the species
#'   is gregarious (school members share a distance draw); `FALSE` for solitary
#'   (independent distances). Pass `NULL` (default) for all species solitary.
#' @param size_range A length-2 numeric vector: min and max organism body size
#'   (arbitrary units, affects eDNA shedding). Default `c(1, 10)`.
#' @param distance_range A length-2 numeric vector: min and max distance from
#'   organism to sampler (same units as `decay_rate` in [generate_eDNA()]).
#'   Default `c(0, 100)`.
#' @param group_size_range A length-2 integer vector: min and max school size
#'   for grouped species. Default `c(2, 10)`.
#' @param n_distribution A string: `"Random"` (uniform between
#'   `n_contributors_range`) or `"Negative Binomial"` (overdispersed, more
#'   realistic). Default `"Negative Binomial"`.
#' @param seed An integer random seed, or `NULL`. Default `NULL`.
#'
#' @return A list with elements:
#' \describe{
#'   \item{`contributors_list`}{A list of data frames, one per sample, each with
#'     columns `SampleID`, `TrueCommunity`, `ContributorID`, `Species`,
#'     `Distance`, `Size`.}
#'   \item{`community_compositions`}{The input composition matrix (passed through).}
#'   \item{`grouped_species_vec`}{The logical grouping vector (passed through).}
#' }
#'
#' @seealso [generate_community_compositions()], [generate_eDNA()]
#'
#' @examples
#' cm   <- generate_community_compositions(n_communities = 2, n_species = 20, seed = 1)
#' cont <- generate_contributors(cm, samples_per_community = 3, seed = 1)
#' length(cont$contributors_list)   # 6 samples (2 communities x 3)
#'
#' @export
generate_contributors <- function(
    community_compositions,
    samples_per_community = 3,
    n_contributors_range  = c(20, 200),
    grouped_species_vec   = NULL,
    size_range            = c(1, 10),
    distance_range        = c(0, 100),
    group_size_range      = c(2, 10),
    n_distribution        = c("Negative Binomial", "Random"),
    seed                  = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  n_distribution <- match.arg(n_distribution)

  K <- nrow(community_compositions)
  S <- ncol(community_compositions)

  if (is.null(grouped_species_vec)) {
    grouped_species_vec <- rep(FALSE, S)
  } else if (length(grouped_species_vec) != S) {
    rlang::abort("`grouped_species_vec` must have length equal to the number of species (columns).")
  }

  sample_list <- list()
  counter     <- 0L
  mu_val      <- n_contributors_range[2] / 5

  for (k in seq_len(K)) {
    comm_freq <- community_compositions[k, ]
    comm_freq <- comm_freq / sum(comm_freq)

    for (s_idx in seq_len(samples_per_community)) {
      counter <- counter + 1L

      n_contributors <- if (n_distribution == "Random") {
        sample(n_contributors_range[1]:n_contributors_range[2], 1)
      } else {
        repeat {
          cand <- stats::rnbinom(1, size = 1, mu = mu_val)
          if (cand >= n_contributors_range[1] && cand <= n_contributors_range[2]) break
        }
        cand
      }

      species_vec <- sample(seq_len(S), n_contributors, replace = TRUE, prob = comm_freq)

      distances <- numeric(n_contributors)
      for (sp in unique(species_vec)) {
        idxs <- which(species_vec == sp)
        n_sp <- length(idxs)
        if (grouped_species_vec[sp]) {
          remaining <- n_sp; taken <- 0L
          while (remaining > 0) {
            max_gs <- min(group_size_range[2], remaining)
            min_gs <- min(group_size_range[1], remaining)
            grp <- if (max_gs == min_gs) min_gs else sample(min_gs:max_gs, 1)
            d   <- stats::runif(1, distance_range[1], distance_range[2])
            distances[idxs[(taken + 1L):(taken + grp)]] <- d
            taken     <- taken + grp
            remaining <- remaining - grp
          }
        } else {
          distances[idxs] <- stats::runif(n_sp, distance_range[1], distance_range[2])
        }
      }

      sizes <- stats::runif(n_contributors, size_range[1], size_range[2])

      sample_list[[counter]] <- data.frame(
        SampleID      = counter,
        TrueCommunity = k,
        ContributorID = seq_len(n_contributors),
        Species       = species_vec,
        Distance      = distances,
        Size          = sizes
      )
    }
  }

  list(
    contributors_list      = sample_list,
    community_compositions = community_compositions,
    grouped_species_vec    = grouped_species_vec
  )
}


#' Simulate eDNA shedding, environmental decay, and bottle sub-sampling
#'
#' @description
#' Each contributor sheds eDNA proportional to its body size raised to an
#' allometric exponent, with optional lognormal multiplicative noise. eDNA
#' decays exponentially with distance to the sampler. The resulting molecular
#' pool is sub-sampled (with replacement) to represent one or more biological
#' replicates (bottles).
#'
#' @param contributors_list A list of contributor data frames from
#'   [generate_contributors()] (the `$contributors_list` element).
#' @param shedding_error A non-negative number: standard deviation of lognormal
#'   multiplicative noise on per-contributor shedding. `0` = no noise (default).
#' @param decay_rate A positive number: exponential decay rate with distance.
#'   Default `0.1`. Higher values = eDNA signal drops off faster with distance.
#' @param shedding_rate A positive number: baseline shedding multiplier.
#'   Default `1`. Scale up (e.g., `1000`) to get more molecules and thus
#'   more stable sub-sampling.
#' @param beta A positive number: allometric exponent relating body size to
#'   eDNA shedding. Default `0.75` (standard metabolic scaling). `1.0` = linear.
#' @param bio_reps A positive integer: number of independent bottle sub-samples
#'   drawn from each local eDNA pool. Default `1`.
#' @param bottle_volume A number in (0, 1]: fraction of the local eDNA pool
#'   captured in each bottle. Default `0.1`.
#'
#' @return A list with:
#' \describe{
#'   \item{`Summary`}{Data frame with one row per (sample, bio_rep) combination:
#'     `SampleID`, `TrueCommunity`, `BioRep`, `Total_eDNA`, `Total_contributors`,
#'     `Avg_distance`, `Closest_contributor`, `SD_distance`.}
#'   \item{`Species_frequencies`}{Data frame with relative species frequencies
#'     per (sample, bio_rep): `SampleID`, `TrueCommunity`, `BioRep`, `Species`,
#'     `Relative_frequency`.}
#' }
#'
#' @seealso [generate_contributors()], [simulate_metabarcoding()]
#'
#' @export
generate_eDNA <- function(
    contributors_list,
    shedding_error = 0,
    decay_rate     = 0.1,
    shedding_rate  = 1,
    beta           = 0.75,
    bio_reps       = 1,
    bottle_volume  = 0.1
) {
  summary_list <- list()
  species_list <- list()

  for (i in seq_along(contributors_list)) {
    sample_df <- contributors_list[[i]]
    sample_id <- unique(sample_df$SampleID)
    true_comm <- unique(sample_df$TrueCommunity)
    n_contrib <- nrow(sample_df)

    base_shed  <- shedding_rate * sample_df$Size^beta
    mult_err   <- if (shedding_error == 0) rep(1, n_contrib) else
      stats::rlnorm(n_contrib, meanlog = 0, sdlog = shedding_error)
    edna_counts <- round(base_shed * mult_err * exp(-decay_rate * sample_df$Distance))

    local_pool <- data.frame(ContributorID = integer(0), Species = integer(0))
    for (j in seq_len(n_contrib)) {
      cnt <- edna_counts[j]
      if (cnt > 0) {
        local_pool <- rbind(local_pool,
          data.frame(ContributorID = rep(sample_df$ContributorID[j], cnt),
                     Species       = rep(sample_df$Species[j],       cnt)))
      }
    }
    total_eDNA <- nrow(local_pool)

    for (r in seq_len(bio_reps)) {
      if (total_eDNA == 0) {
        sub_pool <- data.frame(ContributorID = integer(0), Species = integer(0))
      } else {
        sub_n <- max(0L, floor(total_eDNA * bottle_volume))
        if (sub_n == 0L) {
          sub_pool <- data.frame(ContributorID = integer(0), Species = integer(0))
        } else {
          sel      <- sample.int(total_eDNA, sub_n, replace = TRUE)
          sub_pool <- local_pool[sel, ]
        }
      }
      replicate_eDNA <- nrow(sub_pool)
      uniq_contrib   <- unique(sub_pool$ContributorID)

      if (replicate_eDNA == 0) {
        sp_freq <- numeric(0); sp_names <- character(0)
        avg_d <- NA_real_; min_d <- NA_real_; sd_d <- NA_real_
      } else {
        tbl      <- table(sub_pool$Species)
        sp_freq  <- tbl / sum(tbl)
        sp_names <- names(tbl)
        d_sub    <- sample_df$Distance[sample_df$ContributorID %in% uniq_contrib]
        avg_d    <- mean(d_sub); min_d <- min(d_sub); sd_d <- stats::sd(d_sub)
      }

      summary_list[[length(summary_list) + 1L]] <- data.frame(
        SampleID            = sample_id,
        TrueCommunity       = true_comm,
        BioRep              = r,
        Total_eDNA          = replicate_eDNA,
        Total_contributors  = length(uniq_contrib),
        Avg_distance        = avg_d,
        Closest_contributor = min_d,
        SD_distance         = sd_d
      )

      species_list[[length(species_list) + 1L]] <- if (length(sp_freq) > 0) {
        data.frame(SampleID           = rep(sample_id, length(sp_freq)),
                   TrueCommunity      = rep(true_comm, length(sp_freq)),
                   BioRep             = rep(r,         length(sp_freq)),
                   Species            = sp_names,
                   Relative_frequency = as.numeric(sp_freq))
      } else {
        data.frame(SampleID = integer(0), TrueCommunity = integer(0),
                   BioRep = integer(0), Species = character(0),
                   Relative_frequency = numeric(0))
      }
    }
  }

  list(
    Summary            = do.call(rbind, summary_list),
    Species_frequencies = do.call(rbind, species_list)
  )
}


#' Simulate metabarcoding sequencing from eDNA frequencies
#'
#' @description
#' Applies optional per-species amplification bias to the bottle's species
#' frequencies, adds symmetric Gaussian noise, renormalizes, then draws read
#' counts from a multinomial distribution with lognormally-variable read depth.
#'
#' @param eDNA_data The output of [generate_eDNA()].
#' @param mean_read_depth A positive number: mean total reads per sample.
#'   Default `10000`.
#' @param read_depth_sd A non-negative number: standard deviation of the
#'   lognormal distribution for read depth. Default `0.2`.
#' @param error_sd A non-negative number: standard deviation of Gaussian noise
#'   added to species frequencies before multinomial sampling. Default `0.1`.
#' @param rep A positive integer: number of sequencing technical replicates per
#'   biological replicate. Default `1`.
#' @param amplification_bias A numeric vector of length S (number of species):
#'   per-species amplification multipliers. Values > 1 inflate that species'
#'   apparent frequency; values < 1 deflate it. Pass `NULL` (default) for
#'   unbiased amplification.
#'
#' @return A data frame with columns: `SampleID`, `BioRep`, `SeqRep`, `Species`
#'   (integer index), `Counts`.
#'
#' @seealso [generate_eDNA()], [simulate_eDNA_survey()]
#'
#' @export
simulate_metabarcoding <- function(
    eDNA_data,
    mean_read_depth    = 10000,
    read_depth_sd      = 0.2,
    error_sd           = 0.1,
    rep                = 1,
    amplification_bias = NULL
) {
  summary_df <- eDNA_data$Summary
  sp_df      <- eDNA_data$Species_frequencies

  if (!"BioRep" %in% colnames(summary_df)) summary_df$BioRep <- 1L
  if (!"BioRep" %in% colnames(sp_df))      sp_df$BioRep      <- 1L

  pairs   <- unique(summary_df[, c("SampleID", "BioRep")])
  results <- list()

  for (i in seq_len(nrow(pairs))) {
    sid <- pairs$SampleID[i]
    bio <- pairs$BioRep[i]
    sub <- sp_df[sp_df$SampleID == sid & sp_df$BioRep == bio, ]
    if (nrow(sub) == 0) next

    theta_o <- sub$Relative_frequency
    sp_ids  <- as.integer(as.character(sub$Species))

    if (!is.null(amplification_bias)) {
      bias_vec <- amplification_bias[sp_ids]
      theta_o  <- theta_o * bias_vec
      if (sum(theta_o) > 0) theta_o <- theta_o / sum(theta_o)
    }

    for (seq_r in seq_len(rep)) {
      R <- round(stats::rlnorm(1, meanlog = log(mean_read_depth), sdlog = read_depth_sd))
      theta <- numeric(length(theta_o))
      for (j in seq_along(theta_o)) {
        repeat {
          val <- theta_o[j] + stats::rnorm(1, 0, error_sd)
          if (val >= 0) { theta[j] <- val; break }
        }
      }
      theta <- if (sum(theta) > 0) theta / sum(theta) else rep(1 / length(theta), length(theta))
      counts <- if (R <= 0) rep(0L, length(theta)) else as.vector(stats::rmultinom(1, R, theta))

      results[[length(results) + 1L]] <- data.frame(
        SampleID = sid, BioRep = bio, SeqRep = seq_r,
        Species  = sp_ids, Counts = counts
      )
    }
  }
  do.call(rbind, results)
}


#' Generate environmental covariate metadata for simulated samples
#'
#' @description
#' For each simulated sample, draws environmental covariate values from a
#' community-specific multivariate normal distribution. This encodes the
#' ecological assumption that community membership is predictable from
#' environmental conditions.
#'
#' @param contributors_list A list of contributor data frames from
#'   [generate_contributors()] (the `$contributors_list` element).
#' @param community_covariates A [K × P] numeric matrix of community-specific
#'   covariate means. Row k gives the mean covariate values for community k.
#'   Column names become covariate names.
#' @param covariate_sds A [K × P] numeric matrix of standard deviations, or
#'   `NULL` (default) to use SD = 5 for all covariates and communities.
#' @param seed An integer random seed, or `NULL`. Default `NULL`.
#'
#' @return A data frame with one row per sample: `SampleID`, `TrueCommunity`,
#'   and one column per covariate (named after the columns of
#'   `community_covariates`).
#'
#' @seealso [generate_contributors()], [simulate_eDNA_survey()]
#'
#' @export
generate_sample_covariates <- function(
    contributors_list,
    community_covariates,
    covariate_sds = NULL,
    seed          = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  P <- ncol(community_covariates)
  if (is.null(covariate_sds)) {
    covariate_sds <- matrix(5, nrow = nrow(community_covariates), ncol = P,
                            dimnames = dimnames(community_covariates))
  }

  rows <- lapply(contributors_list, function(df) {
    sid <- unique(df$SampleID)
    tc  <- unique(df$TrueCommunity)
    cov_row <- stats::rnorm(P, mean = community_covariates[tc, ],
                            sd   = covariate_sds[tc, ])
    out <- data.frame(SampleID = sid, TrueCommunity = tc)
    out[colnames(community_covariates)] <- as.list(cov_row)
    out
  })
  do.call(rbind, rows)
}


#' Run the full eDNA survey simulation pipeline in one call
#'
#' @description
#' Convenience wrapper that runs all five simulation steps in sequence:
#' [generate_community_compositions()], [generate_contributors()],
#' [generate_eDNA()], [simulate_metabarcoding()], and
#' [generate_sample_covariates()]. Returns all intermediate outputs plus
#' a ready-to-use count matrix and covariate data frame formatted for
#' [eDNA_dmm()].
#'
#' @param n_communities Integer. Number of communities K. Default `4`.
#' @param n_species Integer. Number of species S. Default `40`.
#' @param samples_per_community Integer. Samples per community. Default `3`.
#' @param community_covariate_means A [K × P] matrix of covariate means per
#'   community. Rows = communities, columns = covariates. If `NULL` (default),
#'   a two-covariate (Depth × Distance_shore) design with four well-separated
#'   communities is used.
#' @param covariate_sds A [K × P] matrix of covariate SDs, or `NULL` for
#'   SD = 5 on all. Default `NULL`.
#' @param samples_per_community Integer. Default `3`.
#' @param mean_read_depth Numeric. Mean reads per sample. Default `10000`.
#' @param bio_reps Integer. Biological replicates per sample. Default `1`.
#' @param seq_reps Integer. Sequencing technical replicates per bio rep.
#'   Default `1`.
#' @param spillover Numeric [0,1]. Community spillover. Default `0.15`.
#' @param shedding_error Numeric. Lognormal shedding noise SD. Default `0.3`.
#' @param decay_rate Numeric. eDNA decay rate with distance. Default `0.1`.
#' @param seed Integer random seed. Default `42`.
#' @param ... Additional arguments passed to [generate_community_compositions()].
#'
#' @return A list with:
#' \describe{
#'   \item{`counts`}{Integer matrix [N × S] ready for [eDNA_dmm()].}
#'   \item{`covariates`}{Data frame with `sample_id`, `TrueCommunity`, and
#'     covariate columns, ready for [eDNA_dmm()].}
#'   \item{`community_compositions`}{The true [K × S] composition matrix.}
#'   \item{`metab_df`}{Raw metabarcoding data frame.}
#'   \item{`sample_metadata`}{Full sample metadata data frame.}
#'   \item{`contributors`}{Output of [generate_contributors()].}
#' }
#'
#' @seealso [eDNA_dmm()], [generate_community_compositions()],
#'   [plot_true_compositions()]
#'
#' @examples
#' \dontrun{
#' sim <- simulate_eDNA_survey(
#'   n_communities         = 3,
#'   n_species             = 30,
#'   samples_per_community = 5,
#'   seed                  = 42
#' )
#' dim(sim$counts)       # N x 30
#' head(sim$covariates)
#'
#' fit <- eDNA_dmm(sim$counts, sim$covariates[, c("Depth", "Distance_shore")], K = 3)
#' }
#'
#' @export
simulate_eDNA_survey <- function(
    n_communities            = 4,
    n_species                = 40,
    samples_per_community    = 3,
    community_covariate_means = NULL,
    covariate_sds            = NULL,
    mean_read_depth          = 10000,
    bio_reps                 = 1,
    seq_reps                 = 1,
    spillover                = 0.15,
    shedding_error           = 0.3,
    decay_rate               = 0.1,
    seed                     = 42,
    ...
) {
  set.seed(seed)

  # Default covariate design for 4 communities (2x2: depth x shore)
  if (is.null(community_covariate_means)) {
    if (n_communities == 4) {
      community_covariate_means <- matrix(
        c(80, 200,
          10, 200,
          80,  20,
          10,  20),
        nrow = 4, byrow = TRUE,
        dimnames = list(NULL, c("Depth", "Distance_shore"))
      )
    } else {
      # Generic: space communities evenly along depth axis
      depths  <- seq(10, 200, length.out = n_communities)
      community_covariate_means <- matrix(
        c(depths, rep(100, n_communities)),
        nrow = n_communities,
        dimnames = list(NULL, c("Depth", "Distance_shore"))
      )
    }
  }

  # 1. Compositions
  community_mat <- generate_community_compositions(
    n_communities = n_communities,
    n_species     = n_species,
    spillover     = spillover,
    seed          = NULL,   # already set globally
    ...
  )

  # 2. Contributors
  grouped_vec <- as.logical(sample(c(TRUE, FALSE), n_species, replace = TRUE))
  contrib_obj <- generate_contributors(
    community_compositions = community_mat,
    samples_per_community  = samples_per_community,
    grouped_species_vec    = grouped_vec,
    seed                   = NULL
  )

  # 3. eDNA
  eDNA_obj <- generate_eDNA(
    contributors_list = contrib_obj$contributors_list,
    shedding_error    = shedding_error,
    decay_rate        = decay_rate,
    shedding_rate     = 1000,
    bio_reps          = bio_reps
  )

  # 4. Metabarcoding
  amplification_bias <- stats::rlnorm(n_species, meanlog = 0, sdlog = 0.01)
  metab_df <- simulate_metabarcoding(
    eDNA_obj,
    mean_read_depth = mean_read_depth,
    rep             = seq_reps,
    amplification_bias = amplification_bias
  )

  # 5. Covariates
  sample_metadata <- generate_sample_covariates(
    contributors_list    = contrib_obj$contributors_list,
    community_covariates = community_covariate_means,
    covariate_sds        = covariate_sds
  )

  # 6. Build count matrix [N x S]
  sample_counts <- metab_df |>
    dplyr::group_by(.data$SampleID, .data$Species) |>
    dplyr::summarise(Counts = sum(.data$Counts), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "Species", values_from = "Counts",
                       values_fill = 0L) |>
    dplyr::arrange(.data$SampleID)

  X <- as.matrix(sample_counts[, -1])
  storage.mode(X) <- "integer"
  rownames(X) <- paste0("STN_", sprintf("%03d", sample_counts$SampleID))
  colnames(X) <- paste0("Sp_", colnames(X))

  # 7. Covariates data frame aligned to X
  covariates_df <- sample_metadata |>
    dplyr::arrange(.data$SampleID) |>
    dplyr::mutate(sample_id = paste0("STN_", sprintf("%03d", .data$SampleID)))

  list(
    counts                = X,
    covariates            = covariates_df,
    community_compositions = community_mat,
    metab_df              = metab_df,
    sample_metadata       = sample_metadata,
    contributors          = contrib_obj
  )
}
