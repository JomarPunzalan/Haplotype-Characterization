group_review_umap <- function(calls, review, min_comparable = 20L, seed = 2026L) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop('Install the optional uwot package with install.packages("uwot") to build the UMAP.', call. = FALSE)
  }
  required <- c("germplasmName", "Group", "Needs_Review")
  if (!is.matrix(calls) || !is.character(calls) || is.null(rownames(calls)) ||
      anyDuplicated(rownames(calls)) || !is.data.frame(review) ||
      !all(required %in% names(review))) {
    stop("UMAP requires named SNP calls and the current group review table.", call. = FALSE)
  }
  if (!is.numeric(min_comparable) || length(min_comparable) != 1L ||
      !is.finite(min_comparable) || min_comparable < 1 || min_comparable != floor(min_comparable)) {
    stop("Minimum comparable SNPs must be a positive whole number.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      abs(seed) > .Machine$integer.max || seed != floor(seed)) {
    stop("UMAP seed must be a whole number.", call. = FALSE)
  }
  ids <- as.character(review$germplasmName)
  rows <- match(ids, rownames(calls))
  if (anyNA(ids) || anyNA(rows) || anyDuplicated(ids)) {
    stop("Each review germplasmName needs one matching raw SNP row.", call. = FALSE)
  }

  varies <- function(first, second) {
    pattern <- 4L * pmin(first, second) + pmax(first, second)
    length(unique(pattern)) > 1L
  }
  comparable <- integer(length(ids))
  available <- logical(ncol(calls))
  for (col in seq_len(ncol(calls))) {
    decoded <- haplotype_review_calls(calls[rows, col])
    known <- !is.na(decoded$first) & !is.na(decoded$second)
    if (sum(known) >= 2L) {
      available[col] <- varies(decoded$first[known], decoded$second[known])
      comparable[known] <- comparable[known] + 1L
    }
  }
  included <- comparable >= min_comparable
  if (sum(included) < 3L) {
    stop(sprintf("UMAP needs at least three lines with %d comparable SNPs. Only %d meet that threshold.",
                 min_comparable, sum(included)), call. = FALSE)
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(previous_kind))
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  markers <- which(available)
  possible_features <- length(markers) * 4L
  # Keep the feature matrix bounded for large uploads by projecting allele
  # features into 64 signed bins. Smaller inputs use the original features.
  projected <- possible_features > 64L && as.double(sum(included)) * possible_features > 8e6
  if (projected) {
    features <- matrix(0, sum(included), 64L)
    bins <- sample.int(ncol(features), size = possible_features, replace = TRUE)
    signs <- sample(c(-1, 1), size = possible_features, replace = TRUE)
  } else {
    features <- matrix(0, sum(included), possible_features)
  }

  eligible_rows <- which(included)
  observed_variants <- logical(length(eligible_rows))
  retained <- 0L
  used_columns <- 0L
  for (index in seq_along(markers)) {
    decoded <- haplotype_review_calls(calls[rows[eligible_rows], markers[index]])
    known <- !is.na(decoded$first) & !is.na(decoded$second)
    first <- decoded$first[known]
    second <- decoded$second[known]
    if (!varies(first, second)) next

    observed_variants <- observed_variants | known
    retained <- retained + 1L
    for (allele in sort(unique(c(first, second)))) {
      dosage <- 0.5 * ((first == allele) + (second == allele))
      if (all(dosage == dosage[1L])) next
      # Mean imputation becomes zero after centering each allele's dosage.
      dosage <- dosage - mean(dosage)
      if (projected) {
        feature <- (index - 1L) * 4L + allele
        bin <- bins[feature]
        features[known, bin] <- features[known, bin] + signs[feature] * dosage
      } else {
        used_columns <- used_columns + 1L
        features[known, used_columns] <- dosage
      }
    }
  }
  if (retained == 0L) {
    stop("No comparable SNPs vary among the lines eligible for this UMAP.", call. = FALSE)
  }

  if (sum(observed_variants) < 3L) {
    stop("UMAP needs at least three eligible lines with observed calls at variable SNPs.", call. = FALSE)
  }
  if (!all(observed_variants)) features <- features[observed_variants, , drop = FALSE]
  if (!projected) features <- features[, seq_len(used_columns), drop = FALSE]
  if (!any(features != 0)) {
    stop("The SNP features have no variation for UMAP.", call. = FALSE)
  }
  dimensions <- if (!projected && ncol(features) > 50L && nrow(features) > 51L) 50L
  coordinates <- uwot::umap(
    features, n_neighbors = min(15L, nrow(features) - 1L), n_components = 2L,
    metric = "euclidean", init = "random", nn_method = "annoy", pca = dimensions,
    n_threads = 1L, n_sgd_threads = 1L, verbose = FALSE
  )
  if (any(!is.finite(coordinates))) {
    stop("UMAP did not produce finite coordinates for the comparable SNPs.", call. = FALSE)
  }

  points <- review[eligible_rows[observed_variants], required, drop = FALSE]
  points$UMAP1 <- coordinates[, 1L]
  points$UMAP2 <- coordinates[, 2L]
  rownames(points) <- NULL
  list(points = points, excluded = nrow(review) - nrow(points), markers = retained)
}
