haplotype_review_calls <- function(values) {
  bases <- c("A", "C", "G", "T")
  calls <- toupper(trimws(as.character(values)))
  missing <- is.na(calls) | calls %in% c("", "FAIL", "99", "N", "NA", "NOCALL", "NO CALL", "-", ".")
  first <- match(calls, bases)
  second <- first
  diploid <- !missing & grepl("^[ACGT]\\s*[/|*]\\s*[ACGT]$", calls)
  first[diploid] <- match(substr(calls[diploid], 1L, 1L), bases)
  second[diploid] <- match(substr(calls[diploid], nchar(calls[diploid]), nchar(calls[diploid])), bases)
  type <- rep("Unresolved", length(calls))
  type[missing] <- "Missing"
  type[!is.na(first) & first == second] <- "Homozygous"
  type[!is.na(first) & first != second] <- "Heterozygous"
  list(first = first, second = second, type = type)
}

haplotype_group_review <- function(data, assignments, qc, min_reference_lines = 2L,
                                   min_consistency = 0.80, min_reference_coverage = 0.50,
                                   min_supporting_markers = 2L, min_support = 0.60,
                                   min_margin = 0.15, min_comparable_coverage = 0.25) {
  settings <- list(
    min_reference_lines = min_reference_lines, min_consistency = min_consistency,
    min_reference_coverage = min_reference_coverage, min_supporting_markers = min_supporting_markers,
    min_support = min_support, min_margin = min_margin,
    min_comparable_coverage = min_comparable_coverage
  )
  proportions <- settings[c("min_consistency", "min_reference_coverage", "min_support", "min_margin",
                            "min_comparable_coverage")]
  valid_proportion <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0 && x <= 1
  positive_integer <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  if (!all(vapply(proportions, valid_proportion, logical(1))) ||
      !positive_integer(min_reference_lines) || !positive_integer(min_supporting_markers)) {
    stop("Group review needs positive whole count limits and proportions between zero and one.",
         call. = FALSE)
  }
  keys <- c("germplasmName", "Concatenated_Grouping", "Concatenated_Specific")
  if (!is.data.frame(data) || !"germplasmName" %in% names(data) ||
      !is.data.frame(assignments) || !all(keys %in% names(assignments)) ||
      !is.data.frame(qc) || !all(c("Marker", "Used_For_Splitting") %in% names(qc))) {
    stop("Group review needs raw SNP calls, original split assignments, and marker QC results.",
         call. = FALSE)
  }

  ids <- as.character(assignments$germplasmName)
  raw_ids <- as.character(data$germplasmName)
  raw_rows <- match(ids, raw_ids)
  original <- as.character(assignments$Concatenated_Grouping)
  specific <- as.character(assignments$Concatenated_Specific)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) || anyDuplicated(raw_ids) ||
      anyNA(raw_rows) || anyNA(original) || any(!nzchar(original))) {
    stop("Group review needs unique matching germplasm IDs and original groups.", call. = FALSE)
  }
  markers <- as.character(qc$Marker[qc$Used_For_Splitting %in% TRUE])
  if (anyDuplicated(markers) || !all(markers %in% names(data))) {
    stop("The QC-passed markers do not match the raw SNP columns.", call. = FALSE)
  }
  raw_calls <- matrix(character(), nrow = length(ids), ncol = length(markers),
                      dimnames = list(ids, markers))
  if (length(ids) && length(markers)) {
    raw_calls[] <- as.character(as.matrix(data[raw_rows, markers, drop = FALSE]))
  }

  groups <- sort(unique(original), method = "radix")
  group_index <- match(original, groups)
  flagged <- grepl("(^|--)U(--|$)", specific)
  # Missing paths also cannot supply reference calls.
  reference <- !flagged & !is.na(specific) & nzchar(specific)
  reference_count <- tabulate(group_index[reference], nbins = length(groups))
  signature <- matrix(NA_integer_, nrow = length(groups), ncol = length(markers))
  strength <- matrix(0, nrow = length(groups), ncol = length(markers))

  fit_signature <- function(rows) {
    parsed <- haplotype_review_calls(raw_calls[rows, , drop = FALSE])
    typed <- parsed$first
    typed[parsed$type != "Homozygous"] <- NA_integer_
    dim(typed) <- c(length(rows), length(markers))
    counts <- do.call(cbind, lapply(seq_len(4L), function(allele) {
      colSums(typed == allele, na.rm = TRUE)
    }))
    called <- rowSums(counts)
    modal <- max.col(counts, ties.method = "first")
    modal_count <- counts[cbind(seq_along(markers), modal)]
    coverage <- called / length(rows)
    consistency <- modal_count / pmax(called, 1L)
    usable <- called >= min_reference_lines & coverage >= min_reference_coverage &
      consistency >= min_consistency & rowSums(counts == modal_count) == 1L
    modal[!usable] <- NA_integer_
    list(allele = modal, strength = ifelse(usable, coverage * consistency, 0))
  }
  if (length(markers)) {
    for (group in which(reference_count >= min_reference_lines)) {
      fit <- fit_signature(which(reference & group_index == group))
      signature[group, ] <- fit$allele
      strength[group, ] <- fit$strength
    }
  }
  diagnostic <- vapply(seq_along(markers), function(j) {
    length(unique(stats::na.omit(signature[, j]))) > 1L
  }, logical(1))
  group_weights <- strength
  group_weights[, !diagnostic] <- 0
  viable <- which(rowSums(group_weights > 0) > 0)

  empty_summary <- data.frame(
    germplasmName = character(), Original_Group = character(), Suggested_Group = character(),
    Best_Group = character(), Runner_Up_Group = character(), Best_Score = numeric(),
    Runner_Up_Score = numeric(), Score_Margin = numeric(), Comparable_Markers = integer(),
    Comparable_Coverage = numeric(),
    Supporting_Markers = integer(), Conflicting_Markers = integer(), Homozygous = integer(),
    Missing = integer(), Explicit_Heterozygous = integer(), Unresolved = integer(),
    Status = character(), stringsAsFactors = FALSE
  )

  score_line <- function(row) {
    parsed <- haplotype_review_calls(raw_calls[row, ])
    summary <- list(
      germplasmName = ids[row], Original_Group = original[row], Suggested_Group = NA_character_,
      Best_Group = NA_character_, Runner_Up_Group = NA_character_, Best_Score = NA_real_,
      Runner_Up_Score = NA_real_, Score_Margin = NA_real_, Comparable_Markers = 0L,
      Comparable_Coverage = NA_real_,
      Supporting_Markers = 0L, Conflicting_Markers = 0L, Homozygous = sum(parsed$type == "Homozygous"),
      Missing = sum(parsed$type == "Missing"), Explicit_Heterozygous = sum(parsed$type == "Heterozygous"),
      Unresolved = sum(parsed$type == "Unresolved"),
      Status = if (length(markers)) "Insufficient reference groups" else "No QC-passed markers"
    )
    if (length(viable) < 2L) return(summary)

    original_supported <- group_index[row] %in% viable
    multiplier <- ifelse(parsed$type == "Homozygous", 1, ifelse(parsed$type == "Heterozygous", 0.5, 0))
    observed <- multiplier > 0
    match_first <- signature == rep(parsed$first, each = length(groups))
    match_second <- signature == rep(parsed$second, each = length(groups))
    hit <- match_first | match_second
    hit[is.na(hit)] <- FALSE
    matches <- sweep(hit, 2L, multiplier, `*`)
    overall <- rowSums(matches * group_weights) /
      rowSums(group_weights[, observed, drop = FALSE])
    ranked <- viable[is.finite(overall[viable])]
    if (!length(ranked)) {
      summary$Status <- if (original_supported) "No comparable called markers" else
        "Original group lacks reference support"
      return(summary)
    }
    best <- ranked[which.max(overall[ranked])]
    alternatives <- setdiff(viable, best)

    # Compare both groups on the same typed loci and the weaker signature weight.
    compare <- function(other) {
      weight <- pmin(group_weights[best, ], group_weights[other, ])
      distinguishing <- weight > 0 & signature[best, ] != signature[other, ]
      common <- observed & distinguishing
      if (!any(common)) {
        return(c(best = NA, other = NA, compared = 0, coverage = 0, support = 0, conflict = 0))
      }
      best_hits <- matches[best, common]
      other_hits <- matches[other, common]
      c(best = stats::weighted.mean(best_hits, weight[common]),
        other = stats::weighted.mean(other_hits, weight[common]), compared = sum(common),
        coverage = sum(common) / sum(distinguishing),
        support = sum(best_hits > other_hits), conflict = sum(best_hits < other_hits))
    }
    comparisons <- vapply(alternatives, compare, numeric(6L))
    margins <- comparisons["best", ] - comparisons["other", ]
    runner <- which.min(ifelse(is.na(margins), -Inf, margins))
    summary$Best_Group <- groups[best]
    summary$Runner_Up_Group <- groups[alternatives[runner]]
    summary$Best_Score <- comparisons["best", runner]
    summary$Runner_Up_Score <- comparisons["other", runner]
    summary$Score_Margin <- margins[runner]
    summary$Comparable_Markers <- as.integer(comparisons["compared", runner])
    summary$Comparable_Coverage <- comparisons["coverage", runner]
    summary$Supporting_Markers <- as.integer(comparisons["support", runner])
    summary$Conflicting_Markers <- as.integer(comparisons["conflict", runner])
    summary$Status <- if (!original_supported) {
      "Original group lacks reference support"
    } else if (anyNA(margins)) {
      "Insufficient comparable markers"
    } else if (any(comparisons["coverage", ] < min_comparable_coverage)) {
      "Too few distinguishing markers called"
    } else if (any(margins <= sqrt(.Machine$double.eps) | margins < min_margin)) {
      "Tie or small score margin"
    } else if (any(comparisons["support", ] < min_supporting_markers)) {
      "Too few supporting markers"
    } else if (any(comparisons["best", ] < min_support)) {
      "Low weighted support"
    } else {
      summary$Suggested_Group <- groups[best]
      if (identical(groups[best], original[row])) "Original group supported" else "Suggested"
    }
    summary
  }

  suggestions <- empty_summary
  if (any(flagged)) {
    suggestions <- do.call(rbind, lapply(which(flagged), function(row) {
      as.data.frame(score_line(row), stringsAsFactors = FALSE)
    }))
    rownames(suggestions) <- NULL
  }

  evidence <- function(germplasm_name) {
    if (length(germplasm_name) != 1L || is.na(germplasm_name) || !germplasm_name %in% ids) {
      stop("Choose a germplasm ID from this haplotype analysis.", call. = FALSE)
    }
    row <- match(germplasm_name, ids)
    review <- match(germplasm_name, suggestions$germplasmName)
    best <- if (is.na(review)) original[row] else suggestions$Best_Group[review]
    runner <- if (is.na(review)) NA_character_ else suggestions$Runner_Up_Group[review]
    bases <- c("A", "C", "G", "T")
    raw <- unname(raw_calls[row, ])
    parsed <- haplotype_review_calls(raw)
    best_codes <- signature[match(best, groups), ]
    runner_codes <- signature[match(runner, groups), ]
    label_match <- function(allele) {
      hit <- !is.na(allele) & !is.na(parsed$first) & (allele == parsed$first | allele == parsed$second)
      result <- rep("No evidence", length(markers))
      result[!is.na(allele) & !is.na(parsed$first)] <- "Mismatch"
      result[hit & parsed$type == "Homozygous"] <- "Match"
      result[hit & parsed$type == "Heterozygous"] <- "Partial"
      result
    }
    data.frame(
      Marker = markers, Raw_Call = raw, Call_Type = parsed$type,
      Best_Group = rep(best, length(markers)), Best_Allele = bases[best_codes],
      Best_Result = label_match(best_codes),
      Runner_Up_Group = rep(runner, length(markers)), Runner_Up_Allele = bases[runner_codes],
      Runner_Up_Result = label_match(runner_codes),
      Distinct_Signatures = !is.na(best_codes) & !is.na(runner_codes) & best_codes != runner_codes,
      stringsAsFactors = FALSE
    )
  }

  group_data <- data.frame(
    Group = groups, Lines = tabulate(group_index, nbins = length(groups)),
    Reference_Lines = reference_count, Consistent_Markers = rowSums(strength > 0),
    Diagnostic_Markers = rowSums(group_weights > 0), stringsAsFactors = FALSE
  )
  list(suggestions = suggestions, evidence = evidence, raw_calls = raw_calls,
       groups = group_data, settings = settings, markers = markers)
}
