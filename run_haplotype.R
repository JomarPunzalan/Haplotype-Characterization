# Run from RStudio with source("run_haplotype.R"), then call run_haplotype().
# In a terminal, use Rscript run_haplotype.R <snp-csv> [output-folder].

read_haplotype_input <- function(path, marker_pattern) {
  data <- utils::read.csv(
    path, colClasses = "character", check.names = FALSE,
    na.strings = "", row.names = NULL, fileEncoding = "UTF-8-BOM"
  )
  if (!nrow(data) || anyDuplicated(names(data)) || any(!nzchar(names(data)))) {
    stop("The CSV must have data rows and unique, nonempty column names.", call. = FALSE)
  }
  if (!"Line" %in% names(data) && "germplasmName" %in% names(data)) {
    names(data)[names(data) == "germplasmName"] <- "Line"
  }
  if (!"Line" %in% names(data)) {
    stop("The SNP CSV must contain a Line or germplasmName column.", call. = FALSE)
  }
  data$Line <- trimws(data$Line)
  if (anyNA(data$Line) || any(!nzchar(data$Line)) || anyDuplicated(data$Line)) {
    stop("Line IDs must be present and unique.", call. = FALSE)
  }

  markers <- setdiff(grep(marker_pattern, names(data), ignore.case = TRUE, value = TRUE), "Line")
  if (!length(markers)) {
    stop("No SNP columns match marker_pattern: ", marker_pattern, call. = FALSE)
  }
  calls <- as.matrix(data[markers])
  rownames(calls) <- data$Line
  calls[] <- toupper(trimws(calls))
  heterozygous <- !is.na(calls) & (calls == "1" | grepl("[/*]", calls))
  known_missing <- is.na(calls) | calls %in% c("", "FAIL", "99", "N", "NA", "NOCALL", "NO CALL", "-", ".")
  homozygous <- !is.na(calls) & calls %in% c("A", "C", "G", "T")
  unrecognized <- colSums(!homozygous & !heterozygous & !known_missing)
  calls[!homozygous & !heterozygous] <- NA_character_
  calls[heterozygous] <- "1"

  list(data = data, calls = calls, unrecognized = unrecognized)
}

haplotype_marker_qc <- function(calls, unrecognized, max_missing, max_heterozygous,
                                report_max_heterozygous) {
  bases <- c("A", "C", "G", "T")
  counts <- vapply(bases, function(base) colSums(calls == base, na.rm = TRUE), numeric(ncol(calls)))
  if (ncol(calls) == 1L) counts <- matrix(counts, nrow = 1L, dimnames = list(NULL, bases))
  called <- rowSums(counts)
  missing <- colSums(is.na(calls))
  het <- colSums(calls == "1", na.rm = TRUE)
  maf <- vapply(seq_len(nrow(counts)), function(i) {
    positive <- counts[i, counts[i, ] > 0]
    if (length(positive) > 1) return(min(positive) / called[i])
    if (length(positive) == 1) 0 else NA_real_
  }, numeric(1))
  missing_frequency <- missing / nrow(calls)
  het_frequency <- het / nrow(calls)
  data.frame(
    Marker = colnames(calls), counts, Called = called,
    Missing = missing, Missing_Frequency = missing_frequency,
    Het = het, Het_Frequency = het_frequency, Unrecognized = unrecognized, MAF = maf,
    Used_For_Splitting = called > 0 & missing_frequency <= max_missing &
      het_frequency <= max_heterozygous,
    Included_In_Filtered_Report = called > 0 & missing_frequency <= max_missing &
      het_frequency <= report_max_heterozygous,
    row.names = NULL, check.names = FALSE
  )
}

split_haplotypes <- function(calls, max_depth = 3L,
                             split_rule = c("highest_maf", "min_maf", "min_maf_fallback"),
                             min_maf = 0.35, fallback_rank = 5L) {
  split_rule <- match.arg(split_rule)
  positive_integer <- function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value) &&
      value >= 1 && value <= .Machine$integer.max && value == floor(value)
  }
  if (!positive_integer(max_depth)) {
    stop("Maximum split depth must be a positive whole number.", call. = FALSE)
  }
  if (!is.numeric(min_maf) || length(min_maf) != 1L || !is.finite(min_maf) ||
      min_maf < 0 || min_maf > 0.5) {
    stop("Minimum MAF must be between 0 and 0.5.", call. = FALSE)
  }
  if (!positive_integer(fallback_rank)) {
    stop("MAF fallback rank must be a positive whole number.", call. = FALSE)
  }

  alleles <- c("A", "C", "G", "T")
  if (!is.matrix(calls) || !is.character(calls) ||
      any(!calls %in% c(alleles, "1", NA_character_))) {
    stop("Expected a character SNP matrix containing A, C, G, T, 1, or NA.", call. = FALSE)
  }
  valid_names <- function(names, size) {
    length(names) == size && !anyNA(names) && all(nzchar(names)) && !anyDuplicated(names)
  }
  if (!valid_names(rownames(calls), nrow(calls)) || !valid_names(colnames(calls), ncol(calls))) {
    stop("The SNP matrix needs unique, nonempty line IDs and marker names.", call. = FALSE)
  }

  depths <- seq_len(max_depth)
  grouping <- matrix(0L, nrow(calls), max_depth)
  specific <- matrix("0", nrow(calls), max_depth)
  split_markers <- matrix(NA_character_, nrow(calls), max_depth)
  last_group <- integer(max_depth)
  records <- list()

  visit <- function(rows, unused, depth, parent) {
    if (depth > max_depth || length(rows) < 2L || !length(unused)) return(invisible(NULL))

    node_calls <- calls[rows, unused, drop = FALSE]
    counts <- do.call(cbind, lapply(alleles, function(allele) {
      colSums(node_calls == allele, na.rm = TRUE)
    }))
    eligible <- which(rowSums(counts > 0) == 2L)
    if (!length(eligible)) return(invisible(NULL))

    called <- rowSums(counts)[eligible]
    majority <- apply(counts[eligible, , drop = FALSE], 1L, max)
    maf <- (called - majority) / called
    if (split_rule == "highest_maf") {
      choice <- which.max(maf)
    } else {
      pool <- which(maf >= min_maf)
      if (!length(pool) && split_rule == "min_maf_fallback") {
        pool <- utils::head(order(maf, decreasing = TRUE, method = "radix"), fallback_rank)
      }
      if (!length(pool)) return(invisible(NULL))
      # The historical threshold option selects the fewest called alleles.
      choice <- pool[which.min(called[pool])]
    }

    position <- eligible[choice]
    marker <- unused[position]
    observed <- which(counts[position, ] > 0)
    # In a tie, the first allele in A/C/G/T order is the major allele.
    ranked_alleles <- observed[order(-counts[position, observed], observed)]
    major <- alleles[ranked_alleles[1L]]
    minor <- alleles[ranked_alleles[2L]]
    typed <- node_calls[, position]
    minor_rows <- typed %in% minor
    uncertain <- !typed %in% c(major, minor)
    group_ids <- last_group[depth] + c(1L, 2L)
    last_group[depth] <<- group_ids[2L]

    # Uncertain lines follow the major branch but remain visible as U.
    grouping[rows, depth] <<- ifelse(minor_rows, group_ids[2L], group_ids[1L])
    specific[rows, depth] <<- ifelse(uncertain, "U", as.character(grouping[rows, depth]))
    split_markers[rows, depth] <<- colnames(calls)[marker]
    records[[length(records) + 1L]] <<- data.frame(
      Depth = depth, Parent = parent, Marker = colnames(calls)[marker], MAF = maf[choice],
      Major = major, Minor = minor, N_Major = sum(typed %in% major),
      N_Minor = sum(minor_rows), N_Uncertain = sum(uncertain),
      row.names = NULL, stringsAsFactors = FALSE
    )

    child_paths <- if (depth == 1L) as.character(group_ids) else paste(parent, group_ids, sep = "--")
    remaining <- unused[-position]
    visit(rows[!minor_rows], remaining, depth + 1L, child_paths[1L])
    visit(rows[minor_rows], remaining, depth + 1L, child_paths[2L])
    invisible(NULL)
  }
  visit(seq_len(nrow(calls)), seq_len(ncol(calls)), 1L, "Root")

  # Present the log in depth order, preserving traversal order within a depth.
  splits <- if (length(records)) {
    out <- do.call(rbind, records)
    out[order(out$Depth), , drop = FALSE]
  } else {
    data.frame(Depth = integer(), Parent = character(), Marker = character(), MAF = numeric(),
               Major = character(), Minor = character(), N_Major = integer(),
               N_Minor = integer(), N_Uncertain = integer())
  }
  rownames(splits) <- NULL

  colnames(grouping) <- paste0("Split_Grouping", depths)
  colnames(split_markers) <- paste0("SNP_Split", depths)
  concatenate <- function(paths) {
    if (!nrow(paths)) return(character())
    apply(paths, 1L, paste, collapse = "--")
  }
  assignments <- data.frame(
    Line = if (nrow(calls)) rownames(calls) else character(), grouping, split_markers,
    Concatenated_Grouping = concatenate(grouping),
    Concatenated_Specific = concatenate(specific),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  list(assignments = assignments, splits = splits)
}

haplotype_allele_similarity <- function(consensus) {
  labels <- rownames(consensus)
  bases <- c("A", "C", "G", "T")
  homozygous <- matrix(
    consensus %in% bases, nrow = nrow(consensus), ncol = ncol(consensus)
  )
  unresolved <- !is.na(consensus) & !homozygous
  mixed <- tcrossprod(1 * homozygous, 1 * unresolved)
  similarity <- 0.5 * (mixed + t(mixed)) + 0.25 * tcrossprod(1 * unresolved)

  for (base in bases) {
    called <- 1 * (!is.na(consensus) & consensus == base)
    similarity <- similarity + tcrossprod(called)
  }
  similarity <- similarity / ncol(consensus)
  diag(similarity) <- NA_real_
  dimnames(similarity) <- list(labels, labels)
  similarity
}

haplotype_report <- function(data, calls, assignments, splits, qc, settings, workbook) {
  markers <- colnames(calls)
  ids <- as.character(assignments$Line)
  raw_rows <- match(ids, data$Line)
  call_rows <- match(ids, rownames(calls))
  if (!length(ids) || !length(markers) || anyNA(raw_rows) || anyNA(call_rows)) {
    stop("Reporting requires line assignments and matching raw and normalized SNP calls.", call. = FALSE)
  }

  data <- data[raw_rows, , drop = FALSE]
  calls <- calls[call_rows, , drop = FALSE]
  groups <- sort(unique(assignments$Concatenated_Grouping), method = "radix")
  group_ids <- match(assignments$Concatenated_Grouping, groups)
  membership <- split(seq_along(ids), factor(group_ids, levels = seq_along(groups)))
  group_summary <- data.frame(
    Grouping_ID = seq_along(groups),
    Grouping_Count = lengths(membership),
    Concatenated_Grouping = groups,
    row.names = NULL
  )

  marker_mode <- function(values) {
    counts <- table(values)
    if (!length(counts)) return(NA_character_)
    paste(names(counts)[counts == max(counts)], collapse = " and ")
  }
  consensus <- matrix(
    vapply(membership, function(rows) {
      vapply(seq_along(markers), function(col) marker_mode(calls[rows, col]), character(1))
    }, character(length(markers))),
    nrow = length(groups), byrow = TRUE, dimnames = list(groups, markers)
  )
  similarity <- haplotype_allele_similarity(consensus)

  split_markers <- grep("^SNP_Split[0-9]+$", names(assignments), value = TRUE)
  representatives <- data.frame(
    group_summary,
    assignments[match(groups, assignments$Concatenated_Grouping), split_markers, drop = FALSE],
    consensus,
    check.names = FALSE,
    row.names = NULL
  )

  metadata <- setdiff(names(data), c(names(assignments), "Grouping_ID", markers))
  line_assignments <- data.frame(
    Line = ids,
    Grouping_ID = group_ids,
    assignments[, setdiff(names(assignments), "Line"), drop = FALSE],
    data[, c(metadata, markers), drop = FALSE],
    check.names = FALSE,
    row.names = NULL
  )
  filtered <- markers[markers %in% qc$Marker[qc$Included_In_Filtered_Report %in% TRUE]]
  filtered_samples <- line_assignments[, c(setdiff(names(line_assignments), markers), filtered), drop = FALSE]
  filtered_representatives <- representatives[, c(setdiff(names(representatives), markers), filtered), drop = FALSE]

  preliminary <- unique(splits[, c("Marker", "Depth"), drop = FALSE])
  preliminary_markers <- data.frame(
    Number = match(preliminary$Marker, unique(preliminary$Marker)),
    Preliminary_SNP_Set = preliminary$Marker,
    Split_Included_In = sprintf("SNP_Split%d", preliminary$Depth)
  )
  description <- data.frame(
    Setting = as.character(names(settings)),
    Value = vapply(settings, function(value) paste(value, collapse = ", "), character(1)),
    row.names = NULL
  )
  workflow <- data.frame(
    Stage = c("SNP quality", "Splits", "Group names", "Representatives", "Allele similarity"),
    Description = c(
      "Quality metrics use normalized calls. Sample sheets retain the original SNP values.",
      "The split decisions identify the selected SNP and allele counts for each branch.",
      "A zero marks an unused split. U marks an uncertain SNP call at a split.",
      "Each group has one representative. Missing calls are excluded from consensus and ties are retained.",
      paste(
        "Across the unfiltered representative SNPs, identical homozygous calls contribute 1,",
        "different homozygous calls contribute 0, an unresolved call with a homozygous call contributes 0.5,",
        "two unresolved calls contribute 0.25, and a missing call contributes 0.",
        "The sum is divided by the number of SNPs. The diagonal is blank."
      )
    )
  )

  sheets <- list(
    "Description" = description,
    "Workflow" = workflow,
    "Preliminary SNP Set" = preliminary_markers,
    "Split Decisions" = splits,
    "Main Data Unfiltered SNPs" = line_assignments,
    "Main Data Filtered SNPs" = filtered_samples,
    "Rep_HaplotypeUnfiltered" = representatives,
    "Rep_HaplotypeFiltered" = filtered_representatives,
    "SNP Frequency" = qc
  )
  if ("general_grouping" %in% names(data)) {
    known <- as.character(data$general_grouping)
    known[is.na(known) | !nzchar(trimws(known))] <- "(Unassigned)"
    counts <- table(factor(group_ids, levels = seq_along(groups)), known)
    sheets[["Cross-Tabulation"]] <- data.frame(
      group_summary, as.data.frame.matrix(counts), check.names = FALSE, row.names = NULL
    )
  }
  sheets[["Allele Similarity"]] <- data.frame(
    Concatenated_Grouping = groups, similarity, check.names = FALSE, row.names = NULL
  )

  wb <- openxlsx::createWorkbook()
  header <- openxlsx::createStyle(
    fgFill = "#461D7C", fontColour = "#FFFFFF", textDecoration = "bold", wrapText = TRUE
  )
  split_header <- openxlsx::createStyle(
    fgFill = "#FDD023", fontColour = "#461D7C", textDecoration = "bold"
  )
  allele_styles <- list(
    A = openxlsx::createStyle(fontColour = "#9C6500", bgFill = "#FFEB9C"),
    C = openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE"),
    G = openxlsx::createStyle(fontColour = "#000080", bgFill = "#B8CCE4"),
    T = openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
  )
  uncertain <- openxlsx::createStyle(fontColour = "#FFFFFF", bgFill = "#3E424B")
  marker_sheets <- c(
    "Main Data Unfiltered SNPs", "Main Data Filtered SNPs",
    "Rep_HaplotypeUnfiltered", "Rep_HaplotypeFiltered"
  )

  for (sheet in names(sheets)) {
    contents <- sheets[[sheet]]
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, contents, headerStyle = header)
    openxlsx::freezePane(wb, sheet, firstRow = TRUE, firstCol = TRUE)
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(contents)), widths = 18)
    if (sheet %in% marker_sheets) {
      cols <- which(names(contents) %in% markers)
      selected_cols <- which(names(contents) %in% splits$Marker)
      if (length(selected_cols)) {
        openxlsx::addStyle(wb, sheet, split_header, rows = 1, cols = selected_cols, stack = TRUE)
      }
      if (length(cols) && nrow(contents)) {
        rows <- seq_len(nrow(contents)) + 1L
        cell <- paste0("UPPER(TRIM(", openxlsx::int2col(cols[1L]), "2))")
        rules <- sprintf('%s="%s"', cell, names(allele_styles))
        for (i in seq_along(rules)) {
          openxlsx::conditionalFormatting(
            wb, sheet, cols, rows, rule = rules[i], style = allele_styles[[i]], type = "expression"
          )
        }
        openxlsx::conditionalFormatting(
          wb, sheet, cols, rows,
          rule = sprintf('AND(%s<>"",NOT(OR(%s)))', cell, paste(rules, collapse = ",")),
          style = uncertain, type = "expression"
        )
      }
    }
  }

  if ("Cross-Tabulation" %in% names(sheets)) {
    cross_tab <- sheets[["Cross-Tabulation"]]
    for (col in seq.int(ncol(group_summary) + 1L, ncol(cross_tab))) {
      if (length(unique(cross_tab[[col]])) > 1L) {
        openxlsx::conditionalFormatting(
          wb, "Cross-Tabulation", cols = col, rows = seq_len(nrow(cross_tab)) + 1L,
          type = "colorScale", style = c("#F8696B", "#FEEB84", "#63BE7B")
        )
      }
    }
  }

  openxlsx::setColWidths(wb, "Description", cols = 1:2, widths = c(26, 70))
  openxlsx::setColWidths(wb, "Workflow", cols = 1:2, widths = c(24, 95))
  openxlsx::addStyle(
    wb, "Workflow", openxlsx::createStyle(wrapText = TRUE, valign = "top"),
    rows = seq_len(nrow(workflow)) + 1L, cols = 2
  )
  openxlsx::setRowHeights(wb, "Workflow", rows = seq_len(nrow(workflow)) + 1L,
                         heights = c(30, 30, 30, 30, 65))

  openxlsx::saveWorkbook(wb, workbook, overwrite = TRUE)
  list(line_assignments = line_assignments, representatives = representatives, similarity = similarity)
}

run_haplotype <- function(
    input_csv, output_dir = file.path(dirname(input_csv), "haplotype_results"),
    max_depth = 3L, max_missing = 0.10, max_heterozygous = 0.10,
    report_max_heterozygous = 0.20, split_rule = "highest_maf",
    min_maf = 0.35, fallback_rank = 5L, marker_pattern = "^SNP_chr",
    overwrite = FALSE) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop('Install the workbook package once with install.packages("openxlsx").', call. = FALSE)
  }
  split_rule <- match.arg(split_rule, c("highest_maf", "min_maf", "min_maf_fallback"))
  rates <- list(max_missing = max_missing, max_heterozygous = max_heterozygous,
                report_max_heterozygous = report_max_heterozygous, min_maf = min_maf)
  valid_rates <- vapply(rates, function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value) && value >= 0 && value <= 1
  }, logical(1))
  if (!all(valid_rates)) {
    stop("Enter values between 0 and 1 for: ", paste(names(rates)[!valid_rates], collapse = ", "),
         call. = FALSE)
  }

  files <- c(workbook = "output.xlsx", line_assignments = "final_merged_with_specific2.csv",
             representatives = "filtered_commonsnps.csv", similarity = "grm.csv")
  paths <- stats::setNames(file.path(output_dir, files), names(files))
  existing <- file.exists(paths)
  if (any(existing) && !isTRUE(overwrite)) {
    stop("Results already exist in ", output_dir,
         ". Choose another output folder or set overwrite = TRUE in R.", call. = FALSE)
  }
  if (normalizePath(input_csv, mustWork = TRUE) %in% normalizePath(paths[existing])) {
    stop("Choose an output folder that keeps the input CSV separate from the results.", call. = FALSE)
  }

  input <- read_haplotype_input(input_csv, marker_pattern)
  qc <- haplotype_marker_qc(input$calls, input$unrecognized, max_missing,
                            max_heterozygous, report_max_heterozygous)
  split <- split_haplotypes(input$calls[, qc$Used_For_Splitting, drop = FALSE],
                           max_depth, split_rule, min_maf, fallback_rank)
  settings <- c(list(
    input = basename(input_csv), lines = nrow(input$data), markers = nrow(qc),
    max_depth = max_depth, split_rule = split_rule,
    fallback_rank = fallback_rank, marker_pattern = marker_pattern
  ), rates)

  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop("Could not create the output folder: ", output_dir, call. = FALSE)
  }
  staging <- tempfile(".haplotype-", tmpdir = output_dir)
  if (!dir.create(staging)) stop("Could not stage the output files.", call. = FALSE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)

  report <- haplotype_report(input$data, input$calls, split$assignments, split$splits, qc,
                             settings, file.path(staging, files[["workbook"]]))
  utils::write.csv(report$line_assignments, file.path(staging, files[["line_assignments"]]),
                   row.names = FALSE, na = "")
  utils::write.csv(report$representatives, file.path(staging, files[["representatives"]]),
                   row.names = FALSE, na = "")
  utils::write.csv(report$similarity, file.path(staging, files[["similarity"]]), na = "")

  copied <- file.copy(file.path(staging, files), paths, overwrite = isTRUE(overwrite))
  if (!all(copied)) {
    stop("Could not save: ", paste(paths[!copied], collapse = ", "), call. = FALSE)
  }
  message("Wrote the haplotype workbook and CSVs to ", normalizePath(output_dir))
  invisible(paths)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 1L && args[1L] %in% c("-h", "--help")) {
    cat("Usage: Rscript run_haplotype.R <snp-csv> [output-folder]\n")
  } else if (!length(args) || length(args) > 2L) {
    stop("Usage: Rscript run_haplotype.R <snp-csv> [output-folder]", call. = FALSE)
  } else {
    do.call(run_haplotype, as.list(args))
  }
}
