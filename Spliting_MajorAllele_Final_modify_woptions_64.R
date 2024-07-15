library(dplyr)
rm(list=ls())

setwd('C:/Users/JPunzalan/OneDrive - LSU AgCenter/Desktop/Final_Haplotype')

# Function to get the next folder name
get_next_folder_name <- function(prefix) {
  i <- 1
  while (dir.exists(paste0(prefix, i))) {
    i <- i + 1
  }
  return(paste0(prefix, i))
}

# Read the SNP data from the file
snps <- read.csv("snp_haplotype2_525.csv")

# Check the structure of the snps data frame
str(snps)

# Replace "FAIL" with 99 in the snps data frame
snps[snps == "FAIL"] <- 99

# Filter columns where the frequency of cells with 99 is greater than 10%
threshold <- 0.10
columns_to_keep <- sapply(snps, function(col) mean(col == 99, na.rm = TRUE) <= threshold)
snps <- snps[, columns_to_keep]

# Remove columns where more than 20% of the cells contain heterozygous genotypes (identified by `*/*` or `/` inside cells)
hets_threshold <- 0.10
columns_to_keep <- sapply(snps, function(col) mean(grepl("/|\\*", col), na.rm = TRUE) <= hets_threshold)
snps <- snps[, columns_to_keep]

# Save the filtered SNP data to a CSV file
write.csv(snps, "snps_filtered.csv", row.names = FALSE)


# Function to calculate the Minor Allele Frequency (MAF)
calculate_maf <- function(genotype_vector) {
  valid_genotypes <- genotype_vector[genotype_vector %in% c("A", "T", "C", "G")]
  if (length(valid_genotypes) == 0) return(1) # If no valid genotypes, return highest MAF
  allele_counts <- table(valid_genotypes)
  allele_frequencies <- allele_counts / sum(allele_counts)
  minor_allele_frequency <- min(allele_frequencies)
  return(minor_allele_frequency)
}

# Function to find the major and minor alleles in a column
find_major_minor_alleles <- function(genotype_vector) {
  valid_genotypes <- genotype_vector[genotype_vector %in% c("A", "T", "C", "G")]
  allele_counts <- table(valid_genotypes)
  major_allele <- names(allele_counts)[which.max(allele_counts)]
  minor_allele <- names(allele_counts)[which.min(allele_counts)]
  return(list(major = major_allele, minor = minor_allele))
}

# Function to split the data based on the highest MAF
split_samples_by_maf <- function(snps, current_groups) {
  highest_maf <- 0
  split_column <- NA
  
  for (col in 2:ncol(snps)) { # Skip the first column which is 'Line'
    maf <- calculate_maf(snps[, col])
    if (maf > highest_maf && maf < 1) { # Ensure MAF is not 1
      highest_maf <- maf
      split_column <- col
    }
  }
  
  if (is.na(split_column)) {
    return(list(groups = current_groups, split_col_name = NA))
  }
  
  # Determine the major and minor alleles for the selected SNP column
  alleles <- find_major_minor_alleles(snps[, split_column])
  major_allele <- alleles$major
  minor_allele <- alleles$minor
  
  # Treat all heterozygous genotypes and NAs as major allele
  major_allele <- c(major_allele, "1", "99", NA)
  
  # Debugging: Print the selected column and alleles
  cat("Splitting based on column:", colnames(snps)[split_column], "\n")
  cat("Major allele(s):", paste(major_allele, collapse = ", "), "\n")
  cat("Minor allele:", minor_allele, "\n")
  
  # Split the samples into two groups based on the major and minor alleles
  group1 <- which(snps[, split_column] %in% major_allele | is.na(snps[, split_column]) | grepl("/", snps[, split_column]))
  group2 <- which(snps[, split_column] == minor_allele)
  
  # Debugging: Print the number of samples in each group
  cat("Group 1 size:", length(group1), "\n")
  cat("Group 2 size:", length(group2), "\n")
  
  # Update current groups
  current_groups[group1] <- 1
  current_groups[group2] <- 2
  
  return(list(groups = current_groups, split_col_name = colnames(snps)[split_column]))
}

# Function to split the data based on MAF with ranking and minimum MAF consideration
split_samples_by_maf_rank <- function(snps, current_groups, min_maf, rank_to_consider, option_2) {
  selected_maf_snps <- list()
  maf_values <- numeric()
  
  for (col in 2:ncol(snps)) { # Skip the first column which is 'Line'
    maf <- calculate_maf(snps[, col])
    maf_values <- c(maf_values, maf)
    if (maf >= min_maf && maf < 1) { # Ensure MAF is at least the specified minimum and not 1
      total_allele_notconsidered <- sum(snps[, col] != "99" & snps[, col] != "1")
      selected_maf_snps[[as.character(col)]] <- total_allele_notconsidered
    }
  }
  
  if (length(selected_maf_snps) == 0 && option_2 == "2b") {
    # Order the MAF values in decreasing order and exclude MAF == 1
    maf_values[maf_values == 1] <- NA
    top_maf_indices <- order(maf_values, decreasing = TRUE, na.last = NA)[1:rank_to_consider]
    selected_maf_snps <- list()
    for (index in top_maf_indices) {
      if (!is.na(maf_values[index])) {
        total_allele_notconsidered <- sum(snps[, index] != "99" & snps[, index] != "1")
        selected_maf_snps[[as.character(index)]] <- total_allele_notconsidered
      }
    }
  }
  
  if (length(selected_maf_snps) == 0) {
    return(list(groups = current_groups, split_col_name = NA))
  }
  
  # Convert the list to a numeric vector for which.min()
  selected_maf_snps_numeric <- unlist(selected_maf_snps)
  
  # Select the SNP with the lowest total_allele_notconsidered
  split_column <- as.numeric(names(which.min(selected_maf_snps_numeric)))
  
  # Determine the major and minor alleles for the selected SNP column
  alleles <- find_major_minor_alleles(snps[, split_column])
  major_allele <- alleles$major
  minor_allele <- alleles$minor
  
  # Treat all heterozygous genotypes and NAs as major allele
  major_allele <- c(major_allele, "1", "99", NA)
  
  # Debugging: Print the selected column and alleles
  cat("Splitting based on column:", colnames(snps)[split_column], "\n")
  cat("Major allele(s):", paste(major_allele, collapse = ", "), "\n")
  cat("Minor allele:", minor_allele, "\n")
  
  # Split the samples into two groups based on the major and minor alleles
  group1 <- which(snps[, split_column] %in% major_allele | is.na(snps[, split_column]) | grepl("/", snps[, split_column]))
  group2 <- which(snps[, split_column] == minor_allele)
  
  # Debugging: Print the number of samples in each group
  cat("Group 1 size:", length(group1), "\n")
  cat("Group 2 size:", length(group2), "\n")
  
  # Update current groups
  current_groups[group1] <- 1
  current_groups[group2] <- 2
  
  return(list(groups = current_groups, split_col_name = colnames(snps)[split_column]))
}

# Function to set target levels for splitting
set_target_levels <- function() {
  cat("Set the levels for MAF, Rank, and Options:\n")
  
  # Prompt the user to choose Option 1 or 2
  cat("Choose Option 1 or 2 (1: highest MAF excluding 1, 2: at least MAF and lowest total_allele_notconsidered): ")
  option_1 <- readLines(n = 1)
  
  if (option_1 == "2") {
    # Prompt the user to set the MAF level
    cat("Enter the MAF threshold (e.g., 0.35): ")
    maf_threshold <- as.numeric(readLines(n = 1))
    
    # Prompt the user to set the rank
    cat("Enter the Rank (number of top SNPs to consider, e.g., 5): ")
    rank <- as.integer(readLines(n = 1))
    
    # Prompt the user to choose Option 2a or 2b
    cat("Choose Option 2a or 2b (2a: stock setting, 2b: top rank with lowest total_allele_notconsidered if no SNP for stock setting): ")
    option_2 <- readLines(n = 1)
    
    return(list(maf_threshold = maf_threshold, rank = rank, option_1 = option_1, option_2 = option_2))
  } else {
    return(list(option_1 = option_1))
  }
}

# Main script execution
options <- set_target_levels()

# Step 4: Get the next folder name and create the directory
output_folder <- get_next_folder_name("repository_snp_")
dir.create(output_folder)

# Save the initial split result in the specified folder
original_working_directory <- getwd()
setwd(output_folder)
write.csv(snps, file = "snps_filtered.csv", row.names = FALSE)

# Function to perform a split on the provided SNP data and save the results
perform_split_and_save <- function(snps, group_tracking, split_level, parent_label, options, output_folder) {
  if (options$option_1 == "1") {
    result <- split_samples_by_maf(snps, rep(1, nrow(snps)))
  } else if (options$option_1 == "2") {
    result <- split_samples_by_maf_rank(snps, rep(1, nrow(snps)), options$maf_threshold, options$rank, options$option_2)
  } else {
    stop("Invalid option selected.")
  }
  
  split_column_name <- result$split_col_name
  
  if (is.na(split_column_name)) {
    return(NULL) # Stop if no split is possible
  }
  
  # Update the group labels
  group_tracking <- data.frame(SampleID = snps$Line, Group = result$groups)
  
  # Save the result to a CSV file
  file_name <- paste0(output_folder, "/", parent_label, "_Split", split_level, "_", split_column_name, ".csv")
  write.csv(group_tracking, file = file_name, row.names = FALSE)
  
  # Extract indices for each group
  group1_indices <- which(group_tracking$Group == 1)
  group2_indices <- which(group_tracking$Group == 2)
  
  # Extract the SNP data for each group
  group1_snps <- snps[group1_indices, ]
  group2_snps <- snps[group2_indices, ]
  
  # Remove the SNP column used for the split
  group1_snps <- group1_snps[, !colnames(group1_snps) %in% split_column_name]
  group2_snps <- group2_snps[, !colnames(group2_snps) %in% split_column_name]
  
  # Save the SNP data for each group
  group1_file_name <- paste0(output_folder, "/", parent_label, "_Split", split_level, "_Group1_", split_column_name, ".csv")
  group2_file_name <- paste0(output_folder, "/", parent_label, "_Split", split_level, "_Group2_", split_column_name, ".csv")
  write.csv(group1_snps, file = group1_file_name, row.names = FALSE)
  write.csv(group2_snps, file = group2_file_name, row.names = FALSE)
  
  return(list(
    group1_snps = group1_snps,
    group2_snps = group2_snps,
    group1_label = paste0(parent_label, "_Split", split_level, "_Group1"),
    group2_label = paste0(parent_label, "_Split", split_level, "_Group2")
  ))
}

# Step 1: Initializing the tracking data frame
group_tracking <- data.frame(SampleID = snps$Line, Split1 = rep(1, nrow(snps)))

# Step 2: Perform the first split
result <- if (options$option_1 == "1") {
  split_samples_by_maf(snps, group_tracking$Split1)
} else {
  split_samples_by_maf_rank(snps, group_tracking$Split1, options$maf_threshold, options$rank, options$option_2)
}

# Update the group_tracking with the result
group_tracking$Split1 <- result$groups
split_column_name <- result$split_col_name

# Step 3: Rename the columns to reflect the split
colnames(group_tracking)[2] <- paste0("Split1_", split_column_name)

# Step 4: Save the initial split result in the specified folder
write.csv(snps, file = paste0("Split1_", split_column_name, ".csv"), row.names = FALSE)

# Extract the indices for each group
group1_indices <- which(group_tracking$Split1 == 1)
group2_indices <- which(group_tracking$Split1 == 2)

# Extract the SNP data for each group
group1_snps <- snps[group1_indices, ]
group2_snps <- snps[group2_indices, ]

# Remove the SNP column used for the split
group1_snps <- group1_snps[, !colnames(group1_snps) %in% split_column_name]
group2_snps <- group2_snps[, !colnames(group2_snps) %in% split_column_name]

# Save the SNP data for each group
write.csv(group1_snps, file = paste0("Split1_Group1_", split_column_name, ".csv"), row.names = FALSE)
write.csv(group2_snps, file = paste0("Split1_Group2_", split_column_name, ".csv"), row.names = FALSE)

# Return to the original working directory
setwd(original_working_directory)

# Continue splitting the initial groups from Split1
split_level <- 2

# Perform the first split and save the groups if they exist
if (exists("group1_snps") && exists("group2_snps")) {
  group1_results <- perform_split_and_save(group1_snps, group_tracking[group1_indices, ], split_level, "Split1_Group1", options, output_folder)
  group2_results <- perform_split_and_save(group2_snps, group_tracking[group2_indices, ], split_level, "Split1_Group2", options, output_folder)
}

# Recursively split each group further
split_recursive <- function(snps, group_tracking, split_level, parent_label, max_splits, options, output_folder) {
  if (split_level > max_splits) {
    return(NULL)
  }
  
  result <- perform_split_and_save(snps, group_tracking, split_level, parent_label, options, output_folder)
  
  if (is.null(result)) {
    return(NULL) # Stop if no further split is possible
  }
  
  # Recursively split each subgroup
  split_recursive(result$group1_snps, group_tracking[group_tracking$Group == 1, ], split_level + 1, result$group1_label, max_splits, options, output_folder)
  split_recursive(result$group2_snps, group_tracking[group_tracking$Group == 2, ], split_level + 1, result$group2_label, max_splits, options, output_folder)
}

# Set the maximum number of splits
max_splits <- 3

if (!is.null(group1_results)) {
  split_recursive(group1_results$group1_snps, group_tracking[group1_indices, ], split_level + 1, group1_results$group1_label, max_splits, options, output_folder)
  split_recursive(group1_results$group2_snps, group_tracking[group1_indices, ], split_level + 1, group1_results$group2_label, max_splits, options, output_folder)
}

if (!is.null(group2_results)) {
  split_recursive(group2_results$group1_snps, group_tracking[group2_indices, ], split_level + 1, group2_results$group1_label, max_splits, options, output_folder)
  split_recursive(group2_results$group2_snps, group_tracking[group2_indices, ], split_level + 1, group2_results$group2_label, max_splits, options, output_folder)
}

# Print the final group tracking data frame for verification
print(group_tracking)

