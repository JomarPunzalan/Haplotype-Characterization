library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)

rm(list = ls())

setwd('C:/Users/JPunzalan/OneDrive - LSU AgCenter/Desktop/LSU_Haplotypes_Jen_64/repository_snp_6/Naming')

# Read the CSV files
split1 <- read.csv("grouping_snp_split1.csv")
split2 <- read.csv("grouping_snp_split2.csv")
split2 <- split2[,-c(4:ncol(split2))]
split3 <- read.csv("grouping_snp_split3.csv")
split3 <- split3[,-c(4:ncol(split3))]

# Process split1
split_1 <- split1[, 1:3]
split_1removed <- split1[, -c(2:3)]
snps <- split_1removed[, -1]

for (col in 1:ncol(snps)) {
  for (row in 1:nrow(snps)) {
    # Check if the value is not one of "A", "T", "C", "G", or NA
    if (!(snps[row, col] %in% c("A", "T", "C", "G", "FAIL"))) {
      # Replace with "Hets"
      snps[row, col] <- "1"
    }
  }
}

snps[snps == "FAIL"] <- 99

# Create a data frame to store the frequencies
frequency_df <- data.frame(
  Column = names(snps),
  A = numeric(ncol(snps)),
  T = numeric(ncol(snps)),
  C = numeric(ncol(snps)),
  G = numeric(ncol(snps)),
  Total_excludinghetsmissing = numeric(ncol(snps)),
  MAF = numeric(ncol(snps)),
  Het = numeric(ncol(snps)),
  Het_Frequency = numeric(ncol(snps)),
  Missing = numeric(ncol(snps))
)

# Calculate frequencies
for (col in seq_along(snps)) {
  total <- sum(!is.na(snps[, col]))
  total_allele <- sum(snps[, col] != "99" & snps[, col] != "1" & !is.na(snps[, col]))
  
  frequency_df$A[col] <- sum(snps[, col] == "A")
  frequency_df$T[col] <- sum(snps[, col] == "T")
  frequency_df$C[col] <- sum(snps[, col] == "C")
  frequency_df$G[col] <- sum(snps[, col] == "G")
  frequency_df$Total_excludinghetsmissing[col] <- total_allele
  frequency_df$Het[col] <- sum(snps[, col] == "1")
  frequency_df$Het_Frequency[col] <- frequency_df$Het[col] / total
  frequency_df$Missing[col] <- sum(snps[, col] == "99" | is.na(snps[, col]))
  
  total_alleles <- 2 * total_allele
  allele_frequencies <- c(2 * frequency_df$A[col], 2 * frequency_df$T[col], 2 * frequency_df$C[col], 2 * frequency_df$G[col]) / total_alleles
  frequency_df$MAF[col] <- min(allele_frequencies[allele_frequencies > 0])
}

# Transpose and prepare split_1removed
transpose <- as.data.frame(t(frequency_df))
transpose <- cbind(Row.Label = rownames(transpose), transpose)
rownames(transpose) <- NULL
colnames(transpose) <- transpose[1,]
transpose <- transpose[-1,]
names(transpose)[names(transpose) == 'Column'] <- 'Line'
split_1removed <- rbind(split_1removed, transpose)

# Extract SNP parts
extract_snp_part <- function(filename, pattern) {
  matches <- regmatches(filename, regexpr(pattern, filename))
  return(matches)
}

split_1$SNP_Split1 <- vapply(split_1$FileName_SNP_Grouping1, extract_snp_part, character(1), pattern = "SNP_chr[0-9]+_[0-9]+")
split_1 <- split_1[,-2]

split2$SNP_Split2 <- vapply(split2$FileName_SNP_Grouping1, extract_snp_part, character(1), pattern = "SNP_chr[0-9]+_[0-9]+")
split2 <- split2[,-2]
names(split2)[2] <- "Split_Grouping2"

split3$SNP_Split3 <- vapply(split3$FileName_SNP_Grouping1, extract_snp_part, character(1), pattern = "SNP_chr[0-9]+_[0-9]+")
split3 <- split3[,-2]
names(split3)[2] <- "Split_Grouping3"

# Perform the joins
merged_data <- left_join(split2, split3, by = "Line")
merged_data2 <- left_join(split_1, merged_data, by = "Line")
final_merged <- left_join(merged_data2, split_1removed, by = "Line")

# Handle excess rows
excess_rows <- anti_join(split_1removed, merged_data2, by = "Line")
final_merged <- bind_rows(final_merged, excess_rows)

final_merged <- final_merged %>%
  mutate(across(starts_with("Split_Grouping"), ~replace_na(.x, 0))) %>%
  mutate(Concatenated_Grouping = paste(Split_Grouping1, Split_Grouping2, Split_Grouping3, sep = "--"))

final_merged <- final_merged %>%
  select(everything(), Concatenated_Grouping) %>%
  relocate(Concatenated_Grouping, .after = SNP_Split3)

# Define the function to convert numeric values to 'U' based on SNP checks
convert_to_U <- function(grouping_value, snp_marker_values) {
  if (any(snp_marker_values == 'FAIL') || any(grepl("/", snp_marker_values))) {
    return('U')
  } else {
    return(as.character(grouping_value))
  }
}

# Function to generate the Concatenated_Specific column
generate_concatenated_specific <- function(row) {
  split_values <- c(
    convert_to_U(row['Split_Grouping1'], row[str_split(row['SNP_Split1'], "--")[[1]]]),
    convert_to_U(row['Split_Grouping2'], row[str_split(row['SNP_Split2'], "--")[[1]]]),
    convert_to_U(row['Split_Grouping3'], row[str_split(row['SNP_Split3'], "--")[[1]]])
  )
  return(paste(split_values, collapse = '--'))
}

# Apply the function to create the Concatenated_Specific column
final_merged$Concatenated_Specific <- apply(final_merged, 1, function(row) generate_concatenated_specific(as.list(row)))

# Insert the new column after Concatenated_Grouping
final_merged <- final_merged %>%
  select(1:which(names(final_merged) == 'Concatenated_Grouping'),
         Concatenated_Specific,
         everything())

# Save the updated DataFrame to a new CSV file
write.csv(final_merged, "final_merged_with_specific2.csv", row.names = FALSE)

# Print the resulting dataframe
print(final_merged)
