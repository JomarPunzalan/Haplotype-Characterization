library(openxlsx)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(magick)

rm(list = ls())

setwd('C:/Users/JPunzalan/OneDrive - LSU AgCenter/Desktop/Final_Haplotype/repository_snp_1/Naming')

load("desc.RData")

# Load the data
data <- read.csv("final_merged_with_specific2.csv")

# Define the most_common_value function
most_common_value <- function(x) {
  x <- x[x != 99]
  if(length(x) == 0) return(NA)
  counts <- sort(table(x), decreasing = TRUE)
  top_values <- names(counts[counts == counts[1]])
  if (length(top_values) > 1) {
    return(paste(top_values, collapse = " and "))
  } else {
    return(top_values[1])
  }
}

# Remove rows with '0--0--0' in the 'Concatenated_Grouping' column
data <- data %>% filter(!str_detect(Concatenated_Grouping, "0--0--0"))

# Identify the columns to process (SNP columns)
start_col <- which(names(data) == "Concatenated_Grouping") + 1
snp_columns <- names(data)[start_col:ncol(data)]

# Function to get the most common SNP for each group
get_most_common_snps <- function(df) {
  snp_values <- sapply(df[snp_columns], most_common_value)
  return(as.data.frame(t(snp_values)))
}

# Apply the function to each group and compile the results
compiled_results <- data %>%
  group_by(Concatenated_Grouping) %>%
  do(get_most_common_snps(.)) %>%
  ungroup() %>%
  mutate(Grouping_ID = row_number()) %>%
  left_join(data %>% select(Concatenated_Grouping, SNP_Split1, SNP_Split2, SNP_Split3), by = "Concatenated_Grouping") %>%
  select(SNP_Split1, SNP_Split2, SNP_Split3, Grouping_ID, Concatenated_Grouping, everything())

# Add a column for the count of each Concatenated_Grouping
compiled_results <- compiled_results %>%
  group_by(Concatenated_Grouping) %>%
  mutate(Grouping_Count = n()) %>%
  ungroup() %>%
  select(SNP_Split1, SNP_Split2, SNP_Split3, Grouping_ID, Grouping_Count, Concatenated_Grouping, everything())

# Print the compiled results
print(compiled_results)

# Filter out duplicates across columns using the Grouping_ID
filtered_results <- compiled_results %>%
  distinct(Grouping_ID, .keep_all = TRUE)

# Print the filtered results
print(filtered_results)

# --- New code to reorder columns based on SNP_Split1, SNP_Split2, and SNP_Split3 ---

# Extract the list of non-NA SNP columns under SNP_Split1, SNP_Split2, and SNP_Split3
snps_to_include <- unique(c(
  na.omit(data$SNP_Split1),
  na.omit(data$SNP_Split2),
  na.omit(data$SNP_Split3)
))

# Create a data frame for Preliminary SNP set
preliminary_snp <- data.frame(Number = seq_along(snps_to_include), Preliminary_SNP_Set = snps_to_include)

# Initialize an empty dataframe for the inclusion table
inclusion_result <- data.frame(SNP = character(), Split_Included_In = character())

# Function to check inclusion and append to result
check_inclusion <- function(snp, column) {
  if (snp %in% na.omit(data[[column]])) {
    inclusion_result <<- rbind(inclusion_result, data.frame(SNP = snp, Split_Included_In = column))
  }
}

# Loop through the unique SNPs and check their inclusion in each column
for (snp in snps_to_include) {
  check_inclusion(snp, "SNP_Split1")
  check_inclusion(snp, "SNP_Split2")
  check_inclusion(snp, "SNP_Split3")
}

# Merge while preserving the order of preliminary_snp_set
preliminary_snp_set <- preliminary_snp %>%
  left_join(inclusion_result, by = c("Preliminary_SNP_Set" = "SNP"))

# Display the final result
print(preliminary_snp_set)

# Aggregate data to avoid many-to-many relationship
aggregated_data <- data %>%
  group_by(Concatenated_Grouping) %>%
  summarize(across(all_of(snps_to_include), ~ most_common_value(.x), .names = "agg_{.col}")) %>%
  ungroup()

# Rename aggregated columns to match original names
names(aggregated_data) <- c("Concatenated_Grouping", snps_to_include)

print(aggregated_data)

filtered_results <- filtered_results[, setdiff(names(filtered_results), snps_to_include)]

filtered_results <- left_join(aggregated_data, filtered_results, by = "Concatenated_Grouping")

print(filtered_results)

write.csv(filtered_results, "delete_concatenate.csv", row.names = FALSE)

# Reorder columns to place SNPs after Grouping_ID and Grouping_Count
final_results <- filtered_results %>%
  select(Grouping_ID, Grouping_Count, Concatenated_Grouping, all_of(snps_to_include), everything())

final_results <- final_results %>% select(-Concatenated_Specific)
# Print the final results
print(final_results)

write.csv(final_results, "filtered_commonsnps.csv", row.names = FALSE)

grouping <- filtered_results %>%
  select(Concatenated_Grouping, Grouping_ID)

main_data <- left_join(data, grouping, by = "Concatenated_Grouping")

main_data <- main_data %>%
  select(Grouping_ID, Concatenated_Grouping, Concatenated_Specific, all_of(snps_to_include), everything())

main_data <- main_data[order(main_data$Grouping_ID, decreasing = FALSE),]

snps <- main_data

# Check the structure of the snps data frame
str(snps)

# Filter columns where the frequency of cells with 99 is greater than 10%
threshold <- 0.10
columns_to_keep <- sapply(snps, function(col) mean(col == "FAIL", na.rm = TRUE) <= threshold)
snps <- snps[, columns_to_keep]

# Remove columns where more than 20% of the cells contain heterozygous genotypes (identified by `*/*` or `/` inside cells)
hets_threshold <- 0.20
columns_to_keep <- sapply(snps, function(col) mean(grepl("/|\\*", col), na.rm = TRUE) <= hets_threshold)
main_data_filtered <- snps[, columns_to_keep]

final_results_filtered <- final_results %>%
  select(Grouping_ID, Grouping_Count, Concatenated_Grouping, all_of(snps_to_include), any_of(names(main_data_filtered)))

final_results_filtered <- final_results_filtered %>% select(-Concatenated_Specific)

data <- final_results[, grepl("^SNP_chr", names(final_results))]

snps <- data

for (col in 1:ncol(snps)) {
  for (row in 1:nrow(snps)) {
    # Check if the value is not one of "A", "T", "C", "G", or NA
    if (!(snps[row, col] %in% c("A", "T", "C", "G"))) {
      # Replace with "Hets"
      snps[row, col] <- "1"
    }
  }
}

count_nucleotides_per_column <- function(df_col) {
  # Convert "Hets" values to 1 for processing
  df_col[df_col == "Hets"] <- 1
  
  A_count <- sum(df_col == "A" | is.na(df_col))
  T_count <- sum(df_col == "T" | is.na(df_col))
  C_count <- sum(df_col == "C" | is.na(df_col))
  G_count <- sum(df_col == "G" | is.na(df_col))
  hets_count <- sum(df_col == 1)  # Count "Hets"
  allele_count <- sum(df_col == "A" | df_col == "T" | df_col == "C" | df_col == "G")
  
  sorted_counts <- sort(c(A_count, T_count, C_count, G_count), decreasing = TRUE)
  first_highest <- sorted_counts[1]
  
  # Check if first_highest is not 1 or 99
  if (first_highest != 1 && first_highest != 99) {
    if (A_count == first_highest) {
      df_col <- ifelse(df_col %in% c("A", NA), 2, 0)  # Change highest frequency to 2, others to 0
    } else if (T_count == first_highest) {
      df_col <- ifelse(df_col %in% c("T", NA), 2, 0)
    } else if (C_count == first_highest) {
      df_col <- ifelse(df_col %in% c("C", NA), 2, 0)
    } else if (G_count == first_highest) {
      df_col <- ifelse(df_col %in% c("G", NA), 2, 0)
    }
  }
  
  # Convert back "Hets" values to 1
  df_col[df_col == 1] <- "Hets"
  
  return(df_col)
}

# Apply the function to each column of the data frame
nucleotide_snp_predictions <- as.data.frame(lapply(snps, count_nucleotides_per_column))

# Convert predicted SNP values to character type
nucleotide_snp_predictions <- as.data.frame(lapply(nucleotide_snp_predictions, as.character))

# Replace values in snps with the predicted SNP values
for (col in colnames(snps)) {
  snps[snps[[col]] != "1" & snps[[col]] != "99", col] <- nucleotide_snp_predictions[snps[[col]] != "1" & snps[[col]] != "99", col]
}

# Display the updated dataframe
print(snps)

snps <- data.frame(lapply(snps, function(x) as.numeric(as.character(x))))

snps[snps == "0"] <- -1
snps[snps == "1"] <- 0.5
snps[snps == "2"] <- 1
snps[snps == "99"] <- NA

# Assuming snps is already defined

A <- as.matrix(snps)
B <- as.matrix(t(snps))
dim(A)
dim(B)

# Get the dimensions of the matrices
n <- nrow(A)
m <- ncol(A)
p <- ncol(B)

# Create an empty matrix to store the result
C <- matrix(0, nrow = n, ncol = p)
dim(C)

# Perform matrix multiplication with adjusted rules
for (i in 1:n) {
  for (j in 1:p) {
    for (k in 1:m) {
      if (is.na(A[i, k]) || is.na(B[k, j])) {
        # Treat NA as 0
        C[i, j] <- C[i, j] + 0
      } else if ((A[i, k] == 0.5 && B[k, j] == 1) ||
                 (A[i, k] == 1 && B[k, j] == 0.5) ||
                 (A[i, k] == 0.5 && B[k, j] == -1) ||
                 (A[i, k] == -1 && B[k, j] == 0.5)) {
        C[i, j] <- C[i, j] + 0.5
      } else if (A[i, k] == 1 && B[k, j] == -1) {
        C[i, j] <- C[i, j] + 0
      } else if (A[i, k] == -1 && B[k, j] == 1) {
        C[i, j] <- C[i, j] + 0
      } else if (A[i, k] == 0.5 && B[k, j] == 0.5) {
        C[i, j] <- C[i, j] + 0.25
      } else {
        C[i, j] <- C[i, j] + A[i, k] * B[k, j]
      }
    }
  }
}

# Print the result
print(C)

C0.mat <- as.matrix(C)

C0.mat[row(C0.mat) == col(C0.mat)] <- NA

GRM.CO <- C0.mat / ncol(A)

GRM <- as.matrix(GRM.CO)

column_names <- final_results$Concatenated_Grouping

n <- length(column_names)
matrix_data <- matrix(1:(n * n), nrow = n)

colnames(GRM) <- column_names

GRM <- round(GRM, 2)

write.csv(GRM, "grm.csv")


print(colnames(frequency_merged))

original_haplotype <- original_haplotype[, 1:3] 

original_haplotype <- original_haplotype %>%
  select(Line, everything())

if ("specific_grouping" %in% colnames(frequency_merged) & 
    "general_grouping" %in% colnames(frequency_merged) & 
    "Concatenated_Specific" %in% colnames(frequency_merged)) {
  frequency_merged_topivot <- frequency_merged %>%
    relocate(Concatenated_Specific, specific_grouping, general_grouping, .after = Line) %>% 
    filter(!str_detect(Concatenated_Grouping, "0--0--0"))
} else {
  print("Columns specific_grouping and/or general_grouping and/or Concatenated_Specific do not exist in the dataframe.")
}

pivot <- table(frequency_merged_topivot$Concatenated_Grouping, frequency_merged_topivot$general_grouping)

rownames(pivot) <- NULL
pivot <- as.data.frame.matrix(pivot)
pivot <- cbind(final_results[, 1:3], pivot[, 1:ncol(pivot)])
pivot <- as.data.frame(pivot)
pivot <- pivot %>%
  relocate(Grouping_ID, Concatenated_Grouping, Grouping_Count, everything())

pivot_data <- pivot

# Rename columns in each dataframe before coloring
rename_columns <- function(df) {
  colnames(df) <- str_replace_all(colnames(df), c("Concatenated_Grouping" = "Node_ID", "Concatenated_Specific" = "Node_ID2"))
  return(df)
}

main_data <- rename_columns(main_data)
main_data_filtered <- rename_columns(main_data_filtered)
final_results <- rename_columns(final_results)
final_results_filtered <- rename_columns(final_results_filtered)
frequency_merged <- rename_columns(frequency_merged)
pivot_data <- rename_columns(pivot_data)
GRM <- rename_columns(GRM)



style_A <- createStyle(fontColour = "#9C6500", bgFill = "#FFEB9C")
style_T <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
style_G <- createStyle(fontColour = "#000080", bgFill = "#B8CCE4")
style_C <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
style_with_slash <- createStyle(fontColour = "#FFFFFF", bgFill = "#3E424B")
style_FAIL <- createStyle(fontColour = "#FFFFFF", bgFill = "#3E424B")

# New custom style for snps_to_include
style_snp_include <- createStyle(fontColour = "#FF2615")

wb <- createWorkbook()

# Add the Description sheet first
addWorksheet(wb, "Description")
writeData(wb, "Description", description)

# Add the Workflow sheet second
addWorksheet(wb, "Workflow")

# Add other sheets
sheets_data <- list(
  "Preliminary SNP Set" = preliminary_snp_set,
  "Main Data Unfiltered SNPs" = main_data,
  "Main Data Filtered SNPs" = main_data_filtered,
  "Rep_HaplotypeUnfiltered" = final_results,
  "Rep_HaplotypeFiltered" = final_results_filtered,
  "SNP Frequency" = frequency_merged,
  "Cross-Tabulation" = pivot_data,
  "GRM_AlleleSimilarity" = GRM
)

for (sheet in names(sheets_data)) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, sheets_data[[sheet]])
}

# Load the image and add it to the workbook
image_path <- "C:/Users/JPunzalan/OneDrive - LSU AgCenter/Desktop/Final_Haplotype/repository_snp_1/Naming/Haplotype Workflow.jpg"
image <- image_read(image_path)
image_write(image, path = "Haplotype_Workflow.png", format = "png")

insertImage(wb, "Workflow", file = "Haplotype_Workflow.png",  width = 14.51, height = 20)

apply_snp_coloring_bulk <- function(wb, sheet, data, styles) {
  snp_columns <- grep("^SNP_chr", colnames(data))
  if (length(snp_columns) > 0) {
    for (rule in names(styles)) {
      conditionalFormatting(
        wb, sheet = sheet, cols = snp_columns, rows = 2:(nrow(data) + 1),
        type = "contains", rule = rule, style = styles[[rule]]
      )
    }
  }
}

# New function to apply custom color
apply_custom_coloring <- function(wb, sheet, data, snps, style) {
  snp_columns <- which(colnames(data) %in% snps)
  if (length(snp_columns) > 0) {
    addStyle(wb, sheet, style = style, rows = 1, cols = snp_columns, gridExpand = TRUE, stack = TRUE)
    addStyle(wb, sheet, style = style, rows = 2:(nrow(data) + 1), cols = snp_columns, gridExpand = TRUE, stack = TRUE)
  }
}

# Function to apply custom font color to entire columns in the "Preliminary SNP Set" sheet
apply_custom_font_color <- function(wb, sheet, style) {
  col_index <- which(colnames(preliminary_snp_set) == "Preliminary_SNP_Set")
  addStyle(wb, sheet, style = style, rows = 1:(nrow(preliminary_snp_set) + 1), cols = col_index, gridExpand = TRUE, stack = TRUE)
}

snp_styles <- list("A" = style_A, "T" = style_T, "G" = style_G, "C" = style_C, "/" = style_with_slash, "FAIL" = style_FAIL)

apply_snp_coloring_bulk(wb, "Main Data Unfiltered SNPs", main_data, snp_styles)
apply_snp_coloring_bulk(wb, "Main Data Filtered SNPs", main_data_filtered, snp_styles)
apply_snp_coloring_bulk(wb, "Rep_HaplotypeFiltered", final_results_filtered, snp_styles)
apply_snp_coloring_bulk(wb, "Rep_HaplotypeUnfiltered", final_results, snp_styles)
apply_snp_coloring_bulk(wb, "SNP Frequency", frequency_merged, snp_styles)
apply_snp_coloring_bulk(wb, "Cross-Tabulation", pivot_data, snp_styles)
apply_snp_coloring_bulk(wb, "GRM_AlleleSimilarity", GRM, snp_styles)

# Apply custom coloring for snps_to_include
apply_custom_font_color(wb, "Preliminary SNP Set",style_snp_include)
apply_custom_coloring(wb, "Main Data Unfiltered SNPs", main_data, snps_to_include, style_snp_include)
apply_custom_coloring(wb, "Main Data Filtered SNPs", main_data_filtered, snps_to_include, style_snp_include)
apply_custom_coloring(wb, "Rep_HaplotypeFiltered", final_results_filtered, snps_to_include, style_snp_include)
apply_custom_coloring(wb, "Rep_HaplotypeUnfiltered", final_results, snps_to_include, style_snp_include)
apply_custom_coloring(wb, "SNP Frequency", frequency_merged, snps_to_include, style_snp_include)
apply_custom_coloring(wb, "Cross-Tabulation", pivot_data, snps_to_include, style_snp_include)
apply_custom_coloring(wb, "GRM_AlleleSimilarity", GRM, snps_to_include, style_snp_include)

numeric_columns <- sapply(pivot_data, is.numeric)

for (col in which(numeric_columns)[which(numeric_columns) >= 3]) {
  col_data <- pivot_data[, col]
  col_max <- max(col_data, na.rm = TRUE)
  col_min <- min(col_data, na.rm = TRUE)
  
  conditionalFormatting(
    wb, sheet = "Cross-Tabulation", cols = col, rows = 2:(nrow(pivot_data) + 1),
    type = "colorScale", style = c("#F8696B", "#FEEB84", "#63BE7B"),
    values = c(col_min, (col_max + col_min) / 2, col_max)
  )
}

# Save the workbook
saveWorkbook(wb, "output.xlsx", overwrite = TRUE)

