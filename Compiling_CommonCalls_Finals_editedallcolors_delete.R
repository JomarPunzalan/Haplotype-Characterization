install.packages("openxlsx")

library(openxlsx)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
rm(list = ls())

setwd('C:/Users/JPunzalan/OneDrive - LSU AgCenter/Desktop/LSU_Haplotypes_Jen_64/repository_snp_6/Naming')

# Load the data
data <- read.csv("final_merged_with_specific.csv")

#data <- subset(data, Split_Grouping1 != "0" )

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

# Print the final results
print(final_results)

write.csv(final_results, "filtered_commonsnps.csv", row.names = FALSE)


grouping <- filtered_results %>%
  select(Concatenated_Grouping, Grouping_ID,)


main_data <- left_join(data, grouping, by = "Concatenated_Grouping")

main_data <- main_data %>%
  select(Grouping_ID, Concatenated_Grouping, Concatenated_Specific, all_of(snps_to_include), everything())


main_data <- main_data[order(main_data$Grouping_ID, decreasing = FALSE),]



#list_of_datasets <- list("Main Data" = main_data, "Common SNPs" = final_results)
#write.xlsx(list_of_datasets, file = "compiled_results.xlsx")




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

C0.mat<-as.matrix(C)

C0.mat[row(C0.mat) == col(C0.mat)] <- NA

#C0.max <- max(C0.mat[upper.tri(C0.mat, diag=FALSE)], na.rm=TRUE) # take max of upper triangle
#C0.max

GRM.CO<-C0.mat/ncol(A)

GRM <- as.matrix(GRM.CO)


# Retain only the upper triangle
#GRM[!upper.tri(GRM)] <- NA  # or NA


column_names <- final_results$Concatenated_Grouping


# Create a sample matrix
n <- length(column_names)
matrix_data <- matrix(1:(n*n), nrow = n)

# Set column names using the extracted values
colnames(GRM) <- column_names

GRM <- round(GRM,2)

write.csv(GRM, "grm.csv")


frequency_merged <- read.csv("final_merged.csv")

original_haplotype <- read.csv("snp_haplotype2_525.csv")
original_haplotype <- original_haplotype[,1:3] 

original_haplotype <- original_haplotype %>%
  select(Line, everything())

frequency_merged <- left_join(frequency_merged, original_haplotype, by= "Line")%>%
  relocate(specific_grouping,	general_grouping, .after = Line) %>% 
  filter(!str_detect(Concatenated_Grouping, "0--0--0"))


pivot <- table(frequency_merged$Concatenated_Grouping, frequency_merged$general_grouping)

rownames(pivot)<-NULL
pivot <- as.data.frame.matrix(pivot)
pivot<-cbind(final_results[,1:3], pivot[,1:ncol(pivot)])
pivot <- as.data.frame(pivot)
pivot <- pivot %>%
  relocate(Grouping_ID,  Concatenated_Grouping, Grouping_Count, everything())

pivot_data <- pivot

#jeje
# Define cell styles for SNP values and "FAIL"
style_A <- createStyle(fontColour = "#9C6500", bgFill = "#FFEB9C") # Red background, black font
style_T <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE") # Green background, black font
style_G <- createStyle(fontColour = "#000080", bgFill = "#B8CCE4") # Blue background, black font
style_C <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE") # Orange background, black font
style_with_slash <- createStyle(fontColour = "#000000", bgFill = "#CC6666") # Yellow background, black font
style_FAIL <- createStyle(fontColour = "#000000", bgFill = "#996666") # Red background, white font

# Integrate SNP coloring into the output file
wb <- createWorkbook()

# Add worksheets for the datasets
addWorksheet(wb, "Main Data")
addWorksheet(wb, "Representative Haplotype")
addWorksheet(wb, "Frequency")
addWorksheet(wb, "Pivot")
addWorksheet(wb, "GRM_AlleleSimilarity")

# Write data to the worksheets
writeData(wb, "Main Data", main_data)
writeData(wb, "Representative Haplotype", final_results)
writeData(wb, "Frequency", frequency_merged)
writeData(wb, "Pivot", pivot_data)
writeData(wb, "GRM_AlleleSimilarity", GRM)

# Apply conditional formatting to "Main Data"
snp_columns_main_data <- grep("^SNP_chr", colnames(main_data), value = TRUE)
for (col_name in snp_columns_main_data) {
  col_index <- which(colnames(main_data) == col_name)
  conditionalFormatting(wb, sheet = "Main Data", cols = col_index, rows = 2:(nrow(main_data) + 1), type = "contains", rule = "A", style = style_A)
  conditionalFormatting(wb, sheet = "Main Data", cols = col_index, rows = 2:(nrow(main_data) + 1), type = "contains", rule = "T", style = style_T)
  conditionalFormatting(wb, sheet = "Main Data", cols = col_index, rows = 2:(nrow(main_data) + 1), type = "contains", rule = "G", style = style_G)
  conditionalFormatting(wb, sheet = "Main Data", cols = col_index, rows = 2:(nrow(main_data) + 1), type = "contains", rule = "C", style = style_C)
  conditionalFormatting(wb, sheet = "Main Data", cols = col_index, rows = 2:(nrow(main_data) + 1), type = "contains", rule = "/", style = style_with_slash)
  conditionalFormatting(wb, sheet = "Main Data", cols = col_index, rows = 2:(nrow(main_data) + 1), type = "contains", rule = "FAIL", style = style_FAIL)
}

# Apply conditional formatting to "Common SNPs"
snp_columns_final_results <- grep("^SNP_chr", colnames(final_results), value = TRUE)
for (col_name in snp_columns_final_results) {
  col_index <- which(colnames(final_results) == col_name)
  conditionalFormatting(wb, sheet = "Representative Haplotype", cols = col_index, rows = 2:(nrow(final_results) + 1), type = "contains", rule = "A", style = style_A)
  conditionalFormatting(wb, sheet = "Representative Haplotype", cols = col_index, rows = 2:(nrow(final_results) + 1), type = "contains", rule = "T", style = style_T)
  conditionalFormatting(wb, sheet = "Representative Haplotype", cols = col_index, rows = 2:(nrow(final_results) + 1), type = "contains", rule = "G", style = style_G)
  conditionalFormatting(wb, sheet = "Representative Haplotype", cols = col_index, rows = 2:(nrow(final_results) + 1), type = "contains", rule = "C", style = style_C)
  conditionalFormatting(wb, sheet = "Representative Haplotype", cols = col_index, rows = 2:(nrow(final_results) + 1), type = "contains", rule = "/", style = style_with_slash)
  conditionalFormatting(wb, sheet = "Representative Haplotype", cols = col_index, rows = 2:(nrow(final_results) + 1), type = "contains", rule = "FAIL", style = style_FAIL)
}

# Apply conditional formatting to "Frequency"
#snp_columns_frequency <- grep("^SNP_chr", colnames(frequency_merged), value = TRUE)
#for (col_name in snp_columns_frequency) {
#  col_index <- which(colnames(frequency_merged) == col_name)
#  conditionalFormatting(wb, sheet = "Frequency", cols = col_index, rows = 2:(nrow(frequency_merged) + 1), type = "contains", rule = "A", style = style_A)
#  conditionalFormatting(wb, sheet = "Frequency", cols = col_index, rows = 2:(nrow(frequency_merged) + 1), type = "contains", rule = "T", style = style_T)
#  conditionalFormatting(wb, sheet = "Frequency", cols = col_index, rows = 2:(nrow(frequency_merged) + 1), type = "contains", rule = "G", style = style_G)
#  conditionalFormatting(wb, sheet = "Frequency", cols = col_index, rows = 2:(nrow(frequency_merged) + 1), type = "contains", rule = "C", style = style_C)
#  conditionalFormatting(wb, sheet = "Frequency", cols = col_index, rows = 2:(nrow(frequency_merged) + 1), type = "contains", rule = "/", style = style_with_slash)
#  conditionalFormatting(wb, sheet = "Frequency", cols = col_index, rows = 2:(nrow(frequency_merged) + 1), type = "contains", rule = "FAIL", style = style_FAIL)
#}


start_column <- 3

for (col in start_column:ncol(pivot_data)) {
  col_index <- col  # Use the actual column index for the conditional formatting
  col_data <- pivot_data[, col]
  
  # Ensure the column data is numeric
  if (is.numeric(col_data)) {
    col_max <- max(col_data, na.rm = TRUE)
    col_min <- min(col_data, na.rm = TRUE)
    
    # Apply a color scale from red (lowest) to green (highest) within the column
    conditionalFormatting(wb, sheet = "Pivot", cols = col_index, rows = 2:(nrow(pivot_data) + 1), type = "colorScale",
                          style = c("#F8696B", "#FEEB84", "#63BE7B"), values = c(col_min, (col_max + col_min) / 2, col_max))
  } else {
    # Handle non-numeric columns, if necessary
    message(paste("Skipping non-numeric column:", colnames(pivot_data)[col]))
  }
}

# Save the workbook after applying the formatting
saveWorkbook(wb, "your_output_workbook_specific_grouping.xlsx", overwrite = TRUE)





















library(openxlsx)

# Create styles
styles <- list(
  "A" = createStyle(fontColour = "#9C6500", bgFill = "#FFEB9C"),
  "T" = createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE"),
  "G" = createStyle(fontColour = "#000080", bgFill = "#B8CCE4"),
  "C" = createStyle(fontColour = "#006100", bgFill = "#C6EFCE"),
  "/" = createStyle(fontColour = "#000000", bgFill = "#CC6666"),
  "FAIL" = createStyle(fontColour = "#000000", bgFill = "#996666")
)

# Create workbook and add sheets
wb <- createWorkbook()
sheets <- c("Main Data", "Representative Haplotype", "Frequency", "Pivot", "GRM_AlleleSimilarity")
for (sheet in sheets) {
  addWorksheet(wb, sheet)
}

# Write data to sheets
writeData(wb, "Main Data", main_data)
writeData(wb, "Representative Haplotype", final_results)
writeData(wb, "Frequency", frequency_merged)
writeData(wb, "Pivot", pivot_data)
writeData(wb, "GRM_AlleleSimilarity", GRM)

# Function to apply conditional formatting efficiently
apply_conditional_formatting <- function(wb, sheet_name, data) {
  snp_columns <- grep("^SNP_chr", colnames(data), value = TRUE)
  for (col_name in snp_columns) {
    col_index <- which(colnames(data) == col_name)
    rules <- names(styles)
    for (rule in rules) {
      conditionalFormatting(
        wb, sheet = sheet_name, cols = col_index, rows = 2:(nrow(data) + 1), 
        type = "contains", rule = rule, style = styles[[rule]]
      )
    }
  }
}

# Apply conditional formatting to the relevant sheets
apply_conditional_formatting(wb, "Main Data", main_data)
apply_conditional_formatting(wb, "Representative Haplotype", final_results)

# Save the workbook
saveWorkbook(wb, "SNP_Analysis.xlsx", overwrite = TRUE)





























saveWorkbook(wb, "compiled_results_colored_delete.xlsx", overwrite = TRUE)


saveWorkbook(wb, "compiled_results_colored_delete.xlsx", overwrite = TRUE)













install.packages("openxlsx")
library(openxlsx)

# Create a sample data frame
data <- data.frame(Column1 = 1:50)

# Create a new Excel workbook and add the data
wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", data)

# Define the conditional formatting rule
conditionalFormatting(wb, "Sheet1", cols = 1, rows = 1:50, 
                      style = createStyle(fontColour = "#000000"),
                      rule = "1:50",
                      type = "colorScale",
                      colors = c("red", "yellow", "green"))

# Save the workbook
saveWorkbook(wb, "conditional_formatting_example.xlsx", overwrite = TRUE)







install.packages("openxlsx")
library(openxlsx)

# Create a sample data frame
data <- data.frame(Column1 = 1:50)

# Create a new Excel workbook and add the data
wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", data)

# Define the conditional formatting rule
conditionalFormatting(wb, "Sheet1", cols = 1, rows = 1:50, 
                      type = "colorScale", 
                      colors = c("FF0000", "FFFF00", "00FF00"))

# Save the workbook
saveWorkbook(wb, "conditional_formatting_example.xlsx", overwrite = TRUE)

