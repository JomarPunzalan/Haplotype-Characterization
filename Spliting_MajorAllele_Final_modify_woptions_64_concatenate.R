library(openxlsx)
library(dplyr)
library(stringr)

# Load the data
setwd('C:/Users/JPunzalan/OneDrive - LSU AgCenter/Desktop/LSU_Haplotypes_Jen_64')

data <- read.csv("sample_delete.csv")

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
data$Concatenated_Specific <- apply(data, 1, function(row) generate_concatenated_specific(as.list(row)))

# Insert the new column after Concatenated_Grouping
data <- data %>%
  select(1:which(names(data) == 'Concatenated_Grouping'),
         Concatenated_Specific,
         everything())



write.xlsx(data, 'updated_sample_delete.xlsx')
