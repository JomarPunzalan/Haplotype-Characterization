library(dplyr)
rm(list=ls())
# Set the working directory
setwd('C:/Users/JPunzalan/OneDrive - LSU AgCenter/Desktop/Final_Haplotype/repository_snp_1/Split3')

# Get a list of all CSV files in the directory
all_csv_files <- list.files(pattern = "*.csv")

# Exclude files starting with 'snp_haplotype'
csv_files <- all_csv_files[!grepl("^snp_haplotype", all_csv_files, ignore.case = TRUE)]

# Initialize an empty data frame to store the results
all_samples <- data.frame(Line = character(), FileName_SNP_Grouping1 = character(), Split_Grouping1 = integer(), stringsAsFactors = FALSE)

# Loop through each file and extract the SampleID and FileName
for (i in seq_along(csv_files)) {
  file <- csv_files[i]
  
  # Read the current CSV file
  data <- read.csv(file)
  
  # Check if the 'Line' column (SampleID) exists in the data
  if ("Line" %in% colnames(data)) {
    # Extract the SampleID and add the FileName and FileNumber
    file_samples <- data.frame(Line = data$Line, FileName_SNP_Grouping1 = file, Split_Grouping1 = i, stringsAsFactors = FALSE)
    
    # Append to the all_samples data frame
    all_samples <- bind_rows(all_samples, file_samples)
  } else {
    cat("Column 'Line' not found in file:", file, "\n")
  }
}

# Save the results to a CSV file
#write.csv(all_samples, file = "all_samples_with_filenames.csv", row.names = FALSE)

# Print the final data frame for verification
print(all_samples)

snps <- read.csv("snp_haplotype2_525.csv")


#merged_data <- merge(all_samples, snps, by.x = "Line", by.y = "Line", all.x = FALSE)


# Merge the snps data frame with the all_samples data frame using left_join to preserve order
merged_data <- left_join(all_samples, snps, by = c("Line" = "Line"))


path <- "C:/Users/JPunzalan/OneDrive - LSU AgCenter/Desktop/Final_Haplotype/repository_snp_1/Naming"


write.csv(merged_data, file.path(path, "grouping_snp_split3.csv"), row.names=FALSE)


