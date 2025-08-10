
library(tximport)
library(dplyr)
library(purrr)


# Identification of sample names

# Set the directory path (Kallisto files directory)
directory_path <- "/Volumes/Untitled/AsymAD_Project_Rumeysa/Differential Analyses/ROSMAP/DEG/ROSMAP_Kallisto_files_all" 


# Get the file names within the directory
file_names <- list.files(directory_path, full.names = TRUE)

# Extract the sample IDs from the file names
sample_ids <- file_names %>%
  map_chr(~ basename(.x)) %>%
  unique()


# Set the directory paths for the "abundance.tsv" file for each sample

files <- file.path(file_names, "abundance.tsv")

file.exists(files)


# Import transcript-level expression data using tximport

txi <- tximport(files, type = "kallisto", txOut = TRUE) # txOut = TRUE returns transcript-level data

# Transpose the count matrix (samples as rows, transcripts as columns)
count_data_DET <- as.data.frame(txi$counts)

colnames(count_data_DET) <- sample_ids

# Write the result to a csv file
write.csv(count_data_DET, "count_data_DET_ROSMAP.csv", row.names = TRUE)



