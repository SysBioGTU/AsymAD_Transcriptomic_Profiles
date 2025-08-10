
library(dplyr)
library(tidyverse)
library(tximport)
library(rtracklayer)


# Read GTF file
gtf_file <- "Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf"
gtf <- rtracklayer::import(gtf_file)
gtf_df <- as.data.frame(gtf)


# Create a transcript-to-gene mapping file
tx2gene <- gtf_df %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::mutate(TXNAME = transcript_id) %>%
  dplyr::select(TXNAME, GENEID = gene_id)

head(tx2gene)


# Identification of sample names

# Set the directory path
directory_path <- "/Volumes/Untitled/AsymAD_Project_Rumeysa/Differential Analyses/ROSMAP/DEG" # Kallisto files directory

# Get the file names within the directory
file_names <- list.files(directory_path, full.names = TRUE)

# Extract the sample IDs from the file names
sample_ids <- file_names %>%
  map_chr(~ str_extract(basename(.x), "\\d+_\\d+")) %>%
  unique()

sample_ids[635] <- "R24_131017" # Correction for the sample name

sample_ids[636] <- "redo4_140501" # Correction for the sample name again


# Set the directory paths for the "abundance.tsv" file for each sample

files <- file.path(file_names, "abundance.tsv")

file.exists(files)



txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, 
                             txOut = FALSE, # Output: gene-level (default FALSE)
                             ignoreTxVersion = TRUE) # Ignoring ID version


head(txi.kallisto.tsv$counts)


gene_counts <- txi.kallisto.tsv[["counts"]]

gene_counts <- as.data.frame(gene_counts)

colnames(gene_counts) <- sample_ids


# Write the result to a CSV file
write.csv(gene_counts, "count_data_DEG_ROSMAP.csv", row.names = TRUE)







