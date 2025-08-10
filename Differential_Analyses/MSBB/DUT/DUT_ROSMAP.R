# DUT

library(DRIMSeq)

library(dplyr)

library(readr)


# Load the metadata

met <- read.csv("met.csv", sep="")

col_Data_out_2 <- met[-c(2, 6, 7, 8, 9)]  # sample_id, msex, age_death, pmi, group, and Sample columns

# Convert categorical variables into factors for model design
col_Data_out_2$msex <- as.factor(col_Data_out_2$msex)
col_Data_out_2$group <- as.factor(col_Data_out_2$group)



# Load the transcript count data

combined_count_data <- read.csv("/Volumes/Untitled/AsymAD_Project_Rumeysa/Differential Analyses/ROSMAP/DET/combined_count_data.csv", row.names = 1, , check.names = FALSE)

combined_count_data <- as.data.frame(combined_count_data)


# making sure the row names in col_Data_2 matches to column names in counts_data
all(colnames(combined_count_data) %in% col_Data_out_2$sample_id)

# are they in the same order?
all(colnames(combined_count_data) == col_Data_out_2$sample_id)


#--------------------------------------------------------------------------------------------------------------------


# Read the reference data (transcripts_to_genes.txt) for "hg19" that includes gene IDs, trancript IDs and gene names ()
# This file was downloada from https://github.com/pachterlab/kallisto-transcriptome-indices/releases

txdf <- read.table("transcripts_to_genes.txt", sep = "\t")

colnames(txdf) <- c("transcript_id", "gene_id", "gene_name")


txdf.sub1 = txdf[match(rownames(combined_count_data),txdf$transcript_id),]
counts_with_IDs = data.frame(gene_id = txdf.sub1$gene_id, feature_id = txdf.sub1$transcript_id, combined_count_data_out)

library(tidyr)

# Remove rows with NA values in the feature_id column
counts_with_IDs <- drop_na(counts_with_IDs, feature_id)



# To create a dmDSdata object
d = dmDSdata(counts = counts_with_IDs, samples = col_Data_out_2)

# head(counts(d), 3)
# head(samples(d), 3)

# This part is filtering out lowly expressed genes like CPM filtration in DEG analysis
d <- dmFilter(d, min_samps_gene_expr = 485, min_samps_feature_expr = 92, min_gene_expr = 10, min_feature_expr = 10)


# Create the design matrix
design_full <- model.matrix(~ -1 + group + msex + age_death + pmi, data = d@samples)

# design_full

# Calculate precision and fit
d <- dmPrecision(d, design = design_full)

# genewise_precision(d)

d <- dmFit(d, design = design_full)


# Control vs AD ------------------------------------------------------------------------------------------------------------------

contrast1 <- c(-1, 0, 1, 0, 0, 0)

d_c1 <- dmTest(d, contrast = contrast1, verbose = TRUE)

# d_c1@design_fit_null



# res_c1.2 for transcripts


res_c1.2 <- d_c1@results_feature
res_c1.2 <- res_c1.2[order(res_c1.2$pvalue, decreasing = FALSE), ]


p_val_2 = res_c1.2$pvalue 

p_g1.1 = which(p_val_2 < 0.05)

p_g1.2 = which(p_val_2 < 0.01)

p_g1.3 = which(p_val_2 < 0.001)


p_adj_2 = res_c1.2$adj_pvalue 

a_g1.1 <- which(p_adj_2  < 0.05)

a_g1.2 <- which(p_adj_2  < 0.01)

a_g1.3 <- which(p_adj_2  < 0.001)



rownames(res_c1.2) <- res_c1.2$feature_id

transcripts_p_g1.1 <- rownames(res_c1.2[p_g1.1,])
transcripts_p_g1.2 <- rownames(res_c1.2[p_g1.2,])
transcripts_p_g1.3 <- rownames(res_c1.2[p_g1.3,])

transcripts_a_g1.1 <- rownames(res_c1.2[a_g1.1,])
transcripts_a_g1.2 <- rownames(res_c1.2[a_g1.2,])
transcripts_a_g1.3 <- rownames(res_c1.2[a_g1.3,])


write.table(transcripts_p_g1.1, "transcripts_AD_Control_pval1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(transcripts_p_g1.2, "transcripts_AD_Control_pval2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(transcripts_p_g1.3, "transcripts_AD_Control_pval3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(transcripts_a_g1.1, "transcripts_AD_Control_padj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(transcripts_a_g1.2, "transcripts_AD_Control_padj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(transcripts_a_g1.3, "transcripts_AD_Control_padj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)




# Control vs AsymAD

contrast2 <- c(0, -1, 1, 0, 0, 0)

d_c2 <- dmTest(d, contrast = contrast2, one_way = FALSE, verbose = TRUE)



# res_c2.2 for transcripts


res_c2.2 <- d_c2@results_feature
res_c2.2 <- res_c2.2[order(res_c2.2$pvalue, decreasing = FALSE), ]


p_val_4 = res_c2.2$pvalue 

p_g4.1 = which(p_val_4 < 0.05)

p_g4.2 = which(p_val_4 < 0.01)

p_g4.3 = which(p_val_4 < 0.001)


p_adj_4 = res_c2.2$adj_pvalue 

a_g4.1 <- which(p_adj_4  < 0.05)

a_g4.2 <- which(p_adj_4  < 0.01)

a_g4.3 <- which(p_adj_4  < 0.001)


rownames(res_c2.2) <- res_c2.2$feature_id

transcripts_p_g4.1 <- rownames(res_c2.2[p_g4.1,])
transcripts_p_g4.2 <- rownames(res_c2.2[p_g4.2,])
transcripts_p_g4.3 <- rownames(res_c2.2[p_g4.3,])

transcripts_a_g4.1 <- rownames(res_c2.2[a_g4.1,])
transcripts_a_g4.2 <- rownames(res_c2.2[a_g4.2,])
transcripts_a_g4.3 <- rownames(res_c2.2[a_g4.3,])


write.table(transcripts_p_g4.1, "transcripts_AsymAD_Control_pval1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(transcripts_p_g4.2, "transcripts_AsymAD_Control_pval2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(transcripts_p_g4.3, "transcripts_AsymAD_Control_pval3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(transcripts_a_g4.1, "transcripts_AsymAD_Control_padj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(transcripts_a_g4.2, "transcripts_AsymAD_Control_padj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(transcripts_a_g4.3, "transcripts_AsymAD_Control_padj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



# AD vs AsymAD

contrast3 <- c(-1, 1, 0, 0, 0, 0)

d_c3 <- dmTest(d, contrast = contrast3, one_way = FALSE, verbose = TRUE)



# res_3.2 for transcripts


res_c3.2 <- d_c3@results_feature
res_c3.2 <- res_c3.2[order(res_c3.2$pvalue, decreasing = FALSE), ]


p_val_6 = res_c3.2$pvalue 

p_g7.1= which(p_val_6 < 0.05)

p_g7.2= which(p_val_6 < 0.01)

p_g7.3 = which(p_val_6 < 0.001)


p_adj_6 = res_c3.2$adj_pvalue 

a_g7.1 <- which(p_adj_6  < 0.05)

a_g7.2 <- which(p_adj_6  < 0.01)

a_g7.3 <- which(p_adj_6  < 0.001)


rownames(res_c3.2) <- res_c3.2$feature_id

transcripts_p_g7.1 <- rownames(res_c3.2[p_g7.1,])
transcripts_p_g7.2 <- rownames(res_c3.2[p_g7.2,])
transcripts_p_g7.3 <- rownames(res_c3.2[p_g7.3,])

transcripts_a_g7.1 <- rownames(res_c3.2[a_g7.1,])
transcripts_a_g7.2 <- rownames(res_c3.2[a_g7.2,])
transcripts_a_g7.3 <- rownames(res_c3.2[a_g7.3,])


write.table(transcripts_p_g7.1, "transcripts_AD_AsymAD_Control_pval1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(transcripts_p_g7.2, "transcripts_AD_AsymAD_Control_pval2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(transcripts_p_g7.3, "transcripts_AD_AsymAD_Control_pval3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(transcripts_a_g7.1, "transcripts_AD_AsymAD_Control_padj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(transcripts_a_g7.2, "transcripts_AD_AsymAD_Control_padj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(transcripts_a_g7.3, "transcripts_AD_AsymAD_Control_padj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)




# Genes

data <- read.csv("mart_export.txt")

transcripts <- rownames(res_c1.2)

# Filter the transcripts from the data based on the ones in the transcripts list
filtered_data <- data[data$Transcript.stable.ID.version %in% transcripts, ]

# Select the relevant columns
filtered_data <- filtered_data[, c("Transcript.stable.ID.version", "Gene.stable.ID")]


library(dplyr)

genes_1 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_p_g1.1) %>%
  pull(Gene.stable.ID) %>% unique()


genes_2 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_p_g1.2) %>%
  pull(Gene.stable.ID) %>% unique()

genes_3 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_p_g1.3) %>%
  pull(Gene.stable.ID) %>% unique()


genes_4 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_a_g1.1) %>%
  pull(Gene.stable.ID) %>% unique()

genes_5 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_a_g1.2) %>%
  pull(Gene.stable.ID) %>% unique()

genes_6 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_a_g1.3) %>%
  pull(Gene.stable.ID) %>% unique()


write.table(genes_1, "genes_AD_Control_pval1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(genes_2, "genes_AD_Control_pval2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(genes_3, "genes_AD_Control_pval3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_4, "genes_AD_Control_padj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(genes_5, "genes_AD_Control_padj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(genes_6, "genes_AD_Control_padj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



#-------------------------------------------------------------------------------------------------------------------


genes_7 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in%  transcripts_p_g4.1) %>%
  pull(Gene.stable.ID) %>% unique()

genes_8 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_p_g4.2) %>%
  pull(Gene.stable.ID) %>% unique()

genes_9 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_p_g4.3) %>%
  pull(Gene.stable.ID) %>% unique()


genes_10 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_a_g4.1) %>%
  pull(Gene.stable.ID) %>% unique()

genes_11 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_a_g4.2) %>%
  pull(Gene.stable.ID) %>% unique()

genes_12 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_a_g4.3) %>%
  pull(Gene.stable.ID) %>% unique()



write.table(genes_7, "genes_AsymAD_Control_pval1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(genes_8, "genes_AsymAD_Control_pval2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(genes_9, "genes_AsymAD_Control_pval3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_10, "genes_AsymAD_Control_padj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(genes_11, "genes_AsymAD_Control_padj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(genes_12, "genes_AsymAD_Control_padj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


#--------------------------------------------------------------------------------------------------------------------

genes_13 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_p_g7.1) %>%
  pull(Gene.stable.ID) %>% unique()


genes_14 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_p_g7.2) %>%
  pull(Gene.stable.ID) %>% unique()

genes_15 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_p_g7.3) %>%
  pull(Gene.stable.ID) %>% unique()


genes_16 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_a_g7.1) %>%
  pull(Gene.stable.ID) %>% unique()

genes_17 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_a_g7.2) %>%
  pull(Gene.stable.ID) %>% unique()

genes_18 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% transcripts_a_g7.3) %>%
  pull(Gene.stable.ID) %>% unique()



write.table(genes_13, "genes_AD_AsymAD_pval1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(genes_14, "genes_AD_AsymAD_pval2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(genes_15, "genes_AD_AsymAD_pval3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_16, "genes_AD_AsymAD_padj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(genes_17, "genes_AD_AsymAD_padj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(genes_18, "genes_AD_AsymAD_padj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



