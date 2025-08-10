
# DET Analysis

# Load count data
count_data_DET <- read.csv("count_data_DET.csv", row.names = 1, check.names = FALSE)

load("/Volumes/Untitled/AsymAD_Project_Rumeysa/Differential Analyses/ROSMAP/DEG/rosmapmetadata.RData") # load the metadata


# At first, the three groups as Consensus Control, Consensus AD and AsymAD groups were determined.


# Consensus Control

cons_cont <- metadata[(metadata$cogdx == 1 | metadata$cogdx == 2 | metadata$cogdx == 3 ) & (metadata$ceradsc == 3 | metadata$ceradsc == 4), ]


# Consensus AD 

cons_AD <- metadata[(metadata$cogdx == 4 | metadata$cogdx == 5) & (metadata$ceradsc == 1 | metadata$ceradsc == 2), ]


# Consensus AsymAD

cons_AsymAD <- metadata[metadata$cogdx == 1 & (metadata$ceradsc == 1 | metadata$ceradsc == 2), ]


# The groups were chosen from countdata.

a <- cons_cont$rnaseq_id  #for the control IDs

b <- cons_AD_1$rnaseq_id  #for the AD IDs

c <- cons_AsymAD_1$rnaseq_id  #for the AsymAD IDs


# Define the order of columns based on 'a', 'b', and 'c'
order_of_columns <- c(a, b, c)


# making sure the row names in colData matches to column names in counts_data
all(colnames(count_data_DET) %in% metadata$rnaseq_id)


# Find the differing column names
difference <- setdiff(colnames(count_data_DET), metadata$rnaseq_id)

# Print the differing column names
print(difference)    # "764_130520" 

# Select columns in the desired order
combined_count_data <- count_data_DET[, order_of_columns]


combined_count_out = combined_count_data[, c(-12, -203, -336)]


# Write the combined_count_out to a csv file
write.csv(combined_count_out, "combined_count_data.csv", row.names = TRUE)


# CPM Normalization----------------------------------------------------------------------------------------------

library(edgeR)

# Calculate CPM values and apply cpm threshold
myCPM <- cpm(combined_count_data,log=F)  
thresh <- myCPM > 0.1 # library size of ROSMAP ~50M. 5 count/ 50 = 0.1 cpm threshold
keep <- rowSums(thresh) >= 92 # Number of samples of the condition with the least number of samples  

summary(keep) 

cpm_filtered_count <- combined_count_data[keep,] 


# Create a metadata table with sample information
col_Data <- data.frame(sample_id = colnames(combined_count_data),
                       group = c(rep("Control", length(a)),
                                 rep("AD", length(b)),
                                 rep("AsymAD", length(c))))

# Set the row names to match the sample IDs
row.names(col_Data) <- col_Data$sample_id

# making sure the row names in colData matches to column names in counts_data
all(colnames(combined_count_data) %in% rownames(col_Data))

# are they in the same order?
all(colnames(combined_count_data) == rownames(col_Data))


# Outlier seperation--------------------------------------------------------------------------------------------

# 380_120503, 367_120502 and 500_120515 samples are the outliers.

# Outlier sample names
outlier_samples <- c("380_120503", "367_120502", "500_120515")

# Get their column indices from the count matrix
outlier_indices <- which(names(cpm_filtered_count) %in% outlier_samples)

# Remove the outliers
cpm_filtered_count_out <- cpm_filtered_count[, -outlier_indices]
col_Data_out <- col_Data[-outlier_indices, ]



# Construct a DESeqDataSet object ----------

library(DESeq2)

cpm_filtered_count_out <- round(cpm_filtered_count_out)  # DESeqDataSetFromMatrix() function accepts only integer values

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cpm_filtered_count_out,
                              colData = col_Data_out,
                              design = ~ group)



# Normalization step---------------------------------------------------------------------------------------------

dds <- estimateSizeFactors(dds)

sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE) 


#----------------COVARIATE ADJUSTMENT-----------------------


data_1 = normalized_counts

data_1 = as.data.frame(data_1)

data_t1 <- as.data.frame(t(data_1))


# Rearrange metadata
sample_names <- colnames(data_1)  # Extract sample names
samples <- as.data.frame(sample_names)

# Perform the merge using the "rnaseq_id" column as the key
met <- merge(samples, metadata, by.x = "sample_names", by.y = "rnaseq_id", sort = FALSE)

met['group'] = col_Data_out['group']


sex.num = met$msex

condition.num = met$group

# Replacing "90+" with "90"
age= as.data.frame(met$age_death)
age[age=="90+"] = "90" 
age <- as.numeric(gsub(",", ".", age$`met$age_death`))

pmi= as.numeric(met$pmi)


# Fit the model after obtaining all variables

fit = lm(as.matrix(data_t1)~ sex.num + age + condition.num + pmi)

coef = fit[["coefficients"]]

#----------------------------------------------------

fitteddata = matrix(0, 104978, 485)  # The size of the count matrix
cins = sex.num
cins[cins=="female"]=0
cins[cins=="male"]=1
cins= as.numeric(cins)

state = condition.num
state[state=="AD"]=0
state[state=="Control"]=1
state[state=="AsymAD"]=2
state = as.numeric(state)


genenames <- character(104978)

for (j in 1:485){
  
  for (i in 1:104978) {
    genename = colnames(coef)[i]
    con = state[j]
    gender = cins[j]
    deathage = age[j]
    pmitime = pmi[j]
    
    count = data_1[i,j]-gender*coef[2,i]-deathage*coef[3,i]-pmitime*coef[6,i] 
    genenames[i] <- genename
    fitteddata[i,j] = count
  }
}


colnames(fitteddata) = colnames(cpm_filtered_count_out)
rownames(fitteddata) = genenames
fitteddata <- as.data.frame(fitteddata)



# Limma 

library(limma)

# AD vs Control

# Convert the group information to a numeric matrix
design <- model.matrix(~ 0 + col_Data_out$group)

# Rename the columns with appropriate names (e.g., Control, AsymAD, AD)
colnames(design) <- c("AD", "AsymAD", "Control")


# Create the contrasts for AD vs Control
contrasts_AD_vs_Control <- makeContrasts(
  AD_vs_Control = AD - Control,
  levels = design
)

# Fit the linear model
fit <- lmFit(fitteddata, design)

# Apply the contrasts for AD vs Control
contrast_fit_AD_vs_Control <- contrasts.fit(fit, contrasts_AD_vs_Control)

# Run eBayes on the contrast fit
contrast_ebayes_AD_vs_Control <- eBayes(contrast_fit_AD_vs_Control)

# Retrieve the results for AD vs Control
results_table_AD_vs_Control <- topTable(contrast_ebayes_AD_vs_Control, coef = "AD_vs_Control",  number = nrow(fitteddata))


# Transcript results

# Filter genes based on adjusted p-value threshold (e.g., 0.05)
significant_transcripts_p_val1 <- subset(results_table_AD_vs_Control, P.Value < 0.05)

# Extract gene names from the results table
significant_transcripts_p_val1 <- rownames(significant_transcripts_p_val1)


significant_transcripts_p_val2 <- subset(results_table_AD_vs_Control, P.Value < 0.01)

significant_transcripts_p_val2 <- rownames(significant_transcripts_p_val2)


significant_transcripts_p_val3 <- subset(results_table_AD_vs_Control, P.Value < 0.001)

significant_transcripts_p_val3 <- rownames(significant_transcripts_p_val3)



significant_transcripts_p_adj1 <- subset(results_table_AD_vs_Control, adj.P.Val < 0.05)

significant_transcripts_p_adj1 <- rownames(significant_transcripts_p_adj1)


significant_transcripts_p_adj2 <- subset(results_table_AD_vs_Control, adj.P.Val < 0.01)

significant_transcripts_p_adj2 <- rownames(significant_transcripts_p_adj2)


significant_transcripts_p_adj3 <- subset(results_table_AD_vs_Control, adj.P.Val < 0.001)

significant_transcripts_p_adj3 <- rownames(significant_transcripts_p_adj3)



# Write the transcript names to a text file
write.table(significant_transcripts_p_val1, file = "transcript_names_AD_vs_Control_limma_p_val1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_transcripts_p_val2, file = "transcript_names_AD_vs_Control_limma_p_val2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_transcripts_p_val3, file = "transcript_names_AD_vs_Control_limma_p_val3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


write.table(significant_transcripts_p_adj1, file = "transcript_names_AD_vs_Control_limma_p_adj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_transcripts_p_adj2, file = "transcript_names_AD_vs_Control_limma_p_adj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_transcripts_p_adj3, file = "transcript_names_AD_vs_Control_limma_p_adj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



# Gene results

# mart_export is a file that includes Gene.stable.IDs, Gene.stable.ID.versions, Transcript.stable.IDs Transcript.stable.ID.versions
# downloaded from Ensembl Biomart

data <- read.csv("mart_export.txt")   

transcripts <- rownames(cpm_filtered_count_out)

# Filter the transcripts in ‘transcripts’ from the data
filtered_data <- data[data$Transcript.stable.ID.version %in% transcripts, ]

# Select the relevant columns
filtered_data <- filtered_data[, c("Transcript.stable.ID.version", "Gene.stable.ID")]



library(dplyr)

genes_1 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts_p_val1) %>%
  pull(Gene.stable.ID) %>% unique()


genes_2 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts_p_val2) %>%
  pull(Gene.stable.ID) %>% unique()

genes_3 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts_p_val3) %>%
  pull(Gene.stable.ID) %>% unique()


genes_4 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts_p_adj1) %>%
  pull(Gene.stable.ID) %>% unique()

genes_5 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts_p_adj2) %>%
  pull(Gene.stable.ID) %>% unique()

genes_6 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts_p_adj3) %>%
  pull(Gene.stable.ID) %>% unique()



# Write the gene names to text files
write.table(genes_1, file = "gene_names_AD_vs_Control_limma_p_val1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_2, file = "gene_names_AD_vs_Control_limma_p_val2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_3, file = "gene_names_AD_vs_Control_limma_p_val3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


write.table(genes_4, file = "gene_names_AD_vs_Control_limma_p_adj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_5, file = "gene_names_AD_vs_Control_limma_p_adj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_6, file = "gene_names_AD_vs_Control_limma_p_adj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)




# AD-AsymAD

# Create the contrasts for AD vs AsymAD
contrasts_AD_vs_AsymAD <- makeContrasts(
  AD_vs_AsymAD = AD - AsymAD,
  levels = design
)

# Apply the contrasts for AD vs AsymAD
contrast_fit_AD_vs_AsymAD <- contrasts.fit(fit, contrasts_AD_vs_AsymAD)

# Run eBayes on the contrast fit
contrast_ebayes_AD_vs_AsymAD <- eBayes(contrast_fit_AD_vs_AsymAD)

# Retrieve the results for AD vs AsymAD
results_table_AD_vs_AsymAD <- topTable(contrast_ebayes_AD_vs_AsymAD, coef = "AD_vs_AsymAD", number = nrow(fitteddata))



# Transcript results

# Filter genes based on adjusted p-value threshold (e.g., 0.05)
significant_transcripts2_p_val1 <- subset(results_table_AD_vs_AsymAD, P.Value < 0.05)

# Extract gene names from the results table
significant_transcripts2_p_val1 <- rownames(significant_transcripts2_p_val1)


significant_transcripts2_p_val2 <- subset(results_table_AD_vs_AsymAD, P.Value < 0.01)

significant_transcripts2_p_val2 <- rownames(significant_transcripts2_p_val2)


significant_transcripts2_p_val3 <- subset(results_table_AD_vs_AsymAD, P.Value < 0.001)

significant_transcripts2_p_val3 <- rownames(significant_transcripts2_p_val3)



significant_transcripts2_p_adj1 <- subset(results_table_AD_vs_AsymAD, adj.P.Val < 0.05)

significant_transcripts2_p_adj1 <- rownames(significant_transcripts2_p_adj1)


significant_transcripts2_p_adj2 <- subset(results_table_AD_vs_AsymAD, adj.P.Val < 0.01)

significant_transcripts2_p_adj2 <- rownames(significant_transcripts2_p_adj2)


significant_transcripts2_p_adj3 <- subset(results_table_AD_vs_AsymAD, adj.P.Val < 0.001)

significant_transcripts2_p_adj3 <- rownames(significant_transcripts2_p_adj3)



# Write the transcript names to a text file
write.table(significant_transcripts2_p_val1, file = "transcript_names_AD_vs_AsymAD_limma_p_val1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_transcripts2_p_val2, file = "transcript_names_AD_vs_AsymAD_limma_p_val2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_transcripts2_p_val3, file = "transcript_names_AD_vs_AsymAD_limma_p_val3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


write.table(significant_transcripts2_p_adj1, file = "transcript_names_AD_vs_AsymAD_limma_p_adj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_transcripts2_p_adj2, file = "transcript_names_AD_vs_AsymAD_limma_p_adj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_transcripts2_p_adj3, file = "transcript_names_AD_vs_AsymAD_limma_p_adj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



# Gene results

genes_7 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts2_p_val1) %>%
  pull(Gene.stable.ID) %>% unique()


genes_8 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts2_p_val2) %>%
  pull(Gene.stable.ID) %>% unique()

genes_9 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts2_p_val3) %>%
  pull(Gene.stable.ID) %>% unique()


genes_10 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts2_p_adj1) %>%
  pull(Gene.stable.ID) %>% unique()

genes_11 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts2_p_adj2) %>%
  pull(Gene.stable.ID) %>% unique()

genes_12 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts2_p_adj3) %>%
  pull(Gene.stable.ID) %>% unique()



# Write the gene names to a text file
write.table(genes_7, file = "gene_names_AD_vs_AsymAD_limma_p_val1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_8, file = "gene_names_AD_vs_AsymAD_limma_p_val2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_9, file = "gene_names_AD_vs_AsymAD_limma_p_val3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


write.table(genes_10, file = "gene_names_AD_vs_AsymAD_limma_p_adj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_11, file = "gene_names_AD_vs_AsymAD_limma_p_adj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_12, file = "gene_names_AD_vs_AsymAD_limma_p_adj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



# Control-AsymAD

# Create the contrasts for AsymAD vs Control
contrasts_AsymAD_vs_Control <- makeContrasts(
  AsymAD_vs_Control = AsymAD - Control,
  levels = design
)


# Apply the contrasts for AsymAD vs Control
contrast_fit_AsymAD_vs_Control <- contrasts.fit(fit, contrasts_AsymAD_vs_Control)

# Run eBayes on the contrast fit
contrast_ebayes_AsymAD_vs_Control <- eBayes(contrast_fit_AsymAD_vs_Control)

# Retrieve the results for AsymAD vs Control
results_table_AsymAD_vs_Control <- topTable(contrast_ebayes_AsymAD_vs_Control, coef = "AsymAD_vs_Control", number = nrow(fitteddata))


# Transcript results

# Filter genes based on adjusted p-value threshold (e.g., 0.05)
significant_transcripts3_p_val1 <- subset(results_table_AsymAD_vs_Control, P.Value < 0.05)

# Extract gene names from the results table
significant_transcripts3_p_val1 <- rownames(significant_transcripts3_p_val1)


significant_transcripts3_p_val2 <- subset(results_table_AsymAD_vs_Control, P.Value < 0.01)

significant_transcripts3_p_val2 <- rownames(significant_transcripts3_p_val2)


significant_transcripts3_p_val3 <- subset(results_table_AsymAD_vs_Control, P.Value < 0.001)

significant_transcripts3_p_val3 <- rownames(significant_transcripts3_p_val3)



significant_transcripts3_p_adj1 <- subset(results_table_AsymAD_vs_Control, adj.P.Val < 0.05)

significant_transcripts3_p_adj1 <- rownames(significant_transcripts3_p_adj1)


significant_transcripts3_p_adj2 <- subset(results_table_AsymAD_vs_Control, adj.P.Val < 0.01)

significant_transcripts3_p_adj2 <- rownames(significant_transcripts3_p_adj2)


significant_transcripts3_p_adj3 <- subset(results_table_AsymAD_vs_Control, adj.P.Val < 0.001)

significant_transcripts3_p_adj3 <- rownames(significant_transcripts3_p_adj3)



# Write the transcript names to a text file
write.table(significant_transcripts3_p_val1, file = "transcript_names_AsymAD_vs_Control_limma_p_val1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_transcripts3_p_val2, file = "transcript_names_AsymAD_vs_Control_limma_p_val2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_transcripts3_p_val3, file = "transcript_names_AsymAD_vs_Control_limma_p_val3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


write.table(significant_transcripts3_p_adj1, file = "transcript_names_AsymAD_vs_Control_limma_p_adj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_transcripts3_p_adj2, file = "transcript_names_AsymAD_vs_Control_limma_p_adj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_transcripts3_p_adj3, file = "transcript_names_AsymAD_vs_Control_limma_p_adj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



# Gene results


genes_13 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts3_p_val1) %>%
  pull(Gene.stable.ID) %>% unique()

genes_14 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts3_p_val2) %>%
  pull(Gene.stable.ID) %>% unique()

genes_15 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts3_p_val3) %>%
  pull(Gene.stable.ID) %>% unique()


genes_16 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts3_p_adj1) %>%
  pull(Gene.stable.ID) %>% unique()

genes_17 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts3_p_adj2) %>%
  pull(Gene.stable.ID) %>% unique()

genes_18 <- filtered_data %>%
  filter(Transcript.stable.ID.version %in% significant_transcripts3_p_adj3) %>%
  pull(Gene.stable.ID) %>% unique()



# Write the gene names to a text file
write.table(genes_13, file = "gene_names_AsymAD_vs_Control_limma_p_val1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_14, file = "gene_names_AsymAD_vs_Control_limma_p_val2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_15, file = "gene_names_AsymAD_vs_Control_p_val3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


write.table(genes_16, file = "gene_names_AsymAD_vs_Control_limma_p_adj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_17, file = "gene_names_AsymAD_vs_Control_limma_p_adj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(genes_18, file = "gene_names_AsymAD_vs_Control_limma_p_adj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
