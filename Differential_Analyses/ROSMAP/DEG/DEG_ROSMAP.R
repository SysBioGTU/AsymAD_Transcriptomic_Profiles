
# DEG

load("/Volumes/Untitled/AsymAD_Project_Rumeysa/Differential Analyses/ROSMAP/DEG/rosmapmetadata.RData") # load the metadata


# Defining groups according to the cogdx and ceradsc scores in the metadata

# Consensus Control

cons_cont_4 <- metadata[(metadata$cogdx == 1 | metadata$cogdx == 2 | metadata$cogdx == 3 ) & (metadata$ceradsc == 3 | metadata$ceradsc == 4), ]


# Consensus AD 

cons_AD <- metadata[(metadata$cogdx == 4 | metadata$cogdx == 5) & (metadata$ceradsc == 1 | metadata$ceradsc == 2), ]


# AsymAD

cons_AsymAD <- metadata[metadata$cogdx == 1 & (metadata$ceradsc == 1 | metadata$ceradsc == 2), ]



# Choose the groups from countdata

a <- cons_cont_4$rnaseq_id  #for the control IDs

b <- cons_AD$rnaseq_id  #for the AD IDs

c <- cons_AsymAD$rnaseq_id  #for the AsymAD IDs


# Define the order of columns based on 'a', 'b', and 'c'
order_of_columns <- c(a, b, c)

# Select columns in the desired order
gene_counts <- gene_counts[, order_of_columns]

# Convert each column of the data frame to numeric

for (i in 1:ncol(combined_count_data)) {
  combined_count_data[, i] <- as.numeric(combined_count_data[, i])
}

# Make the Ensembl gene IDs as row names
row.names(combined_count_data) <- make.unique(rosmapcountdata$ensembl_gene_id)


# Calculate CPM values and apply cpm treshold

library(edgeR)

myCPM <- cpm(gene_counts,log=F)  
thresh <- myCPM > 0.1  # library size of ROSMAP ~50M. 5 count/ 50 = 0.1 cpm treshold
keep <- rowSums(thresh) >= 92  # Number of samples of the condition with the least number of samples  
summary(keep) 
cpm_filtered_count <- gene_counts[keep,] 


# Create a metadata table with sample information
col_Data <- data.frame(sample_id = colnames(cpm_filtered_count),
                       group = c(rep("Control", length(a)),
                                 rep("AD", length(b)),
                                 rep("AsymAD", length(c))))

# Set the row names to match the sample IDs
row.names(col_Data) <- col_Data$sample_id

# making sure the row names in colData matches to column names in counts_data
all(colnames(cpm_filtered_count) %in% rownames(col_Data))  # yes

# are they in the same order?
all(colnames(cpm_filtered_count) == rownames(col_Data)) # yes


# Outlier separation--------------------------------------------------------------------------------------------

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

normalized_counts <- counts(dds, normalized = TRUE) 



# PCA Plot with the groups------------------------------------------------------------------------------

# After the codes here for PCA, outliers were obtained and went back to the outlier separation step above. 

#library(ggplot2)

# Perform PCA on the normalized counts
#pca_result <- prcomp(t(normalized_counts), scale = TRUE)

# Create a data frame for visualization
#pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], group = col_Data$group)

# Create the PCA plot
#ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
#  geom_point() +
#  labs(x = "Principal Component 1", y = "Principal Component 2") +
#  ggtitle("ROSMAP PCA Plot")

#--------------------------------------------------------------------------------------------------------



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

fitteddata = matrix(0, 22303, 485)  # The size of the count matrix

cins = sex.num
cins[cins=="female"]=0
cins[cins=="male"]=1
cins= as.numeric(cins)

state = condition.num
state[state=="AD"]=0
state[state=="Control"]=1
state[state=="AsymAD"]=2
state = as.numeric(state)


genenames <- character(22303)

for (j in 1:485){
  
  for (i in 1:22303) {
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


#-----------------------------------------------------------------------------------------------------------------------------------------------------


# Limma 

library(limma)


# AD vs Control

# Convert the group information into a numeric matrix
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


# Filter genes based on adjusted p-value threshold (e.g., 0.05)
significant_genes_p_val1 <- subset(results_table_AD_vs_Control, P.Value < 0.05)

# Extract gene names from the results table
significant_genes_p_val1 <- rownames(significant_genes_p_val1)


significant_genes_p_val2 <- subset(results_table_AD_vs_Control, P.Value < 0.01)

significant_genes_p_val2 <- rownames(significant_genes_p_val2)


significant_genes_p_val3 <- subset(results_table_AD_vs_Control, P.Value < 0.001)

significant_genes_p_val3 <- rownames(significant_genes_p_val3)



significant_genes_p_adj1 <- subset(results_table_AD_vs_Control, adj.P.Val < 0.05)

significant_genes_p_adj1 <- rownames(significant_genes_p_adj1)


significant_genes_p_adj2 <- subset(results_table_AD_vs_Control, adj.P.Val < 0.01)

significant_genes_p_adj2 <- rownames(significant_genes_p_adj2)


significant_genes_p_adj3 <- subset(results_table_AD_vs_Control, adj.P.Val < 0.001)

significant_genes_p_adj3 <- rownames(significant_genes_p_adj3)



# Write the gene names to a text file
write.table(significant_genes_p_val1, file = "gene_names_AD_vs_Control_limma_p_val1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(significant_genes_p_val2, file = "gene_names_AD_vs_Control_limma_p_val2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(significant_genes_p_val3, file = "gene_names_AD_vs_Control_limma_p_val3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_genes_p_adj1, file = "gene_names_AD_vs_Control_limma_p_adj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(significant_genes_p_adj2, file = "gene_names_AD_vs_Control_limma_p_adj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(significant_genes_p_adj3, file = "gene_names_AD_vs_Control_limma_p_adj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


# AsymAD vs Control


# Create the contrasts for AsymAD vs Control
contrasts_AsymAD_vs_Control <- makeContrasts(
  AsymAD_vs_Control = AsymAD - Control,
  levels = design
)

# The linear model was fitted before as "fit" in the 237th line

# Apply the contrasts for AsymAD vs Control
contrast_fit_AsymAD_vs_Control <- contrasts.fit(fit, contrasts_AsymAD_vs_Control)

# Run eBayes on the contrast fit
contrast_ebayes_AsymAD_vs_Control <- eBayes(contrast_fit_AsymAD_vs_Control)

# Retrieve the results for AsymAD vs Control
results_table_AsymAD_vs_Control <- topTable(contrast_ebayes_AsymAD_vs_Control, coef = "AsymAD_vs_Control", number = nrow(fitteddata))



significant_genes2_p_val1 <- subset(results_table_AsymAD_vs_Control, P.Value < 0.05)

significant_genes2_p_val1 <- rownames(significant_genes2_p_val1)


significant_genes2_p_val2 <- subset(results_table_AsymAD_vs_Control, P.Value < 0.01)

significant_genes2_p_val2 <- rownames(significant_genes2_p_val2)


significant_genes2_p_val3 <- subset(results_table_AsymAD_vs_Control, P.Value < 0.001)

significant_genes2_p_val3 <- rownames(significant_genes2_p_val3)




significant_genes2_p_adj1 <- subset(results_table_AsymAD_vs_Control, adj.P.Val < 0.05)

significant_genes2_p_adj1 <- rownames(significant_genes2_p_adj1)


significant_genes2_p_adj2 <- subset(results_table_AsymAD_vs_Control, adj.P.Val < 0.01)

significant_genes2_p_adj2 <- rownames(significant_genes2_p_adj2)


significant_genes2_p_adj3 <- subset(results_table_AsymAD_vs_Control, adj.P.Val < 0.001)

significant_genes2_p_adj3 <- rownames(significant_genes2_p_adj3)



# Write the gene names to a text file
write.table(significant_genes2_p_val1, file = "gene_names_AsymAD_vs_Control_limma_p_val1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(significant_genes2_p_val2, file = "gene_names_AsymAD_vs_Control_limma_p_val2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(significant_genes2_p_val3, file = "gene_names_AsymAD_vs_Control_limma_p_val3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_genes2_p_adj1, file = "gene_names_AsymAD_vs_Control_limma_p_adj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(significant_genes2_p_adj2, file = "gene_names_AsymAD_vs_Control_limma_p_adj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(significant_genes2_p_adj3, file = "gene_names_AsymAD_vs_Control_limma_p_adj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



# AD vs. AsymAD

# Create the contrasts for AD vs AsymAD
contrasts_AD_vs_AsymAD <- makeContrasts(
  AD_vs_AsymAD = AD - AsymAD,
  levels = design
)

# The linear model was fitted before as "fit" in the 237th line

# Apply the contrasts for AD vs AsymAD
contrast_fit_AD_vs_AsymAD <- contrasts.fit(fit, contrasts_AD_vs_AsymAD)

# Run eBayes on the contrast fit
contrast_ebayes_AD_vs_AsymAD <- eBayes(contrast_fit_AD_vs_AsymAD)

# Retrieve the results for AD vs AsymAD
results_table_AD_vs_AsymAD <- topTable(contrast_ebayes_AD_vs_AsymAD, coef = "AD_vs_AsymAD", number = nrow(fitteddata))



significant_genes3_p_val1 <- subset(results_table_AD_vs_AsymAD, P.Value < 0.05)

significant_genes3_p_val1 <- rownames(significant_genes3_p_val1)


significant_genes3_p_val2 <- subset(results_table_AD_vs_AsymAD, P.Value < 0.01)

significant_genes3_p_val2 <- rownames(significant_genes3_p_val2)


significant_genes3_p_val3 <- subset(results_table_AD_vs_AsymAD, P.Value < 0.001)

significant_genes3_p_val3 <- rownames(significant_genes3_p_val3)



significant_genes3_p_adj1 <- subset(results_table_AD_vs_AsymAD, adj.P.Val < 0.05)

significant_genes3_p_adj1 <- rownames(significant_genes3_p_adj1)


significant_genes3_p_adj2 <- subset(results_table_AD_vs_AsymAD, adj.P.Val < 0.01)

significant_genes3_p_adj2 <- rownames(significant_genes3_p_adj2)


significant_genes3_p_adj3 <- subset(results_table_AD_vs_AsymAD, adj.P.Val < 0.001)

significant_genes3_p_adj3 <- rownames(significant_genes3_p_adj3)



# Write the gene names to a text file
write.table(significant_genes3_p_val1, file = "gene_names_AD_vs_AsymAD_limma_p_val1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(significant_genes3_p_val2, file = "gene_names_AD_vs_AsymAD_limma_p_val2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(significant_genes3_p_val3, file = "gene_names_AD_vs_AsymAD_limma_p_val3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(significant_genes3_p_adj1, file = "gene_names_AD_vs_AsymAD_limma_p_adj1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(significant_genes3_p_adj2, file = "gene_names_AD_vs_AsymAD_limma_p_adj2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(significant_genes3_p_adj3, file = "gene_names_AD_vs_AsymADl_limma_p_adj3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


