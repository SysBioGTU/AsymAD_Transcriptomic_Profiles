
library(readxl)
library(dplyr)
library(writexl)

# Read the first file

KPM_184gene_results <- read_excel("KPM_184gene_results.xlsx")

df1 <- KPM_184_gene_results 


# Read the second file

candidate_genes <- read_excel("candidate_genes.xlsx")

df2 <- candidate_genes

df2 <- df2 %>%
  rename(shared_name = `Shared Genes`, Datasets_Scores = `Datasets & Scores`)


# Fill the resource column by matching 'shared_name' column

df1 <- df1 %>%
  left_join(df2, by = "shared_name") %>%
  mutate(resource = ifelse(is.na(Datasets_Scores), resource, Datasets_Scores)) %>%
  select(shared_name, resource)

# Save the new file
write_xlsx(df1, "KPM_184_gene_results.xlsx")

