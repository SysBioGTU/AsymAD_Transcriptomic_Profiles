

# Fill the resource column based on shared_name and the DGE, DTE, and DTU columns

data <- Common_Genes_Whole_List_ROSMAP

library(dplyr)

library(stringr)

library(writexl)

data <- data %>%
  mutate(resource = str_trim(paste(
    ifelse(genes %in% DGE, "DGE", ""),
    ifelse(genes %in% DTE, "DTE", ""),
    ifelse(genes %in% DTU, "DTU", ""),
    sep = " "
  ))) %>%
  mutate(resource = na_if(resource, ""))



write_xlsx(data, "Common_Genes_Whole_List_ROSMAP.xlsx")



data_2 <- Common_Genes_Whole_List_MSBB

library(dplyr)

library(stringr)

library(writexl)

data_2 <- data_2 %>%
  mutate(resource = str_trim(paste(
    ifelse(genes %in% DGE, "DGE", ""),
    ifelse(genes %in% DTE, "DTE", ""),
    ifelse(genes %in% DTU, "DTU", ""),
    sep = " "
  ))) %>%
  mutate(resource = na_if(resource, ""))



write_xlsx(data, "Common_Genes_Whole_List_MSBB.xlsx")




filtered_genes_all_rosmap <- data %>% filter(resource %in% "DGE DTE DTU") %>% select(genes)

filtered_genes_all_msbb <- data_2 %>% filter(resource %in% "DGE DTE DTU") %>% select(genes)


filtered_genes_DGE_DTE_rosmap <- data %>% filter(resource %in% "DGE DTE") %>% select(genes)

filtered_genes_DGE_DTE_msbb <- data_2 %>% filter(resource %in% "DGE DTE") %>% select(genes)


filtered_genes_DTE_DTU_rosmap <- data %>% filter(resource %in% "DTE DTU") %>% select(genes)

filtered_genes_DTE_DTU_msbb <- data_2 %>% filter(resource %in% "DTE DTU") %>% select(genes)


filtered_genes_DGE_DTU_rosmap <- data %>% filter(resource %in% "DGE DTU") %>% select(genes)

filtered_genes_DGE_DTU_msbb <- data_2 %>% filter(resource %in% "DGE DTU") %>% select(genes)


filtered_genes_DGE_rosmap <- data %>% filter(resource %in% "DGE") %>% select(genes)

filtered_genes_DGE_msbb <- data_2 %>% filter(resource %in% "DGE") %>% select(genes)


filtered_genes_DTE_rosmap <- data %>% filter(resource %in% "DTE") %>% select(genes)

filtered_genes_DTE_msbb <- data_2 %>% filter(resource %in% "DTE") %>% select(genes)


filtered_genes_DTU_rosmap <- data %>% filter(resource %in% "DTU") %>% select(genes)

filtered_genes_DTU_msbb <- data_2 %>% filter(resource %in% "DTU") %>% select(genes)




cat(filtered_genes_DTE_DTU_rosmap$genes, sep = "\n")

cat(filtered_genes_DTE_DTU_msbb$genes, sep = "\n")



cat(filtered_genes_DGE_DTE_rosmap$genes, sep = "\n")

cat(filtered_genes_DGE_DTE_msbb$genes, sep = "\n")




data <- data[, -c(3,4,5)]

data_2 <- data_2[, -c(3,4,5)]


data <- data %>%
  mutate(dataset = "ROSMAP")

data_2 <- data_2 %>%
  mutate(dataset = "MSBB")


# Combining the two datasets
combined_data <- bind_rows(data, data_2)




combined_data <- combined_data %>%
  mutate(
    scores = case_when(
      # Scoring for ROSMAP
      dataset == "ROSMAP" & resource %in% c("DGE", "DTE", "DTU") ~ 1,
      dataset == "ROSMAP" & resource %in% c("DGE DTE", "DTE DTU", "DGE DTU") ~ 2,
      dataset == "ROSMAP" & resource == "DGE DTE DTU" ~ 3,
      
      # Scoring for MSBB 
      dataset == "MSBB" & resource %in% c("DGE", "DTE", "DTU") ~ 1,
      dataset == "MSBB" & resource %in% c("DGE DTE", "DTE DTU", "DGE DTU") ~ 2,
      dataset == "MSBB" & resource == "DGE DTE DTU" ~ 3,
      
      TRUE ~ NA_real_
    )
  )




# Condition 1: R3-M1

# Select gene IDs with scores == 3 in ROSMAP
rosmap_con1 <- combined_data %>%
  filter(dataset == "ROSMAP" & scores == 3) %>%
  select(genes)

# Select gene IDs with scores == 1 in MSBB
msbb_con1 <- combined_data %>%
  filter(dataset == "MSBB" & scores == 1) %>%
  select(genes)

# Find the common gene IDs
common_con1 <- intersect(rosmap_con1$genes, msbb_con1$genes)



# Define the conditions

conditions <- list(
  "R3-M1" = list(rosmap_score = 3, msbb_score = 1),
  "R3-M2" = list(rosmap_score = 3, msbb_score = 2),
  "R3-M3" = list(rosmap_score = 3, msbb_score = 3),
  "R2-M1" = list(rosmap_score = 2, msbb_score = 1),
  "R2-M2" = list(rosmap_score = 2, msbb_score = 2),
  "R2-M3" = list(rosmap_score = 2, msbb_score = 3),
  "R1-M1" = list(rosmap_score = 1, msbb_score = 1),
  "R1-M2" = list(rosmap_score = 1, msbb_score = 2),
  "R1-M3" = list(rosmap_score = 1, msbb_score = 3)
)

# Create a list to store the common genes
results <- list()

# Iterate over each condition using a loop

for (cond in names(conditions)) {
  # Get the conditions for the current case
  rosmap_score <- conditions[[cond]]$rosmap_score
  msbb_score <- conditions[[cond]]$msbb_score
  
  # ROSMAP için genleri filtrele
  rosmap_genes <- combined_data %>%
    filter(dataset == "ROSMAP" & scores == rosmap_score) %>%
    select(genes)
  
  # Filter the genes for ROSMAP
  msbb_genes <- combined_data %>%
    filter(dataset == "MSBB" & scores == msbb_score) %>%
    select(genes)
  
  # Find the common genes
  common_genes <- intersect(rosmap_genes$genes, msbb_genes$genes)
  
  # Add the results to the list
  results[[cond]] <- common_genes
}

results


# Write the common genes found for each condition to .csv files

for (cond in names(results)) {
  write.csv(results[[cond]], paste0(cond, "_common_genes.csv"), row.names = FALSE)
}



