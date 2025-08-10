

library(readxl)

biogrid_interactions <- read_excel("~/Desktop/biogrid_interactions.xlsx")

target_genes <- read.table("target_genes.txt")

target_genes <- target_genes[[1]]  # ilk sütunu vektör olarak çek

target_genes <- unique(target_genes)

target_genes <- target_genes[-697]

write.table(target_genes, "target_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)



library(dplyr) 

filtered_interactions <- biogrid_interactions %>% 
  filter(A %in% genelist_185 & C %in% target_genes)


filtered_interactions_2 <- biogrid_interactions %>% 
  filter(C %in% genelist_185 & A %in% target_genes)

merged_ineractions <- bind_rows(filtered_interactions, filtered_interactions_2)

merged_ineractions <- merged_ineractions %>%
  mutate(
    match = ifelse(A %in% genelist_185, 1, 0)
  )

sum(merged_ineractions$match == 1)

merged_ineractions_unique <- unique(merged_ineractions)

write.table(merged_ineractions, "merged_ineractions.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


n_important_genes_with_interaction <- merged_ineractions_unique %>%
  filter(A %in% genelist_185) %>%
  distinct(A) %>%
  nrow()


n_DEGs_with_interaction <- merged_ineractions_unique %>%
  filter(C %in% target_genes) %>%
  distinct(C) %>%
  nrow()






