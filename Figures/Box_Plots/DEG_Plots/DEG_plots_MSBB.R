

# DGE plots for ADAMTS2
# ENSG00000087116



gene_of_interest <- "ENSG00000087116" 


# AD vs AsymAD countdata

count_AD_Asym <- fitteddata[, c(183:485)]

col_Data_3 <- col_Data_out[ c(183:485), ]


library(dplyr)

library(tibble)


# Filter data
filtered_data_1 <- count_AD_Asym[rownames(count_AD_Asym) %in% gene_of_interest, ]

filtered_data_1 <- as.data.frame(filtered_data_1)


# logsuz

library(tidyverse)

# Convert rownames of filtered_data_GABARAPL2 to a column if needed
filtered_data_1 <- filtered_data_1 %>%
  rownames_to_column(var = "Sample")  #!!!

filtered_data_1 <- as.data.frame(filtered_data_1)

library(dplyr)

library(tidyr)


# Reshape the data to long format for ggplot2
long_data_1 <- filtered_data_1 %>%
  pivot_longer(cols = -Sample, names_to = "sample_id", values_to = "Expression")

# Join long_data_GABARAPL2 with col_Data_3 to get the group information
long_data_1 <- long_data_1 %>%
  left_join(col_Data_3, by = "sample_id")




library(ggsignif)
library(ggplot2)



# Create the plot
ggplot(long_data_1, aes(x = group, y = Expression, color = group, fill = group, shape = group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(size = 2, alpha = 0.7, position = position_jitter(width = 0.2, height = 0)) +
  theme_minimal() +
  labs(title = "ADAMTS2 Expressions in AD vs AsymAD (ROSMAP)",
       x = "Condition", y = "Expression") +
  scale_color_manual(values = c("AD" = "blue", "AsymAD" = "purple")) +
  scale_fill_manual(values = c("AD" = "blue", "AsymAD" = "purple")) +
  scale_shape_manual(values = c("AD" = 16, "AsymAD" = 17)) + # Circles for AD, triangles for AsymAD
  geom_signif(comparisons = list(c("AD", "AsymAD")), 
              map_signif_level = TRUE, textsize = 3, test = "t.test", color = "black") + # Add significance line
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.grid.major = element_blank(), # Remove large grids
        panel.grid.minor = element_blank(),  # Remove small grids
        axis.line = element_line(color = "black") 
  )



path <- "H:/plots_last/DGE_plots"


ggsave(filename = file.path(path, "ADAMTS2_Expressions_AD_vs_AsymAD_ROSMAP.tiff"),
       plot = last_plot(),
       width = 9, height = 7, device = "tiff", dpi = 600)





# MSBB

# Gene: GNG4 ENSG00000168243

gene_of_interest <- "ENSG00000168243" 


# AD vs AsymAD countdata

count_AD_Asym <- fitteddata[, c(21:181)]

col_Data_3 <- col_Data[ c(21:181), ]


library(dplyr)

library(tibble)


# Filter data

filtered_data_2 <- count_AD_Asym[rownames(count_AD_Asym) %in% gene_of_interest, ]

filtered_data_2 <- as.data.frame(filtered_data_2)



#logsuz


# Convert rownames of filtered_data_GABARAPL2 to a column if needed
filtered_data_2 <- filtered_data_2 %>%
  rownames_to_column(var = "Sample")  #!!!

filtered_data_2 <- as.data.frame(filtered_data_2)


library(dplyr)

library(tidyr)


# Reshape the data to long format for ggplot2
long_data_2 <- filtered_data_2 %>%
  pivot_longer(cols = -Sample, names_to = "sample_id", values_to = "Expression")

# Join long_data_GABARAPL2 with col_Data_3 to get the group information
long_data_2 <- long_data_2 %>%
  left_join(col_Data_3, by = "sample_id")




library(ggsignif)
library(ggplot2)

# Create the plot
ggplot(long_data_2, aes(x = group, y = Expression, color = group, fill = group, shape = group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(size = 2, alpha = 0.7, position = position_jitter(width = 0.2, height = 0)) +
  theme_minimal() +
  labs(title = "GNG4 Expressions in AD vs AsymAD (MSBB)",
       x = "Condition", y = "Expression") +
  scale_color_manual(values = c("AD" = "skyblue", "AsymAD" = "darkblue")) +
  scale_fill_manual(values = c("AD" = "skyblue", "AsymAD" = "darkblue")) +
  scale_shape_manual(values = c("AD" = 16, "AsymAD" = 17)) + # Circles for AD, triangles for AsymAD
  geom_signif(comparisons = list(c("AD", "AsymAD")), 
              map_signif_level = TRUE, textsize = 3, test = "t.test", color = "black") + # Add significance line
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.grid.major = element_blank(), # Remove large grids
        panel.grid.minor = element_blank(), # Remove small grids
        axis.line = element_line(color = "black") 
  )


path <- "H:/plots_last/DGE_plots"

ggsave(filename = file.path(path, "GNG4_Expressions_AD_vs_AsymAD_MSBB.tiff"),
       plot = last_plot(),
       width = 9, height = 7, device = "tiff", dpi = 600)





# Gene: MRPL1 ENSG00000169288

gene_of_interest <- "ENSG00000169288" 





# AD vs AsymAD countdata

count_AD_Asym <- fitteddata[, c(21:181)]

col_Data_3 <- col_Data[ c(21:181), ]


library(dplyr)

library(tibble)


# Filter data

filtered_data_2 <- count_AD_Asym[rownames(count_AD_Asym) %in% gene_of_interest, ]

filtered_data_2 <- as.data.frame(filtered_data_2)




# Convert rownames of filtered_data_GABARAPL2 to a column if needed
filtered_data_2 <- filtered_data_2 %>%
  rownames_to_column(var = "Sample")  #!!!

filtered_data_2 <- as.data.frame(filtered_data_2)


library(dplyr)

library(tidyr)


# Reshape the data to long format for ggplot2
long_data_2 <- filtered_data_2 %>%
  pivot_longer(cols = -Sample, names_to = "sample_id", values_to = "Expression")

# Join long_data_GABARAPL2 with col_Data_3 to get the group information
long_data_2 <- long_data_2 %>%
  left_join(col_Data_3, by = "sample_id")




library(ggsignif)
library(ggplot2)

# Create the plot
ggplot(long_data_2, aes(x = group, y = Expression, color = group, fill = group, shape = group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(size = 2, alpha = 0.7, position = position_jitter(width = 0.2, height = 0)) +
  theme_minimal() +
  labs(title = "MRPL1 Expressions in AD vs AsymAD (MSBB)",
       x = "Condition", y = "Expression") +
  scale_color_manual(values = c("AD" = "skyblue", "AsymAD" = "darkblue")) +
  scale_fill_manual(values = c("AD" = "skyblue", "AsymAD" = "darkblue")) +
  scale_shape_manual(values = c("AD" = 16, "AsymAD" = 17)) + # Circles for AD, triangles for AsymAD
  geom_signif(comparisons = list(c("AD", "AsymAD")), 
              map_signif_level = TRUE, textsize = 3, test = "t.test", color = "black") + # Add significance line
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.grid.major = element_blank(), # Remove large grids
        panel.grid.minor = element_blank(), # Remove small grids
        axis.line = element_line(color = "black") 
  )


path <- "H:/plots_last/DGE_plots"

ggsave(filename = file.path(path, "MRPL1_Expressions_AD_vs_AsymAD_MSBB.tiff"),
       plot = last_plot(),
       width = 9, height = 7, device = "tiff", dpi = 600)
