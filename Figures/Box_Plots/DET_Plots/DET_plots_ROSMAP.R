
library(dplyr)
library(tibble)
library(tidyverse)
library(ggsignif)
library(ggplot2)



# Gene: GNG4 ENSG00000168243

transcripts_of_interest <- c("ENST00000391854.7", "ENST00000450593.5", "ENST00000366598.8", "ENST00000366597.5", "ENST00000484517.2")


x <- rownames(fitteddata)


transcripts_of_interest %in% x

enst_598.8 <- combined_count_out["ENST00000366598.8",]

enst_517.2 <- combined_count_out["ENST00000484517.2",]


fitteddata_2 <- rbind(fitteddata, enst_598.8, enst_517.2 )



fitteddata_with_rownames <- data.frame(RowNames = row.names(fitteddata_2), fitteddata_2, check.names = FALSE)



# AD vs AsymAD countdata

count_AD_Asym <- fitteddata_with_rownames[, c(1, 184:486)]

col_Data_2 <- col_Data_out[ 183:485, ]


# Filter data

filtered_data_1 <- count_AD_Asym[count_AD_Asym$RowNames %in% transcripts_of_interest, ]


library(tidyverse)

# Reshape the data to long format for ggplot2
long_data_1 <- filtered_data_1 %>%
  pivot_longer(cols = -RowNames, names_to = "Sample", values_to = "Expression")

# Add group information from col_Data_3
col_Data_2 <- col_Data_2 %>%
  mutate(Sample = col_Data_2$sample_id)

long_data_1 <- long_data_1 %>%
  left_join(col_Data_2, by = c("Sample"))



# Create the plot
ggplot(long_data_1, aes(x = group, y = Expression, color = group, shape = group, fill = group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.7, position = position_jitter(width = 0.2)) +
  theme_minimal() +
  labs(title = "Transcripts of GNG4 Expressions in AD vs AsymAD (ROSMAP)",
       x = "Group", y = "Expression") +
  scale_color_manual(values = c("AD" = "blue", "AsymAD" = "purple")) +
  scale_fill_manual(values = c("AD" = "blue", "AsymAD" = "purple")) +
  geom_signif(comparisons = list(c("AD", "AsymAD")), 
              map_signif_level = TRUE, textsize = 3, test = "t.test", color = "black") +
  facet_wrap(~ RowNames, scales = "free") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +  # Ensure 4 breaks on y-axis
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 4), # Adjust Y-axis tick label size
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black")  # Add black axis lines
  )



path <- "H:/plots_last/DET_plots"

ggsave(filename = file.path(path, "GNG4_Expressions_AD_vs_AsymAD_ROSMAP_DET.tiff"),
       plot = last_plot(),
       width = 9, height = 7, device = "tiff", dpi = 600)




# Gene: ADAMTS2 ENSG00000087116


transcripts_of_interest <- c("ENST00000251582.12", "ENST00000518335.3", "ENST00000274609.5", "ENST00000698889.1", "ENST00000522937.1", "ENST00000523450.1")

x <- rownames(fitteddata)

transcripts_of_interest %in% x


enst_609.5 <- combined_count_out["ENST00000274609.5",]

enst_937.1 <- combined_count_out["ENST00000522937.1",]

enst_450.1 <- combined_count_out["ENST00000523450.1",]


fitteddata_2 <- rbind(fitteddata, enst_609.5, enst_937.1, enst_450.1)


fitteddata_with_rownames <- data.frame(RowNames = row.names(fitteddata_2), fitteddata_2, check.names = FALSE)



# AD vs AsymAD countdata

count_AD_Asym <- fitteddata_with_rownames[, c(1, 184:486)]


col_Data_2 <- col_Data_out[ 183:485, ]


# Filter data

filtered_data_1 <- count_AD_Asym[count_AD_Asym$RowNames %in% transcripts_of_interest, ]



# Reshape the data to long format for ggplot2
long_data_1 <- filtered_data_1 %>%
  pivot_longer(cols = -RowNames, names_to = "Sample", values_to = "Expression")

# Add group information from col_Data_2
col_Data_2 <- col_Data_2 %>%
  mutate(Sample = rownames(col_Data_2))

long_data_1 <- long_data_1 %>%
  left_join(col_Data_2, by = c("Sample"))



# Create the plot
ggplot(long_data_1, aes(x = group, y = Expression, color = group, shape = group, fill = group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.7, position = position_jitter(width = 0.2)) +
  theme_minimal() +
  labs(title = "Transcripts of ADAMTS2 Expressions in AD vs AsymAD (ROSMAP)",
       x = "Group", y = "Expression") +
  scale_color_manual(values = c("AD" = "blue", "AsymAD" = "purple")) +
  scale_fill_manual(values = c("AD" = "blue", "AsymAD" = "purple")) +
  geom_signif(comparisons = list(c("AD", "AsymAD")), 
              map_signif_level = TRUE, textsize = 3, test = "t.test", color = "black") +
  facet_wrap(~ RowNames, scales = "free") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +  # Ensure 4 breaks on y-axis
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 4), # Adjust Y-axis tick label size
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.major = element_blank(), # Remove large grids
    panel.grid.minor = element_blank(),  # Remove small grids
    axis.line = element_line(color = "black")  # Add black axis lines
  )



path <- "H:/plots_last/DET_plots"

ggsave(filename = file.path(path, "ADAMTS2_Expressions_AD_vs_AsymAD_ROSMAP_DET.tiff"),
       plot = last_plot(),
       width = 9, height = 7, device = "tiff", dpi = 600)





# Gene: MRPL1 ENSG00000169288


transcripts_of_interest <- c("ENST00000515625.1", "ENST00000511521.1", "ENST00000506674.1", "ENST00000504901.2", "ENST00000502384.3", "ENST00000315567.13")


x <- rownames(fitteddata)


transcripts_of_interest %in% x



enst_297.5 <- combined_count_out["ENST00000515625.1",]

enst_555.5 <- combined_count_out["ENST00000502384.3",]


fitteddata_2 <- rbind(fitteddata,enst_297.5)


fitteddata_2 <- rbind(fitteddata_2, enst_555.5)


fitteddata_with_rownames <- data.frame(RowNames = row.names(fitteddata_2), fitteddata_2, check.names = FALSE)




# AD vs AsymAD countdata

count_AD_Asym <- fitteddata_with_rownames[, c(1, 184:486)]


# Filter data

filtered_data_1 <- count_AD_Asym[count_AD_Asym$RowNames %in% transcripts_of_interest, ]


library(tidyverse)

# Reshape the data to long format for ggplot2
long_data_1 <- filtered_data_1 %>%
  pivot_longer(cols = -RowNames, names_to = "Sample", values_to = "Expression")

# Add group information from col_Data_2
col_Data_2 <- col_Data_2 %>%
  mutate(Sample = col_Data_2$sample_id)

long_data_1 <- long_data_1 %>%
  left_join(col_Data_2, by = c("Sample"))



# Create the plot
ggplot(long_data_1, aes(x = group, y = Expression, color = group, shape = group, fill = group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.7, position = position_jitter(width = 0.2)) +
  theme_minimal() +
  labs(title = "Transcripts of MRPL1 Expressions in AD vs AsymAD (ROSMAP)",
       x = "Group", y = "Expression") +
  scale_color_manual(values = c("AD" = "blue", "AsymAD" = "purple")) +
  scale_fill_manual(values = c("AD" = "blue", "AsymAD" = "purple")) +
  geom_signif(comparisons = list(c("AD", "AsymAD")), 
              map_signif_level = TRUE, textsize = 3, test = "t.test", color = "black") +
  facet_wrap(~ RowNames, scales = "free") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +  # Ensure 4 breaks on y-axis
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 4), # Adjust Y-axis tick label size
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black")  # Add black axis lines
  )



path <- "H:/plots_last/DET_plots"

ggsave(filename = file.path(path, "MRPL1_Expressions_AD_vs_AsymAD_ROSMAP_DET.tiff"),
       plot = last_plot(),
       width = 9, height = 7, device = "tiff", dpi = 600)




# Gene: PCOLCE ENSG00000106333


transcripts_of_interest <- c("ENST00000223061.6", "ENST00000496269.1", "ENST00000486440.5", "ENST00000472348.1", "ENST00000487172.1", "ENST00000462260.1", "ENST00000482863.1", "ENST00000468214.1", "ENST00000460002.1", "ENST00000490909.1")

x <- rownames(fitteddata)

transcripts_of_interest %in% x


enst_909.1 <- combined_count_out["ENST00000490909.1",]

enst_863.1 <- combined_count_out["ENST00000482863.1",]

enst_260.1 <- combined_count_out["ENST00000462260.1",]


fitteddata_2 <- rbind(fitteddata, enst_909.1 )



fitteddata_with_rownames <- data.frame(RowNames = row.names(fitteddata_2), fitteddata_2, check.names = FALSE)




# AD vs AsymAD countdata

count_AD_Asym <- fitteddata_with_rownames[, c(1, 184:486)]


col_Data_2 <- col_Data_out[ 183:485, ]


# Filter data

filtered_data_1 <- count_AD_Asym[count_AD_Asym$RowNames %in% transcripts_of_interest, ]


library(tidyverse)

# Reshape the data to long format for ggplot2
long_data_1 <- filtered_data_1 %>%
  pivot_longer(cols = -RowNames, names_to = "Sample", values_to = "Expression")

# Add group information from col_Data_2
col_Data_2 <- col_Data_2 %>%
  mutate(Sample = rownames(col_Data_2))

long_data_1 <- long_data_1 %>%
  left_join(col_Data_2, by = c("Sample"))




# Create the plot
ggplot(long_data_1, aes(x = group, y = Expression, color = group, shape = group, fill = group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.7, position = position_jitter(width = 0.2)) +
  theme_minimal() +
  labs(title = "Transcripts of PCOLCE Expressions in AD vs AsymAD (ROSMAP)",
       x = "Group", y = "Expression") +
  scale_color_manual(values = c("AD" = "blue", "AsymAD" = "purple")) +
  scale_fill_manual(values = c("AD" = "blue", "AsymAD" = "purple")) +
  geom_signif(comparisons = list(c("AD", "AsymAD")), 
              map_signif_level = TRUE, textsize = 3, test = "t.test", color = "black") +
  facet_wrap(~ RowNames, scales = "free") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +  # Ensure 4 breaks on y-axis
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 4), # Adjust Y-axis tick label size
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black")  # Add black axis lines
  )



path <- "H:/plots_last/DET_plots"

ggsave(filename = file.path(path, "PCOLCE_Expressions_AD_vs_AsymAD_ROSMAP_DET.tiff"),
       plot = last_plot(),
       width = 9, height = 7, device = "tiff", dpi = 600)





