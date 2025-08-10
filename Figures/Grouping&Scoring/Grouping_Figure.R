
library(ggplot2)

data_2 <- Book2


color_palette <- c(
  "R3-M1" = "#FF7F7F",    
  "R3-M2" = "grey",    
  "R3-M3" = "grey",    
  "R2-M1" = "#9370DB",    
  "R2-M2" = "#647CB0",
  "R2-M3" = "darkgrey", 
  "R1-M1" = "#20B2AA",
  "R1-M2" = "#ADD8E6",
  "R1-M3" = "grey"
)


# Creating a bar plot

ggplot(data_2, aes(x = Scores, y = Genes, fill = Scores)) +
  geom_bar(stat = "identity", color = "black") +  # Defining the bars
  geom_text(aes(label = Genes), vjust = -0.5, size = 5) +  # Adding the numbers
  scale_fill_manual(values = color_palette) +  # Color assignment
  theme_minimal() +
  labs(
    title = "Number of Shared Genes for Each Scored Group",
    x = "Dataset Scores"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14)
  )

path <- "/Volumes/Untitled/Tubitak_Proje_Revizyon"

ggsave(filename = file.path(path, "scoring_plot.tiff"),
       plot = last_plot(),
       width = 10, height = 6, device='tiff', dpi=600)
