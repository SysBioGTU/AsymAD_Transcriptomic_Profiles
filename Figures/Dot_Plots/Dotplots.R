
library(ggplot2)
library(dplyr)

# Create the data frame according to the chosen terms from ROSMAP enrichment

df <- data.frame(
  term = c("regulation of GTPase activity", "neurotransmitter secretion", "synaptic vesicle cycle",
           "synapse organisation", "cognition", "glutamatergic synapse",
           "postsynaptic specialisation", "GABAergic synapse", "transmission across chemical synapses",
           "glutamatergic synapse", "presynaptic membrane", "postsynaptic density",
           "glial cell differentiation", "regulation of synaptic plasticity"),
  log10_p_value = c(3.05868433925078, 3.7161123436789, 1.43547333808953,
                    2.3382032498538, 1.99454285665662, 1.59268809101131,
                    5.52021215546944, 2.23770698265947, 2.0169788462494,
                    1.59268809101131, 1.47149649078908, 3.79435080507582,
                    2.09714718211474, 2.30812637688591),
  Count = c(81, 29, 34, 60, 44, 63, 53, 19, 41, 21, 25, 46, 38, 34)
)

# Create the dot plot

ggplot(df, aes(x = log10_p_value, y = reorder(term, log10_p_value), size = Count, color = log10_p_value)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "red", high = "green", name = "log10(p-value)") +  # Renk skalası
  scale_size(range = c(2, 10), name = "Count") +  # Count için boyut skalası
  theme_minimal() +
  labs(x = "log10(p-value)", y = "", title = "Dot Plot of Enriched Terms") +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"))

path <- "~/Desktop/Figürler"

ggsave(filename = file.path(path, "ROSMAP_all_enrichment.tiff"),
       plot = last_plot(),
       width = 9, height = 6, device='tiff', dpi=600)




# Create the data frame according to the chosen terms from MSBB enrichment

df <- data.frame(
  term = c("ATP binding", "vesicle-mediated transport in synapse", "synapse",
           "microtubule cytoskeleton", "glutamatergic synapse", "heme biosynthesis",
           "lipid metabolic process", "mitochondrial inner membrane",
           "regulation of postsynapse organisation", "cellular response to stress"),
  log10_p_value = c(3.00313163545272, 1.41508876249324, 4.83288946451113,
                    4.77482380028829, 1.63565722043885, 1.90971210639947,
                    2.43736799748681, 2.56166993308191, 1.7797292343038,
                    2.32289652363223),
  Count = c(206, 41, 197, 190, 64, 8, 179, 77, 35, 222)
)

# Create the dot plot

ggplot(df, aes(x = log10_p_value, y = reorder(term, log10_p_value), size = Count, color = log10_p_value)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "red", high = "green", name = "log10(p-value)") +  # Renk skalası
  scale_size(range = c(2, 10), name = "Count") +  # Count için boyut skalası
  theme_minimal() +
  labs(x = "log10(p-value)", y = "", title = "Dot Plot of Enriched Terms") +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"))


path <- "~/Desktop/Figürler"

ggsave(filename = file.path(path, "MSBB_all_enrichment.tiff"),
       plot = last_plot(),
       width = 9, height = 6, device='tiff', dpi=600)



