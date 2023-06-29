require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

dataset_type <- "archaea"
scores_sm <- read_csv(paste(paste("../results/summary_scores_methods", dataset_type, sep = "-"), "csv", sep = "."))
scores_sm <- scores_sm %>%
  filter(Taxonomic_Rank != "superkingdom") %>%
  mutate(Distance_to_closest = factor(Distance_to_closest, levels = unique(Distance_to_closest)))

scores_sm$Taxonomic_Rank <- factor(
  scores_sm$Taxonomic_Rank,
  levels = c("phylum", "class", "order", "family", "genus", "species")
)

ggplot(scores_sm, aes(y = Precision, x = Recall, color = Distance_to_closest, shape = Method)) +
  facet_wrap(vars(Taxonomic_Rank)) +
  geom_path(aes(group = Distance_to_closest), color = "darkgray") +
  geom_point(aes(group = Method), size = 2.5, alpha = 0.85) +
  labs(shape = "Tool", colour = "Distance to closest", y = "Precision", x = "Recall") +
  theme_cowplot(font_size = 16) +
  scale_colour_brewer(palette = "Paired") +
  theme(aspect.ratio = 1.1, panel.spacing.x = unit(1.15, "lines"))
ggsave2(paste(paste("../figures/classification_tool_comparison", dataset_type, sep = "-"), "1.pdf", sep="-"), width=10, height = 6)

ggplot(scores_sm, aes(y = Distance_to_closest, x = F1, color = Method, shape = Method)) +
  facet_wrap(vars(Taxonomic_Rank)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(shape = "Tool", colour = "Tool", y = "Distance to the closest", x = "F1") +
  theme_minimal_vgrid(font_size = 17) +
  theme(strip.background = element_rect(fill = "gray")) +
  scale_colour_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), aspect.ratio = 0.9) +
  theme(panel.spacing.x = unit(1, "lines"))
ggsave2(paste(paste("../figures/classification_tool_comparison", dataset_type, sep = "-"), "2.pdf", sep="-"), width=10, height = 5)