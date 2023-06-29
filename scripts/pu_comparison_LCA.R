require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

summary_scores_pu <- read_csv("../results/summary_scores_pu.csv")
summary_scores_pu <- summary_scores_pu %>%
  filter(Taxonomic_Rank != "superkingdom") %>%
  mutate(Distance_to_closest = factor(Distance_to_closest, levels = unique(Distance_to_closest)))

summary_scores_pu$Taxonomic_Rank <- factor(
  summary_scores_pu$Taxonomic_Rank,
  levels = c("phylum", "class", "order", "family", "genus", "species")
)

ggplot(summary_scores_pu, aes(x = Distance_to_closest, y = F1, color = Method, shape = Method, group = Method)) +
  facet_wrap(vars(Taxonomic_Rank)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line(aes(linetype = Method)) +
  labs(linetype = "LCA", shape = "LCA", colour = "LCA", x = "Distance to the closest", y = "F1") +
  theme_cowplot(font_size = 15) +
  scale_colour_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), aspect.ratio = 0.8)
ggsave2("../figures/pu_comparison_LCA.pdf", width = 8, height = 5)