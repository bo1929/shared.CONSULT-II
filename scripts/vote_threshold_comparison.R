require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

scores_cc <- read_csv("../results/summary_scores_c.csv")
scores_cc <- scores_cc %>%
  filter(Taxonomic_Rank != "superkingdom") %>%
  mutate(Distance_to_closest = factor(Distance_to_closest, levels = unique(Distance_to_closest)))

scores_cc$Taxonomic_Rank <- factor(
  scores_cc$Taxonomic_Rank,
  levels = c("phylum", "class", "order", "family", "genus", "species")
)

ggplot(scores_cc, aes(x = Precision, y = Recall, color = Distance_to_closest, shape = Method, group = Method)) +
  facet_wrap(vars(Taxonomic_Rank)) +
  geom_line(aes(group = Distance_to_closest), color = "gray") +
  geom_point(size = 3, alpha = 0.8) +
  labs(shape = "Vote threshold", colour = "Distance to closest", x = "Precision", y = "Recall") +
  theme_cowplot(font_size = 16.5) +
  scale_colour_brewer(palette = "Paired") +
  theme(aspect.ratio = 1, panel.spacing.x = unit(1.25, "lines"))
ggsave2("../figures/vote_threshold_comparison-2.pdf")

ggplot(scores_cc, aes(x = Distance_to_closest, y = F1)) +
  facet_wrap(vars(Taxonomic_Rank)) +
  geom_line(aes(group = Method, linetype = Method), color="#d95f02") +
  geom_point(size = 1.75, alpha = 0.85, shape="triangle", color="#d95f02") +
  labs(shape = "Vote threshold", linetype = "Vote threshold", x = "Distance to closest", y = "F1") +
  theme_cowplot(font_size = 17) +
  scale_colour_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), aspect.ratio = 0.9) +
  theme(panel.spacing.x = unit(1, "lines"))
ggsave2("../figures/vote_threshold_comparison-1.pdf")