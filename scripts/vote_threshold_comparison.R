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
  facet_wrap(vars(Taxonomic_Rank), nrow = 1) +
  geom_path(aes(group = Distance_to_closest), color = "darkgray") +
  geom_point(aes(group = Method), size = 2.5, alpha = 0.85) +
  labs(shape = "Vote threshold", colour = "MinGND", x = "Precision", y = "Recall") +
  scale_colour_brewer(palette = "Paired") +
  theme_cowplot(font_size = 16) +
  theme(legend.text = element_text(size = 12.75), legend.title = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 12.25), axis.text.x = element_text(size = 12.25)) +
  theme(aspect.ratio = 1.2, panel.spacing.x = unit(1.15, "lines"))
ggsave2("../figures/vote_threshold_comparison-2.pdf", , width = 15, height = 3.7)

ggplot(scores_cc, aes(x = Distance_to_closest, y = F1, shape = Method)) +
  facet_wrap(vars(Taxonomic_Rank), nrow = 1) +
  geom_line(aes(group = Method, linetype = Method, color = Method)) +
  geom_point(size = 1.75, alpha = 0.85, aes(shape = Method, color = Method)) +
  labs(color = "Vote threshold", shape = "Vote threshold", linetype = "Vote threshold", x = "Distance to the closest (MinGND)", y = "F1") +
  scale_colour_brewer(palette = "Set1") +
  scale_linetype_manual(values = c("22", "32", "42")) +
  theme_cowplot(font_size = 17) +
  theme(legend.text = element_text(size = 12.75), legend.title = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 12.25), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12.25)) +
  theme(aspect.ratio = 1.2, panel.spacing.x = unit(1.15, "lines"))
ggsave2("../figures/vote_threshold_comparison-1.pdf", width = 12, height = 3.75)
