require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

dataset_type <- "archaea"
scores_sm <- read_csv(paste(paste("../results/summary_scores_methods", dataset_type, sep = "-"), "csv", sep = "."))
scores_spm <- read_csv(paste(paste("../results/summary_p_scores_methods", dataset_type, sep = "-"), "csv", sep = "."))

scores_sm <- scores_sm %>%
  filter(Taxonomic_Rank != "superkingdom") %>%
  mutate(Distance_to_closest = factor(Distance_to_closest, levels = unique(Distance_to_closest)))
scores_spm <- scores_spm %>%
  filter(Taxonomic_Rank != "superkingdom")

scores_sm$Taxonomic_Rank <- factor(
  scores_sm$Taxonomic_Rank,
  levels = c("phylum", "class", "order", "family", "genus", "species")
)
scores_spm$Taxonomic_Rank <- factor(
  scores_spm$Taxonomic_Rank,
  levels = c("phylum", "class", "order", "family", "genus", "species")
)
scores_spm <- scores_spm %>% filter(Distance_to_closest < 0.9)

ggplot(scores_spm, aes(x = Distance_to_closest, y = F1, color = Method, linetype = Method, shape = Method)) +
  facet_wrap(vars(Taxonomic_Rank), nrow = 1) +
  geom_point(alpha = 0.5) +
  labs(shape = "Tool", colour = "Tool", linetype = "Tool", x = "Distance to the closest", y = "F1") +
  stat_smooth(se = F, span = 0.7, method = "glm", method.args = list(family = binomial), size = 1) +
  scale_colour_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("22","solid", "42")) +
  theme_cowplot(font_size = 17) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), aspect.ratio = 1.25) +
  # theme(legend.position = "bottom", legend.justification = "center") +
  theme(panel.spacing.x = unit(1, "lines"))
ggsave2(paste(paste("../figures/classification_tool_comparison", dataset_type, sep = "-"), "2.pdf", sep="-"), width=15, height = 3.75)

ggplot(scores_spm %>%
         mutate(
           Distance_to_closest = cut(
             Distance_to_closest,
             include.lowest = TRUE,
             breaks = c(0.00, 0.001, 0.02, 0.05, 0.1, 0.125, 0.20, 0.40, 1.0)
           )
         ) %>%
         group_by(Method, Taxonomic_Rank, Distance_to_closest) %>%
         summarise(Recall = mean(Recall), F1 = mean(F1), Precision = mean(Precision))
       , aes(x = Precision, y = Recall, color = Distance_to_closest, shape = Method)) +
  facet_wrap(vars(Taxonomic_Rank), nrow = 1) +
  geom_path(aes(group = Distance_to_closest), color = "darkgray") +
  geom_point(aes(group = Method), size = 2.5, alpha = 0.85) +
  labs(shape = "Tool", colour = "Distance to the closest", x = "Precision", y = "Recall") +
  scale_colour_brewer(palette = "Paired") +
  theme_cowplot(font_size = 16) +
  # theme(legend.position = "bottom", legend.justification = "center") +
  theme(legend.text = element_text(size=12.75), legend.title = element_text(size=14)) +
  theme(axis.text.y = element_text(size=12.25), axis.text.x = element_text(size=12.25)) +
  theme(aspect.ratio = 1.2, panel.spacing.x = unit(1.15, "lines"))
ggsave2(paste(paste("../figures/classification_tool_comparison", dataset_type, sep = "-"), "1.pdf", sep = "-"), width=15, height = 3.7)

ggplot(scores_spm %>%
         mutate(
           Distance_to_closest = cut(
             Distance_to_closest,
             include.lowest = TRUE,
             breaks = c(0.00, 0.001, 0.02, 0.05, 0.1, 0.125, 0.20, 0.40, 1.0)
           )
         ) %>%
         group_by(Method, Taxonomic_Rank, Distance_to_closest) %>%
         summarise(Recall = mean(Recall), F1 = mean(F1), Precision = mean(Precision)),
       aes(x = Distance_to_closest, y = F1, color = Method, shape = Method)) +
  facet_wrap(vars(Taxonomic_Rank), nrow = 1) +
  geom_point(size = 3, alpha = 0.85) +
  geom_line(aes(group=Method)) +
  labs(shape = "Tool", colour = "Tool", x = "Distance to the closest", y = "F1") +
  scale_colour_brewer(palette = "Dark2") +
  theme_cowplot(font_size = 17) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust = 1), aspect.ratio = 1.25) +
  # theme(legend.position = "bottom", legend.justification = "center") +
  theme(panel.spacing.x = unit(1, "lines"))
ggsave2(paste(paste("../figures/classification_tool_comparison", dataset_type, sep = "-"), "3.pdf", sep="-"), width=13.5, height = 3.75)