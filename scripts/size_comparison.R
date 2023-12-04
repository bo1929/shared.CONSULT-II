require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

scores <- read_csv("../results/summary_p_scores_sizes-bacteria.csv")
scores$Taxonomic_Rank[scores$Taxonomic_Rank == "superkingdom"] <- "kingdom"
scores <- scores %>% filter (Taxonomic_Rank != "kingdom")
scores$Taxonomic_Rank <- factor(
  scores$Taxonomic_Rank,
  levels = c("phylum", "class", "order", "family", "genus", "species")
)
scores$Method[scores$Method == "CONSULT-II (0.03)"] <- "CONSULT-II 140Gb (0.03)"
scores$Method[scores$Method == "CONSULT-II h=14, b=10 (0.03)"] <- "CONSULT-II 32Gb (0.03)"
scores$Method[scores$Method == "CONSULT-II h=13, b=16 (0.03)"] <- "CONSULT-II 18Gb (0.03)"
scores$Method[scores$Method == "Kraken-II"] <- "Kraken-II 44Gb"

ggplot(scores %>% filter(Method != "CONSULT-II 18Gb (0.03)") %>%
         mutate(
           Distance_to_closest = cut(
             Distance_to_closest,
             include.lowest = TRUE,
             breaks = c(0, 0.001, 0.02, 0.06, 0.12, 0.16, 0.35)
           )
         ) %>%
         group_by(Method, Taxonomic_Rank, Distance_to_closest) %>%
         summarise(Recall = mean(Recall), F1 = mean(F1), Precision = mean(Precision))
       , aes(x = Precision, y = Recall, color = Distance_to_closest, shape = Method)) +
  geom_point(alpha = 0.8, size=3) +
  facet_wrap(facets = "Taxonomic_Rank", nrow=1) +
  labs(shape = "Tool", colour = "Tool", linetype = "Tool", x = "Precision", y = "Recall") +
  scale_colour_brewer(palette = "Paired") +
  scale_shape_manual(values = c(17, 2, 15)) +
  theme_cowplot(font_size = 17) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), aspect.ratio = 1.25) +
  theme(panel.spacing.x = unit(1, "lines")) +
  theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.box = "vertical")
ggsave2("../figures/classification_size_comparison-bacteria-1.pdf", width = 15, height = 5)

ggplot(scores %>% filter(Method != "CONSULT-II 18Gb (0.03)"), aes(x = Distance_to_closest, y = F1, color = Method, linetype = Method, shape = Method)) +
  facet_wrap(vars(Taxonomic_Rank), nrow = 1) +
  geom_point(alpha = 0.5) +
  labs(shape = "Tool", colour = "Tool", linetype = "Tool", x = "Distance to the closest", y = "F1") +
  stat_smooth(se = F, span = 0.7, method = "glm", method.args = list(family = binomial), size = 1.25) +
  theme_cowplot(font_size = 17) +
  scale_colour_manual(values = c("#d95f02", "#a6761d", "#7570b3")) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), aspect.ratio = 1.25) +
  theme(legend.position = "bottom", legend.justification = "center", legend.direction = "vertical") +
  theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.box = "vertical")
ggsave2("../figures/classification_size_comparison-bacteria-2.pdf", width = 15, height = 4)
