require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

pm <- read_tsv("../results/version_comparison-profiling_evaluation-CAMI1_hc.tsv")
pm <- pm %>% filter(metric == "Bray-Curtis distance" | metric == "Shannon equitability")
pm <- pm %>% filter(rank != "superkingdom") %>% filter(rank != "strain")
pm <- pm %>%
  group_by(metric) %>%
  mutate(dvalue = abs(value - value[tool == "Gold standard"]))
pm <- pm %>% filter(tool != "Gold standard")
pm$rank <- factor(
  pm$rank,
  levels = c("phylum", "class", "order", "family", "genus", "species")
)
pm <- pm %>%
  group_by(rank, tool, metric) %>%
  summarise(lower = min(dvalue), upper = max(dvalue), p = mean(dvalue))

pm[pm == "CONSULT-II (new - all ranks)"] <- "With read-level normalization"
pm[pm == "CONSULT-II (old - all ranks)"] <- "Without read-level normalization"
pm[pm == "CONSULT-II (old - based on species, RECOMB-CG)"] <- "With read-level normalization (based on species profile)"
pm[pm == "CONSULT-II (new - based on species)"] <- "Without read-level normalization (based on species profile)"

pm$tool <- factor(
  pm$tool,
  levels = c("Without read-level normalization", "With read-level normalization", "With read-level normalization (based on species profile)", "Without read-level normalization (based on species profile)")
)

pm_all <- pm %>% filter(tool == "Without read-level normalization" | tool == "With read-level normalization")
fs <- 17
ggplot(pm_all, mapping = aes(x = p, y = rank, shape = tool, color = tool)) +
  # geom_point(size = 3, alpha = 0.9) +
  geom_pointrange(size = 1, linewidth = 1, mapping = aes(xmin = lower, xmax = upper)) +
  xlim(0, NA) +
  facet_wrap(vars(metric), ncol = 1, scales = "free_x") +
  theme_minimal_grid(font_size = fs) +
  labs(shape = "Approach", color = "Approach", y = "", x = "Absolute difference to gold standard") +
  scale_colour_brewer(palette = "Set1") +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(legend.direction = "vertical", legend.position = "bottom", legend.justification = "center", legend.text = element_text(size = fs)) +
  theme(axis.title.y =  element_text(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), aspect.ratio = 1.2)
ggsave2("../figures/profiling_normalization_comparison-all_ranks.pdf", width = 5, height = 10)

ggplot(pm, mapping = aes(x = p, y = rank, shape = tool, color = tool)) +
  # geom_point(size = 3, alpha = 0.9) +
  geom_pointrange(size = 1, linewidth = 1, mapping = aes(xmin = lower, xmax = upper)) +
  xlim(0, NA) +
  facet_wrap(vars(metric), ncol = 1, scales = "free_x") +
  theme_minimal_grid(font_size = fs) +
  labs(shape = "Approach", color = "Approach", y = "", x = "Absolute difference to gold standard") +
  scale_colour_brewer(palette = "Set1") +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(legend.direction = "vertical", legend.position = "bottom", legend.justification = "center", legend.text = element_text(size = fs)) +
  theme(axis.title.y =  element_text(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), aspect.ratio = 1.2)
ggsave2("../figures/profiling_normalization_comparison-from_species.pdf", width = 5, height = 10)