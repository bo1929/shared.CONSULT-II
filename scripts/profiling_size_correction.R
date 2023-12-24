require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
library(tidyr)

pm <- read_tsv("../results/results_profiling-marine-vt1000-norm_genome_sizes-n_f1/results.tsv") %>% 
  filter(tool %in% c("CONSULT-II v0.4.0", "Gold standard")  & rank != "strain") %>% 
  group_by(sample, metric) %>% 
  mutate(dvalue = abs(value - value[tool == "Gold standard"])) %>% 
  filter(tool != "Gold standard") %>%
  group_by(rank, tool, metric) %>%
  summarise(lower = min(dvalue), upper = max(dvalue), p = mean(dvalue))
pm$rank <- factor(
  pm$rank,
  levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
)
pm$metric <- as.factor(pm$metric)
pm$tool[pm$tool == "CONSULT-II v0.4.0"] <- "with genome size correction"

xm <- read_tsv("../results/results_profiling-marine-vt1000-norm_free_baseline-n_f1/results.tsv") %>% 
  filter(tool %in% c("CONSULT-II v0.4.0", "Gold standard")  & rank != "strain") %>% 
  group_by(sample, metric) %>% 
  mutate(dvalue = abs(value - value[tool == "Gold standard"])) %>% 
  filter(tool != "Gold standard") %>%
  group_by(rank, tool, metric) %>%
  summarise(lower = min(dvalue), upper = max(dvalue), p = mean(dvalue))
xm$rank <- factor(
  xm$rank,
  levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
)
xm$metric <- as.factor(xm$metric)
xm$tool[xm$tool == "CONSULT-II v0.4.0"] <- "without genome size correction"

tm <- pivot_wider(
  rbind(xm,pm),
  id_cols = c(rank, tool),
  names_from = metric,
  values_from = p
)
tm <- tm[order(tm$rank),]
tm$`Unweighted UniFrac (CAMI) (unfiltered)` <- rep((tm %>% filter(is.na(rank)))$`Unweighted UniFrac (CAMI) (unfiltered)`, 8)
tm$`Weighted UniFrac (CAMI) (unfiltered)` <- rep((tm %>% filter(is.na(rank)))$`Weighted UniFrac (CAMI) (unfiltered)`, 8)
p1 <- ggplot(tm %>% filter(!rank %in% c(NA, "superkingdom")),
             aes(x = 2-`L1 norm error (unfiltered)`,
                 y = 16-`Weighted UniFrac (CAMI) (unfiltered)`,
                 color = tool)
) +
  facet_wrap("rank", scales = "fixed") +
  geom_point(size=3, alpha=0.7) +
  labs(x="2 - L1 norm error", y="16 - weighted UniFrac error", shape="", color="Tool") +
  ylim(5, 16) +
  xlim(1, 2) +
  scale_color_brewer(palette = "Set1") +
  scale_shape_discrete(guide = "none") +
  theme_cowplot(font_size = 17)
p1

pm <- read_tsv("../results/results_profiling-strain_madness-vt1000-norm_genome_sizes-n_f1/results.tsv") %>% 
  filter(tool %in% c("CONSULT-II v0.4.0", "Gold standard")  & rank != "strain") %>% 
  group_by(sample, metric) %>% 
  mutate(dvalue = abs(value - value[tool == "Gold standard"])) %>% 
  filter(tool != "Gold standard") %>%
  group_by(rank, tool, metric) %>%
  summarise(lower = min(dvalue), upper = max(dvalue), p = mean(dvalue))
pm$rank <- factor(
  pm$rank,
  levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
)
pm$metric <- as.factor(pm$metric)
pm$tool[pm$tool == "CONSULT-II v0.4.0"] <- "with genome size correction"

xm <- read_tsv("../results/results_profiling-strain_madness-vt1000-norm_free_baseline-n_f1/results.tsv") %>% 
  filter(tool %in% c("CONSULT-II v0.4.0", "Gold standard")  & rank != "strain") %>% 
  group_by(sample, metric) %>% 
  mutate(dvalue = abs(value - value[tool == "Gold standard"])) %>% 
  filter(tool != "Gold standard") %>%
  group_by(rank, tool, metric) %>%
  summarise(lower = min(dvalue), upper = max(dvalue), p = mean(dvalue))
xm$rank <- factor(
  xm$rank,
  levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
)
xm$metric <- as.factor(xm$metric)
xm$tool[xm$tool == "CONSULT-II v0.4.0"] <- "without genome size correction"
pm <- rbind(xm,pm)

tm <- pivot_wider(
  pm,
  id_cols = c(rank, tool),
  names_from = metric,
  values_from = p
)
tm <- tm[order(tm$rank),]
tm$`Unweighted UniFrac (CAMI) (unfiltered)` <- rep((tm %>% filter(is.na(rank)))$`Unweighted UniFrac (CAMI) (unfiltered)`, 8)
tm$`Weighted UniFrac (CAMI) (unfiltered)` <- rep((tm %>% filter(is.na(rank)))$`Weighted UniFrac (CAMI) (unfiltered)`, 8)
p2 <- ggplot(tm %>% filter(!rank %in% c(NA, "superkingdom")),
             aes(x = 2-`L1 norm error (unfiltered)`,
                 y = 16-`Weighted UniFrac (CAMI) (unfiltered)`,
                 color = tool)
) +
  facet_wrap("rank", scales = "fixed") +
  geom_point(size=3, alpha=0.7) +
  labs(x="2 - L1 norm error", y="16 - weighted UniFrac error", shape="", color="Tool") +
  ylim(5, 16) +
  xlim(1, 2) +
  scale_color_brewer(palette = "Set1") +
  scale_shape_discrete(guide = "none") +
  theme_cowplot(font_size = 17)
p2

plot_grid(
  p1 + ggtitle("Marine dataset") + theme(legend.box.margin = margin(0, 0, 0, 0), legend.position = "bottom"),
  p2 + ggtitle("Strain-madness dataset") + theme(legend.position = "none"),
  ncol = 1,
  rel_heights = c(1.2, 1)
)
ggsave2("../figures/profiling-size_correction.pdf", width = 9, height = 8)