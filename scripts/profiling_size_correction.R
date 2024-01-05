require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
library(tidyr)

pm <- read_tsv("../results/profile-filtered_tv1000-norm_genome_sizes-marine/results.tsv") %>% 
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

xm <- read_tsv("../results/profile-filtered_tv1000-norm_free_baseline-marine/results.tsv") %>% 
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
  labs(x="2 - L1 norm error", y="16 - weighted UniFrac error", shape="", color="Approach") +
  ylim(5, 16) +
  xlim(1, 2) +
  scale_color_brewer(palette = "Set1") +
  scale_shape_discrete(guide = "none") +
  theme_cowplot(font_size = 17)
p1

pm <- read_tsv("../results/profile-filtered_tv1000-norm_genome_sizes-strain_madness/results.tsv") %>% 
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

xm <- read_tsv("../results/profile-filtered_tv1000-norm_free_baseline-strain_madness/results.tsv") %>% 
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
  labs(x="2 - L1 norm error", y="16 - weighted UniFrac error", shape="", color="Approach") +
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
ggsave2("../figures/profiling-size_correction-CAMI2.pdf", width = 9, height = 8)

pm <- read_tsv("../results/all_tools-profiling_evaluation-CAMI1_hc.tsv")
pm <- pm %>% filter(rank != "superkingdom") %>% filter(rank != "strain")
pm <- pm %>%
  group_by(sample, metric) %>%
  mutate(dvalue = abs(value - value[tool == "Gold standard"]))
pm <- pm %>% filter(tool != "Gold standard")
pm$rank <- factor(
  pm$rank,
  levels = c("phylum", "class", "order", "family", "genus", "species")
)
pm <- pm %>% filter(tool %in% c("CONSULT-II (Eq. 7)", "CONSULT-II (Eq. 6)"))
pm$tool[pm$tool == "CONSULT-II (Eq. 7)"] <- "with genome size correction"
pm$tool[pm$tool == "CONSULT-II (Eq. 6)"] <- "without genome size correction"

pm$tool <- factor(
  pm$tool,
  levels = c("with genome size correction", "without genome size correction")
)
pm <- pm %>%
  group_by(rank, tool, metric) %>%
  summarise(lower = min(dvalue), upper = max(dvalue), p = mean(dvalue))
pm$metric <- as.factor(pm$metric)

bc_plot <- ggplot(pm %>% filter(metric == "Bray-Curtis distance"), mapping = aes(x = p, y = tool, color = tool)) +
  geom_pointrange(size = 1, alpha = 0.85, linewidth = 1, mapping = aes(xmin = lower, xmax = upper)) +
  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 15) +
  labs(color = "Approach", y = "", x = "Bray-Curtis dissimilarity to true profile") +
  scale_color_brewer(palette = "Set1") +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 1, size=9)) +
  theme(aspect.ratio = 0.4, panel.spacing.x = unit(1.5, "lines")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
bc_plot
se_plot <- ggplot(pm %>% filter(metric == "Shannon equitability"), mapping = aes(x = p, y = tool, color = tool)) +
  geom_pointrange(size = 1, alpha = 0.85, linewidth = 1, mapping = aes(xmin = lower, xmax = upper)) +
  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 15) +
  labs(color = "Approach", y = "", x = "Shannon's equitability (absolute difference to gold standard)") +
  scale_color_brewer(palette = "Set1") +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 1, size=9)) +
  theme(aspect.ratio = 0.4, panel.spacing.x = unit(1.5, "lines")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))
se_plot
prow <- plot_grid(
  bc_plot + theme(legend.position = "none"),
  se_plot + theme(legend.position = "none"),
  ncol = 2,
  rel_heights = c(1, 1),
  vjust = 5,
  hjust = 5
)
legend <- get_legend(
  # create some space to the left of the legend
  se_plot + theme(legend.box.margin = margin(-1, 0, 0, 0), legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")
)
plot_grid(legend, prow, rel_heights = c(.15, 3), ncol = 1)

ggsave2("../figures/profiling-size_correction-CAMI1.pdf", width = 14, height = 3.5)