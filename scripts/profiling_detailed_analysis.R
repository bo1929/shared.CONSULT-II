require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
library(tidyr)

pm <- read_tsv("../results/profile-filtered_tv1000-norm_genome_sizes-strain_madness/results.tsv")
# pm <- read_tsv("../results/profile-filtered_tv1000-norm_genome_sizes-marine/results.tsv")
pm <- pm %>% filter(rank != "strain")
pm <- pm %>% group_by(sample, metric) %>% mutate(dvalue = abs(value - value[tool == "Gold standard"]))
pm <- pm %>% filter(tool != "Gold standard")
pm <- pm %>%
  group_by(rank, tool, metric) %>%
  summarise(lower = min(dvalue), upper = max(dvalue), p = mean(dvalue))
pm  <- pm %>% filter(
  tool %in% c(
    "Bracken 2.2",
    "CCMetagen 1.1.3",
    "Centrifuge 1.0.4 beta",
    "CONSULT-II v0.4.0",
    "DUDes 0.08",
    "FOCUS 1.5",
    "Metalign 0.6.2",
    "MetaPhlAn 2.9.22",
    "MetaPhyler 1.25",
    "mOTUs 2.5.1_2",
    "MetaPalette 1.0.0",
    "TIPP 4.3.10",
    "NBC++"
  )
)
pm$tool[pm$tool == "mOTUs 2.5.1_2"] <- "mOTUs 2.5.1"
pm$rank <- factor(
  pm$rank,
  levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
)
pm$metric <- as.factor(pm$metric)
pm$tool <- factor(pm$tool, levels = c(
    "CONSULT-II v0.4.0",
    "Bracken 2.2",
    "CCMetagen 1.1.3",
    "Centrifuge 1.0.4 beta",
    "DUDes 0.08",
    "FOCUS 1.5",
    "Metalign 0.6.2",
    "MetaPhlAn 2.9.22",
    "MetaPhyler 1.25",
    "mOTUs 2.5.1",
    "TIPP 4.3.10",
    "NBC++",
    "MetaPalette 1.0.0"
  )
)

bc_plot <- ggplot(pm %>% filter(metric == "Bray-Curtis distance" & rank != "superkingdom"),
                  mapping = aes(x = p, y = reorder(tool, p), color = tool, shape=grepl("CONS",tool))) +
  geom_pointrange(size = 0.25, linewidth = 0.5, 
                  mapping = aes(xmin = lower, xmax = upper)) +
  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 17) +
  labs(color = "Tool", y = "", x = "Bray-Curtis dissimilarity to true profile") +
  scale_colour_manual(
    values=c("black", RColorBrewer::brewer.pal(12,"Paired"))
  ) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 1, size=9)) +
  theme(aspect.ratio = 0.4, panel.spacing.x = unit(1.5, "lines")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_shape_manual(values=c(1,16), guide = "none")
bc_plot
se_plot <- ggplot(pm %>% filter(metric == "Shannon equitability" & rank != "superkingdom"), mapping = aes(x = p, y = tool, color = tool, shape=grepl("CONS",tool))) +
  geom_pointrange(size = 0.25, linewidth = 0.5, 
                  mapping = aes(xmin = lower, xmax = upper)) +  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 17) +
  labs(color = "Tool", y = "", x = "Shannon's equitability (absolute difference to gold standard)") +
  scale_colour_manual(
    values=c("black", RColorBrewer::brewer.pal(12,"Paired"))
  ) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 1, size=9)) +
  theme(aspect.ratio = 0.4, panel.spacing.x = unit(1.5, "lines")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_shape_manual(values=c(1,16), guide = "none")
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
  se_plot + theme(legend.box.margin = margin(0, 0, 0, 12))
)
plot_grid(prow, legend, rel_widths = c(3, .45))

tm <- pivot_wider(
    pm,
    id_cols = c(rank, tool),
    names_from = metric,
    values_from = p
  )
p1 <- ggplot(tm %>% filter(rank != "na" & rank != "superkingdom"),
       aes(y = 1-`Completeness (unfiltered)`,
           x = 1-`Purity (unfiltered)`, 
           color = tool,
           shape = grepl("CONS", tool))
       ) +
  geom_point(size=3, alpha=0.8) +
  facet_wrap("rank") +
  xlim(0, NA) +
  labs(x="Purity", y="Completeness", shape="", color="Tool") +
  scale_shape_discrete(guide = "none") +
  theme_cowplot(font_size = 17) +
  scale_colour_manual(
    values=c("black", RColorBrewer::brewer.pal(12,"Paired"))
  ) +
  theme(panel.spacing.x = unit(1.5, "lines"))
p1
tm$`Unweighted UniFrac (CAMI) (unfiltered)` <- rep((tm %>% filter(is.na(rank)))$`Unweighted UniFrac (CAMI) (unfiltered)`, 8)
tm$`Weighted UniFrac (CAMI) (unfiltered)` <- rep((tm %>% filter(is.na(rank)))$`Weighted UniFrac (CAMI) (unfiltered)`, 8)
p2 <- ggplot(tm %>% filter(!rank %in% c(NA, "superkingdom")),
       aes(x = 2-`L1 norm error (unfiltered)`,
           y = 16-`Weighted UniFrac (CAMI) (unfiltered)`,
           color = tool,
           shape=grepl("CONS", tool))
       ) +
  facet_wrap("rank", scales = "fixed") +
  geom_point(size=3, alpha=0.8) +
  scale_colour_manual(
    values=c("black", RColorBrewer::brewer.pal(12,"Paired"))
  ) +
  labs(x="2 - L1 norm error", y="16 - weighted UniFrac error", shape="", color="Tool") +
  ylim(0, NA) +
  scale_shape_discrete(guide = "none") +
  theme_cowplot(font_size = 17) +
  theme(panel.spacing.x = unit(1.5, "lines"))
p2
plot_grid(p1+theme(legend.position = "none"), p2+theme(legend.position = "none"), ncol=2)
ggsave2("../figures/profiling-strain_madness-tool_comparison.pdf", width = 15, height = 5)
# ggsave2("../figures/profiling-marine-tool_comparison.pdf", width = 15, height = 5)
# plot_grid(get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 0), legend.position = "bottom", legend.justification = "center")))
# ggsave2("../figures/profiling-legend-CAMI2.pdf", width = 14, height = 1)