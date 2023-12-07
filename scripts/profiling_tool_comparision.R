require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
library(tidyr)

pm <- read_tsv("../results/all_tools-profiling_evaluation-CAMI1_hc.tsv")
# pm <- pm %>% filter(metric == "Bray-Curtis distance" | metric == "Shannon equitability")
pm <- pm %>% filter(rank != "superkingdom") %>% filter(rank != "strain")
pm <- pm %>%
  group_by(sample, metric) %>%
  mutate(dvalue = abs(value - value[tool == "Gold standard"]))
pm <- pm %>% filter(tool != "Gold standard")
pm$rank <- factor(
  pm$rank,
  levels = c("phylum", "class", "order", "family", "genus", "species")
)
pm <- pm %>% filter(tool != "CONSULT-II")
pm$tool[pm$tool == "CONSULT-II (rn)"] <- "CONSULT-II"
pm$tool <- factor(
  pm$tool,
  levels = c("CLARK", "CONSULT-II", "Bracken")
)
pm <- pm %>%
  group_by(rank, tool, metric) %>%
  summarise(lower = min(dvalue), upper = max(dvalue), p = mean(dvalue))
pm$metric <- as.factor(pm$metric)

bc_plot <- ggplot(pm %>% filter(metric == "Bray-Curtis distance"), mapping = aes(x = p, y = tool, color = tool)) +
  geom_pointrange(size = 0.7, linewidth = 1, mapping = aes(xmin = lower, xmax = upper)) +
  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 17) +
  labs(color = "Tool", y = "", x = "Bray-Curtis dissimilarity to true profile") +
  scale_colour_brewer(palette = "Dark2") +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 1, size=9)) +
  theme(aspect.ratio = 0.4, panel.spacing.x = unit(1.5, "lines")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
bc_plot
se_plot <- ggplot(pm %>% filter(metric == "Shannon equitability"), mapping = aes(x = p, y = tool, color = tool)) +
  geom_pointrange(size =  0.7, linewidth = 1, mapping = aes(xmin = lower, xmax = upper)) +
  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 17) +
  labs(color = "Tool", y = "", x = "Shannon's equitability (absolute difference to gold standard)") +
  scale_colour_brewer(palette = "Dark2") +
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
  se_plot + theme(legend.box.margin = margin(0, 0, 0, 12))
)
plot_grid(prow, legend, rel_widths = c(3, .45))

ggsave2("../figures/profiling_tool_comparison.pdf", width = 13.5, height = 4)

tm <- pivot_wider(
  pm, id_cols = c(rank, tool),
  names_from = metric,
  values_from = p)
ggplot(tm %>% filter(rank != "na"), aes(x = Completeness, y = Purity, color = tool)) +
  geom_point(size=3) +
  facet_wrap("rank") +
  xlim(0, NA) +
  theme_cowplot(font_size = 17) +
  scale_colour_brewer(palette = "Dark2")

rm <- (tm %>% filter(rank == "genus"))
rm$`Weighted UniFrac error` <- (tm %>% filter(is.na(rank)))$`Weighted UniFrac error`
ggplot(rm, aes(x = 2-`L1 norm error`, y = 16-`Weighted UniFrac error`, color = tool)) +
  geom_point(size=3) +
  xlim(0, NA) +
  ylim(0, 16) +
  theme_cowplot(font_size = 17) +
  scale_colour_brewer(palette = "Dark2")
