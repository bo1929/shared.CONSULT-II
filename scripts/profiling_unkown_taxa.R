require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
library(tidyr)

pm <- read_tsv("../results/all_tools-profiling_evaluation-CAMI1_hc.tsv")
pm <- pm %>% filter(metric %in% c("Bray-Curtis distance", "Shannon equitability", "Sum of abundances"))
pm <- pm %>% filter(rank != "superkingdom") %>% filter(rank != "strain")
pm <- pm %>%
  group_by(sample, metric) %>%
  mutate(dvalue = abs(value - value[tool == "Gold standard"]))
pm <- pm %>% filter(tool != "Gold standard")
pm$rank <- factor(
  pm$rank,
  levels = c("phylum", "class", "order", "family", "genus", "species")
)
pm <- pm %>% filter(tool %in% c("CONSULT-II (Eq. 7 without unclassified)", "CONSULT-II (Eq. 7)"))
pm$tool[pm$tool == "CONSULT-II (Eq. 7 without unclassified)"] <- "without unclassified abundance"
pm$tool[pm$tool == "CONSULT-II (Eq. 7)"] <- "with unclassified abundance"
pm$tool <- factor(pm$tool)
pm <- pm %>%
  group_by(rank, tool, metric) %>%
  summarise(lower = min(dvalue), upper = max(dvalue), p = mean(dvalue))
pm$metric <- as.factor(pm$metric)

bc_plot <- ggplot(pm %>% filter(metric == "Bray-Curtis distance"),
                  mapping = aes(x = p, y = reorder(tool, p), color = tool)) +
  geom_pointrange(size =  0.7, linewidth = 1, mapping = aes(xmin = lower, xmax = upper)) +
  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 16) +
  labs(color = "Approach", y = "", x = "Bray-Curtis dissimilarity to true profile") +
  scale_colour_brewer(palette = "Set1", type = "seq") +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 1, size=9)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  scale_shape_manual(values=c(1,16))
bc_plot
se_plot <- ggplot(pm %>% filter(metric == "Shannon equitability"), mapping = aes(x = p, y = tool, color = tool)) +
  geom_pointrange(size =  0.7, linewidth = 1, mapping = aes(xmin = lower, xmax = upper)) +
  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 16) +
  labs(color = "Approach", y = "", x = "Shannon's equitability (absolute difference to gold standard)") +
  scale_colour_brewer(palette = "Set1") +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 1, size=9)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))
se_plot
figperf <- plot_grid(
  bc_plot + theme(legend.box.margin = margin(0, 0, 0, 0), legend.position = "top"),
  ggplot() + theme_void(),
  se_plot + theme(legend.position = "none"),
  ncol = 1,
  rel_heights = c(1.2, 0.2, 1)
)
ggsave2("../figures/profiling_unclassified_comparison.pdf", width = 10, height = 6)

ut <- read_tsv("../results/unclassified_abundances.tsv")
ut$rank <- factor(
  ut$rank,
  levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)
ut <- ut %>%
  group_by(rank) %>%
  summarise(lower=min(100-portion), upper=max(100-portion), p=mean(100-portion))

figut <- ggplot(ut) + 
  geom_bar(aes(x=rank, y=p, fill=rank), color="black", stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rank, ymin=lower, ymax=upper), width=0.4,  alpha=0.9, size=1.3) +
  theme_cowplot() +
  labs(x="Taxonomic rank", y="Unclassified taxon abundance") +
  ylim(c(0,100)) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_brewer(palette="Paired", guide="none")
ggsave2("../figures/profiling_unclassified_taxon_proportions.pdf", width = 4, height = 4)

plot_grid(figperf, plot_grid(ggplot() + theme_void(), figut, ggplot() + theme_void(), ncol=1, rel_heights = c(1,2,1)), rel_widths = c(4, 2.5), ncol=2, labels=c("a", "b"))
ggsave2("../figures/profiling_with_unclassified.pdf", width = 13, height = 6)