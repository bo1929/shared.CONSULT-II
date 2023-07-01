require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
require(plotrix)

num_reads <- 66667
L <- 150
k <- 32

ms <- read_tsv("../results/CONSULTII-summary_match_distances-bacteria-w4s5.tsv")
ms$Rank <- tolower(ms$Rank)
ms$Rank <- factor(
  ms$Rank,
  levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)

ms <- ms %>% filter(Rank != "kingdom")

ms$Ratio <- ms$Count / num_reads / (L - k + 1)

ms_summary <- ms %>%
  group_by(Rank, Hamming_distance, Match, Bin) %>%
  summarize(avg = mean(Ratio), std = std.error(Ratio))

ggplot(ms_summary, aes(color = Rank, fill = Rank)) +
  facet_grid(rows = vars(Bin), cols = vars(Match)) +
  geom_line(aes(x = factor(Hamming_distance), y = avg, linetype = Rank, group = Rank), stat = "identity", linewidth = 1) +
  # geom_pointrange(aes(x = factor(Hamming_distance), y = avg, ymin = avg - std, ymax = avg + std), alpha = 0.9, size = 0.3) +
  geom_point(aes(x = factor(Hamming_distance), y = avg), alpha = 0.9, size = 2) +
  scale_y_continuous(trans = "log10", labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_cowplot(font_size = 16) +
  labs(shape = "Rank", colour = "Rank", x = "Hamming distance", y = "Avg. proportion of matched kmers") +
  scale_colour_brewer(palette = "Dark2") +
  # scale_shape_manual(values=c(0, 1, 2, 15, 16, 17))+
  theme(panel.spacing.x = unit(1.15, "lines")) +
  theme(legend.key.width = unit(3.5, "line")) +
  annotation_logticks(sides = "l") +
  theme(aspect.ratio = 1.35)
ggsave2("../figures/number_of_matches.pdf", width = 6.5, height = 8.5)