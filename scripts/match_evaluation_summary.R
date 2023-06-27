require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
require(plotrix)

ms <- read_csv("../results/CONSULTII-summary_matches-bacteria-w4s5.csv")
ms$Rank <- tolower(ms$Rank)
ms$Rank <- factor(
  ms$Rank,
  levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)

ms <- ms %>% filter(Rank != "kingdom")

ms_summary <- ms %>%
  group_by(Rank, Hamming_distance, Match) %>%
  summarize(avg = mean(Count), std = std.error(Count))

ggplot(ms_summary, aes(color = Rank, fill = Rank)) +
  facet_wrap(vars(Match), scales = "free_y") +
  geom_line(aes(x = factor(Hamming_distance), y = avg, linetype = Rank, group = Rank), stat = "identity", linewidth = 1) +
  geom_pointrange(aes(x = factor(Hamming_distance), y = avg, ymin = avg - std, ymax = avg + std), alpha = 0.9, size = 0.3) +
  scale_y_continuous(trans = "log10", labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_cowplot(font_size = 16) +
  labs(shape = "Rank", colour = "Rank", x = "Hamming distance", y = "Avg. number of matches") +
  scale_colour_brewer(palette = "Dark2") +
  theme(panel.spacing.x = unit(1.15, "lines")) +
  theme(legend.key.width = unit(2.5, "line")) +
  annotation_logticks(sides = "l") +
  theme(aspect.ratio = 1.5)
ggsave2("../figures/number_of_matches.pdf")