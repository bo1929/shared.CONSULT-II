require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
require(reshape2)
require(latex2exp)

eval0 <- read_csv("../results/CONSULTII-evaluation-bacteria-w4s5_th05_c00.csv")
eval0 <- arrange(eval0, vote)
eval0 <- melt(eval0[, 2:11], id.vars = 1:4)
eval0_p <- eval0[eval0$value %in% c("TP", "FP"), ]
eval0_n <- eval0[!eval0$value %in% c("TP", "FP"), ]

ggplot(aes(x = vote, color = variable, linetype = value), data = eval0_p) +
  stat_ecdf(n = 1000) +
  scale_y_continuous(name = "ECDF") +
  scale_x_continuous(name = TeX("Total vote ($\\bar{v}$)"), trans = "log10") +
  scale_linetype_manual(name = "", values = c(1, 3)) +
  annotate("rect", xmin = 0.003, xmax = 0.01, ymin = 0, ymax = 1, alpha = .1) +
  annotate("rect", xmin = 0.003, xmax = 0.03, ymin = 0, ymax = 1, alpha = .1) +
  scale_color_brewer(palette = "Dark2", name = "Rank", direction=-1) +
  theme_cowplot(font_size = 17) +
  theme(aspect.ratio = 1)
ggsave2("../figures/total_vote_impact.pdf")

ggplot(aes(x = voteNormalized, color = variable), data = eval0_p) +
  stat_ecdf(n = 1000) +
  scale_y_continuous(name = "ECDF") +
  facet_wrap(vars(value)) +
  scale_x_continuous(name = "Total vote (normalized)", trans = "log10") +
  scale_color_brewer(palette = "Dark2", name = "Rank", direction=-1) +
  theme_cowplot(font_size = 17) +
  theme(legend.position="none")
  theme(aspect.ratio = 1.25)
ggsave2("../figures/total_vote_impact.pdf")