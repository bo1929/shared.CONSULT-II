require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
require(ggthemes)

d2c <- read_tsv("../misc/dist-bacteria-to-closest.txt", col_names = FALSE)
d2c <- d2c %>%
  rowwise() %>%
  mutate(name = substr(X1, start = 1, stop = 10))
d2c <- d2c %>%
  arrange(X3) %>%
  mutate(name = factor(name, levels = name))

bb <- data.frame(X3 = c(0.001, 0.02, 0.06, 0.12, 0.16, 0.22))

ggplot(d2c, aes(x = name, y = X3, color = X3)) +
  geom_point(size = 1) +
  geom_hline(aes(yintercept = X3, color = X3), bb, linewidth = 1, alpha = 0.5) +
  theme_cowplot(font_size = 15) +
  labs(y = "Distance to the closest reference genome", x = "Genome", color = "") +
  theme(axis.text.x = element_text(angle = 90, size = 5), aspect.ratio = 0.45)
  # theme(axis.text.x = element_blank(), aspect.ratio=0.65)
ggsave2("../figures/novelty_bins-wlabels.pdf", width = 11, height = 6)