require(ggplot2)
require(reshape2)
require(readr)
require(cowplot)
require(dplyr)

rc <- read_tsv("../results/resource_tool_comparison.tsv")
rc <- rc %>% filter(Tool == "CLARK" | Tool == "CONSULT-II" | Tool == "Kraken-II")

ggplot(melt(rc[, 0:4], id.vars = 1:2)) +
  geom_col(aes(x=Tool, colour=variable, fill=variable, y=value), position = position_dodge2(padding = 0.1)) +
  theme_cowplot(font_size = 17) +
  labs(fill = "", colour = "", x = "") +
  theme(aspect.ratio = 0.9) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_y_continuous(name = "Memory (GB)", c(16, 32, 64, 128), sec.axis = sec_axis(~., name="Run-time (minutes)", breaks = c(15, 60, 120))) +
  scale_colour_brewer(palette = "Accent", labels=c("Memory_GB" = "Memory", "Time_min"="Time")) +
  scale_fill_brewer(palette = "Accent", labels=c("Memory_GB" = "Memory", "Time_min"="Time"))
ggsave2("../figures/resource_tool_comparison.pdf")