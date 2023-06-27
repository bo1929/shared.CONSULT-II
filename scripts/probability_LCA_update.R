require(rlang)
require(ggplot2)
require(cowplot)

theme_set(theme_minimal_grid(font_size = 10))

w <- 4
s <- 5

f <- as_function(~ pmin(w / (pmax(.x + w - s, w)) + 1 / s**2, 1))
pu <- function(x, w, s) pmin(w / (pmax(x + w - s, w)) + 1 / s**2, 1)

ggplot() +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(1, 10**5)
  ) +
  labs(y = "Probability of LCA update", x = "Number of genomes with that k-mer") +
  geom_function(fun = function(x, w, s) pmin(w / (pmax(x + w - s, w)) + 1 / s**2, 1), args = list(w = 4, s = 5), linewidth = 1) +
  theme_minimal_grid() +
  theme(aspect.ratio = 1)

ggsave2("../figures/probability_LCA_update.pdf")