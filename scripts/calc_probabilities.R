require(rlang)
require(ggplot2)
require(cowplot)

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
  geom_function(fun = function(x, w, s) pmin(w / (pmax(x + w - s, w)) + 1 / s**2, 1), args = list(w = 4, s = 5), linewidth = 1, aes(color = "old soft-LCA function")) +
  geom_function(fun = function(x, s) 1 / log2(((x - 1) / s)**2 + 2), args = list(s = 6), linewidth = 1, aes(color = "new soft-LCA function")) +
  scale_color_discrete(name = "functions") +
  theme_minimal_grid(font_size = 15) +
  theme(aspect.ratio = 0.8)
ggsave2("../figures/probability_LCA_update.pdf", width = 4.25, height = 3)



require(scales)
SL <- 32
p <- 3

f <- function(x, k, l, SL) {
  1 - (1 - (1 - x / SL)^k)^l
}

f(3, 15, 2, 32)

K <- function(l, p, alpha, L) {
  round(log(1 - (1 - alpha)^(1 / l)) / log(1 - p / L))
}

K(4, 2, 0.95, 10)

alpha <- 0.1
ggplot(data.frame(x = c(1, 16)), aes(x)) +
  theme_classic() +
  mapply(FUN = function(L, K) {
    stat_function(fun = function(x, k, l, SL) (1 + 150 - SL) * f(x, k, l, SL), args = list(k = K, l = L, SL = SL), aes(color = as.factor(K), linetype = as.factor(L)))
  }, L = rep(c(2, 4, 1), each = 5), K = c(6:10) * 2 - 1) +
  scale_y_log10(lim = c(0.1, 100)) +
  scale_linetype_manual(values = c(2, 1, 3), name = expression(l)) +
  scale_x_continuous(labels = function(x) paste(x, percent(x / SL), sep = "\n")) +
  scale_color_brewer(palette = "Dark2", name = "h") +
  ylab("Expected number of matches") +
  xlab("Hamming distance") +
  stat_function(fun = function(x) (1 + 150 - SL) * f(x, 35 - 7, 1, 35) / 2, color = "black", aes(linetype = as.factor(1))) +
  theme_minimal_grid(font_size = 15) +
  theme(aspect.ratio = 1)
ggsave2("../figures/expected_num_matches.pdf", width = 4.25, height = 4)
