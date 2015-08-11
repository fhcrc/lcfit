#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)

theme_set(theme_bw(16) + theme(strip.background = element_blank()))

agg <- read.csv('agg_mltol.csv', as.is = TRUE)

to_plot <- transform(agg, tol = factor(tolerance),
                     tol_str = sprintf('10^%d', log10(tolerance)),
                     rel_err = abs(ml_brent - ml_est) / ml_brent,
                     abs_err = abs(ml_brent - ml_est),
                     n_eval_ratio = n_eval / n_eval_brent)

x_lev <- as.character(to_plot$tolerance)
x_lab <- sprintf('10^%d', log10(to_plot$tolerance))
x_lab <- parse(text = x_lab)

p1 <- ggplot(to_plot, aes(x = tol, y = n_eval_ratio)) +
  geom_boxplot() +
  xlab("Tolerance") +
  scale_x_discrete(breaks = x_lev, labels = x_lab) +
  ylab("Branch log-likelihood calls: lcfit / Brent") +
  ylim(0, 1)
ggsave('branch_likelihood_calls_by_tol.svg', width = 8, height = 5)

p2 <- ggplot(to_plot, aes(x = ml_brent, y = ml_est)) +
  geom_abline(slope = 1, color = 'grey', linetype = 'dashed') +
  geom_point() +
  coord_equal() +
  xlab("ML branch length (Brent)") +
  ylab("ML branch length (lcfit)") +
  facet_grid(~ tol_str, labeller = label_parsed)
ggsave('ml_err_by_tol.svg', width = 5, height = 4, dpi = 72)

