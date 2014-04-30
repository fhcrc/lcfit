#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)

theme_set(theme_bw(16) + theme(strip.background = element_blank()))

agg <- read.csv('agg_mltol.csv', as.is = TRUE)

to_plot <- transform(agg, tol = factor(tolerance),
                     n_eval_diff = n_eval - n_eval_brent)

x_lev <- as.character(to_plot$tolerance)
x_lab <- sprintf('10^%d', log10(to_plot$tolerance))
x_lab <- parse(text = x_lab)


p1 <- ggplot(to_plot, aes(x = tol, y = n_eval_diff)) +
  geom_boxplot() +
  xlab("Tolerance") +
  geom_hline(yintercept = 0.0, linetype = 'dashed') +
  scale_x_discrete(breaks = x_lev, labels = x_lab) +
  ylab("Branch log-likelihood calls: lcfit - Brent")
ggsave('branch_likelihood_calls_by_tol.svg', width = 8, height = 5)

p2 <- ggplot(agg, aes(x = ml_brent, y = ml_est)) +
  geom_abline(slope = 1) +
  geom_point() +
  coord_equal() +
  facet_wrap(~tolerance)
p2
