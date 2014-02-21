#!/usr/bin/env Rscript

library(ggplot2)
agg <- read.delim('lcfit_plots_agg.txt', sep='\t')
p <- ggplot(agg, aes(x=t, y=rss, shape=model_name, color=model_name)) +
  geom_point(size=1.2) +
  ggtitle("RSS by branch length") +
  xlab("ML Branch length") +
  ylab("RSS") +
  theme_bw() +
  xlim(0,2)
ggsave('rss_bl.svg', width=7,height=5)

maxima <- read.csv('lcfit_maxima.csv', as.is=TRUE)
p <- ggplot(maxima, aes(x=t, y=t_hat, color=model_name)) +
    geom_point(size=1.2) +
    geom_smooth(method='lm', se=FALSE) +
    geom_abline(slope=1, intercept=0, linetype='dashed', color='grey') +
    xlab('t') +
    ylab(expression(hat(t))) +
    ggtitle("lcfit versus ML branch length") +
    facet_wrap(~model_name) +
    theme(legend.position='none')
ggsave('ml_lcfit_comparison.svg')

