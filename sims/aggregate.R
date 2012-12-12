#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
#agg <- read.delim('lcfit_plots_agg.txt', sep='\t')
#p <- ggplot(agg, aes(x=t, y=rss, shape=model_name, color=model_name)) +
  #geom_point(size=1.2) +
  #ggtitle("RSS by branch length") +
  #xlab("ML Branch length") +
  #ylab("RSS") +
  #theme_bw() +
  #xlim(0,2)
#ggsave('rss_bl.svg', width=7,height=5)

maxima <- read.csv('lcfit_maxima.csv', as.is=TRUE)

m <- melt(maxima, id.vars=c('node_id', 'model_name', 'ml_tree'), measure.vars=c('brent_n', 'lcfit_n'))
m <- transform(m, variable=ifelse(variable=='brent_n', 'Brent', 'lcfit'))

p <- ggplot(m, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  labs(x='Fit method', y='# of LL calls') +
  theme(legend.position='none') +
  facet_wrap(~model_name) +
ggsave('ml_lcfit_nsteps.svg')

p <- ggplot(maxima, aes(x=lcfit_t, y=brent_t, color=model_name)) +
  geom_abline(slope=1, intercept=0) +
  geom_point(size=1.2) +
  ggtitle("ML branch length estimates") +
  facet_wrap(~model_name) +
  theme(legend.position='none')

ggsave('ml_lcfit_comparison.svg')
