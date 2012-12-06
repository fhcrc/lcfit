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
