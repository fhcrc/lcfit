#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(TRUE)

stopifnot(length(args) == 2)

bls <- read.csv(args[1], as.is=TRUE)

theme_set(theme_bw(12))
theme_update(axis.line=element_blank(),
             axis.text.x=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks=element_blank(),
             axis.title.x=element_blank(),
             axis.title.y=element_blank(),
             legend.position="none",
             panel.background=element_blank(),
             panel.border=element_blank(),
             panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             plot.background=element_blank(),
             strip.background=element_blank(),
             strip.text=element_blank())

p <- ggplot(bls, aes(x=branch_length, y=bpp_ll)) +
  geom_line() +
  facet_wrap(~node_id, scales='free', ncol=6)

ggsave(args[2], p, width=8, height=6)
