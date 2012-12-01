#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(reshape2)

args <- commandArgs(TRUE)
input <- args[1]
output <- args[2]

d <- read.csv(input, as.is=TRUE)
m <- melt(d, id.vars=1:2)
pdf(output)
d_ply(m, .(node_id), function(piece) {
  p <- ggplot(piece, aes(x=branch_length, y=value, color=variable)) +
    geom_line() +
    opts(title=piece$node_id[1])
  print(p)
})
dev.off()
