#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(reshape2)

main <- function(infile, outfile) {
  d <- read.csv(input, as.is=TRUE)
  m <- melt(d, id.vars=1:2)
  pdf(output)
  d_ply(m, .(node_id), function(piece) {
    p <- ggplot(piece, aes(x=branch_length,
                           y=value,
                           color=variable,
                           linetype=variable)) +
      geom_line() +
        opts(title=piece$node_id[1])
    print(p)
  })
  dev.off()

  # Norms
  norms <- ddply(d, .(node_id), function(piece) {
    node_id <- piece$node_id[1]
    bl <- piece$branch_length[1]
    l1 <- with(piece, sum(abs(bpp_ll - fit_ll)))
    l2 <- with(piece, sum((bpp_ll - fit_ll)^2))
    data.frame(node_id=node_id, branch_length=bl,
               l1=l1, l2=l2)
  })
}

if(!interactive()) {
  args <- commandArgs(TRUE)
  if(length(args) != 2) {
    stop("usage: plot_fits.R <infile> <outfile>")
  }
  input <- args[1]
  output <- args[2]
  main(input, output)
}
