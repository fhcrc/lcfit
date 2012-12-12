#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(reshape2)

main <- function(input_bls, input_maxima, input_fit, outfile) {
  d <- read.csv(input_bls, as.is=TRUE)
  maxima <- read.csv(input_maxima, as.is=TRUE)
  m <- melt(d, id.vars=1:2)

  level_names <- c('Bio++', 'lcfit')
  bl_translate <- c('bpp_ll', 'fit_ll')
  max_translate <- c('t', 't_hat')

  m <- transform(m, name=level_names[match(variable, bl_translate)])
  fit <- read.csv(input_fit, as.is=TRUE)

  pdf(outfile)
  p <- ggplot(maxima, aes(x=brent_t, y=lcfit_t)) +
    geom_abline(intercept=0, slope=1, linetype='dashed') +
    geom_point() +
    ggtitle("ML branch length: lcfit vs. Brent") +
    xlab(expression(t[brent])) +
    ylab(expression(t[lcfit]))
  print(p)
  p <- ggplot(melt(maxima, id.vars='node_id', measure.vars=c('brent_n', 'lcfit_n')),
              aes(x=variable, y=value, fill=variable)) +
    geom_boxplot() +
    ggtitle('Number of peeling recursions to fit ML branch length') +
    xlab('') +
    ylab('# of peeling recursions')
  print(p)


  d_ply(m, .(node_id), function(piece) {
    node_id <- piece$node_id[1]
    f <- subset(fit, node_id == piece$node_id[1])
    f$name <- 'Bio++'

    rss <- with(dcast(piece, node_id+branch_length~variable), sum((fit_ll - bpp_ll)^2))

    p <- ggplot(piece, aes(color=name, linetype=name)) +
        geom_line(aes(x=branch_length,
                      y=value), data=piece) +
        ggtitle(sprintf("Node #%s\nRSS=%f", piece$node_id[1], rss)) +
        geom_point(aes(x=branch_length, y=ll), data=f) +
        xlim(0, max(c(max(f$branch_length), 1)))
    print(p)
  })
  dev.off()

  ## Norms
  #norms <- ddply(d, .(node_id), function(piece) {
    #node_id <- piece$node_id[1]
    #bl <- piece$branch_length[1]
    #l1 <- with(piece, sum(abs(bpp_ll - fit_ll)))
    #l2 <- with(piece, sum((bpp_ll - fit_ll)^2))
    #data.frame(node_id=node_id, branch_length=bl,
               #l1=l1, l2=l2)
  #})
}

if(!interactive()) {
  args <- commandArgs(TRUE)
  if(length(args) != 4) {
    stop("usage: plot_fits.R <input_bls> <input_maxima> <input_fit> <outfile>")
  }
  input_bls <- args[1]
  input_maxima <- args[2]
  input_fit <- args[3]
  output <- args[4]
  main(input_bls, input_maxima, input_fit, output)
}
