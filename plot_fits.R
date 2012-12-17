#!/usr/bin/env Rscript

library(animation)
library(ggplot2)
library(plyr)
library(reshape2)

animate_node_fit <- function(fit_log, fit_points, bpp_ll, outfile, title="Fit process") {
  cfn_loglike <- function(t, c, m, r, b) c*log((1+exp(-r*(t+b)))/2)+m*log((1-exp(-r*(t+b)))/2)
  x <- seq(0, 1, length.out=300)
  if(nrow(fit_log) > 100)
    fit_log <- fit_log[1:nrow(fit_log) %% 5 == 0,]
  fname <- basename(outfile)
  saveMovie({
      for (i in 1:nrow(fit_log)) {
         row <- fit_log[i, ]
         message(paste("Iteration", row$iter))
         r <- range(subset(bpp_ll, branch_length <=1, select=value))
         iter_fit <- data.frame(branch_length=x, ll=cfn_loglike(x, row$c, row$m, row$r, row$b), name='lcfit')
         p <- ggplot(bpp_ll, aes(color=name, linetype=name)) +
           geom_line(aes(x=branch_length, y=value), data=bpp_ll) +
           geom_point(aes(x=branch_length, y=ll), data=fit_points) +
           geom_line(aes(x=branch_length, y=ll), data=iter_fit) +
           theme_bw() + xlim(0, 1) + ylim(r[1]-50, r[2]+50) +
           ggtitle(paste(title, ": Iteration", row$iter))
         print(p)
      }
  }, interval=0.4, movie.name=fname, ani.width = 600, ani.height = 600, outdir=getwd())
}

read_fit_log <- function(path) {
  log <- read.csv(path, as.is=TRUE)
  node <- -1
  last <- 1000000
  log$node <- NA
  # Number log with node #
  for(i in 1:nrow(log)) {
    if(log$iter[i] != last + 1)
      node <- node + 1
    log$node[i] <- node
    last <- log$iter[i]
  }
  log
}


main <- function(input_bls, input_maxima, input_fit, input_fit_log, outfile) {
  d <- read.csv(input_bls, as.is=TRUE)
  maxima <- read.csv(input_maxima, as.is=TRUE)
  fit_log <- read_fit_log(input_fit_log)
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

    lg <- subset(fit_log, node==node_id)
    anim_out <- sprintf("node%03d.gif", node_id)
    animate_node_fit(lg, f, subset(piece, variable=='bpp_ll'), anim_out,
                     paste("Fit for node", node_id))
    file.rename(anim_out, paste(dirname(outfile), anim_out, sep='/'))
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
  if(length(args) != 5) {
    stop("usage: plot_fits.R <input_bls> <input_maxima> <input_fit> <fit_log> <outfile>")
  }
  input_bls <- args[1]
  input_maxima <- args[2]
  input_fit <- args[3]
  input_fit_log <- args[4]
  output <- args[5]
  main(input_bls, input_maxima, input_fit, input_fit_log, output)
}
