#!/usr/bin/env Rscript

library(animation)
library(ggplot2)
library(plyr)
library(reshape2)

# Setup theme
setup_theme <- function(size=19) {
  theme_set(theme_grey(size) %+replace%
      theme(panel.background  = element_rect(fill = "white", colour = NA),
        panel.border      = element_rect(fill = NA, colour = "grey50"),
        panel.grid.major  = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor  = element_line(colour = "grey98", size = 0.5)))
  theme <- theme_update(legend.key=element_blank(),
          legend.title=element_text(size=rel(1.1)),
          strip.background=element_blank())
  #theme <- modifyList(theme, list(
          #legend.key=element_blank(),
          #strip.background=element_blank(),
          #legend.title=element_text(size=20)))
  theme_set(theme)
}

animate_node_fit <- function(fit_log, fit_points, bpp_ll, outfile) {
  cfn_loglike <- function(t, c, m, r, b) c*log((1+exp(-r*(t+b)))/2)+m*log((1-exp(-r*(t+b)))/2)
  x <- seq(0, max(bpp_ll$branch_length), length.out=300)
  if(nrow(fit_log) > 100)
    fit_log <- fit_log[1:nrow(fit_log) %% 5 == 0,]
  fname <- basename(outfile)
  saveMovie({
      for (i in 1:nrow(fit_log)) {
         row <- fit_log[i, ]
         message(paste("Iteration", row$iter))
         r <- range(subset(bpp_ll, branch_length <=1, select=value))
         r_diff = diff(r) * 0.05
         iter_fit <- data.frame(branch_length=x, ll=cfn_loglike(x, row$c, row$m, row$r, row$b), name='lcfit')
         p <- ggplot(bpp_ll, aes(color=name, linetype=name)) +
           geom_line(aes(x=branch_length, y=value), data=bpp_ll) +
           geom_point(aes(x=branch_length, y=ll), data=fit_points) +
           geom_line(aes(x=branch_length, y=ll), data=iter_fit) +
           theme_bw() + ylim(r[1]-r_diff, r[2]+r_diff) +
           xlab("t") + ylab("Log-likelihood") +
           ggtitle(paste("Iteration", row$iter))
         print(p)
         # Double down on the last frame
         if(i == nrow(fit_log)) print(p)
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
    geom_abline(intercept=0, slope=1, linetype='dashed', color='grey80') +
    geom_point() +
    #ggtitle("ML branch length: lcfit vs. Brent") +
    xlab(expression(t[brent])) +
    ylab(expression(t[lcfit]))

  svg(sub('.pdf', '_ml_len.svg', outfile))
  print(p)
  dev.off()

  p <- ggplot(maxima, aes(x=brent_t, y=lcfit_fit_t)) +
    geom_abline(intercept=0, slope=1, linetype='dashed', color='grey80') +
    geom_point() +
    #ggtitle("ML branch length: lcfit vs. Brent") +
    xlab(expression(t[brent])) +
    ylab(expression(hat(t)[lcfit]))
  svg(sub('.pdf', '_fit_ml_est.svg', outfile))
  print(p)
  dev.off()

  melted_maxima <- melt(maxima, id.vars='node_id', measure.vars=c('brent_n', 'lcfit_n'))
  melted_maxima <- transform(melted_maxima, variable=gsub('_n', '', variable))

  p <- ggplot(melted_maxima, aes(x=variable, y=value, fill=variable)) +
    geom_boxplot() +
    #ggtitle('Number of peeling recursions to fit ML branch length') +
    xlab('Fit method') +
    ylab('# of peeling recursions') +
    theme(legend.position='none') +
    ylim(0, max(melted_maxima$value + 5))
  svg(sub('.pdf', '_peels.svg', outfile))
  print(p)
  dev.off()

  d_ply(m, .(node_id), function(piece) {
    node_id <- piece$node_id[1]
    f <- subset(fit, node_id == piece$node_id[1])
    f$name <- 'Bio++'

    rss <- with(dcast(piece, node_id+branch_length~variable), sum((fit_ll - bpp_ll)^2))

    p <- ggplot(piece, aes(color=name, linetype=name)) +
        geom_line(aes(x=branch_length,
                      y=value), data=piece) +
        geom_point(aes(x=branch_length, y=ll), data=f) +
        xlim(0, max(c(max(f$branch_length), 1)))
    print(p)

    # Animate
    lg <- subset(fit_log, node==node_id)
    anim_out <- sprintf("node%03d.gif", node_id)
    animate_node_fit(lg, f, subset(piece, variable=='bpp_ll'), anim_out)

    file.rename(anim_out, paste(dirname(outfile), anim_out, sep='/'))
  })
  dev.off()
}

if(!interactive()) {
  setup_theme()
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
