#!/usr/bin/env Rscript

##' Aggregate error mesures for various types of fittig procedures across all the simulated trees.
##' 
##' usage: aggfit.R <control_file> ...
##' e.g.   aggfit.R runs/10/0/JTT92/gamma4-0.2/control.json runs/10/0/JTT92/uniform/control.json
##'
##' Usually you would run this 
##' 

library(ggplot2)
library(gridExtra)
library(plyr)
library(reshape2)
library(RJSONIO)
suppressPackageStartupMessages(library(distrEx))
suppressPackageStartupMessages(library(entropy))
source("analysis/utils.R")

theme_set(theme_bw(16) + theme(strip.background=element_blank(),
                               legend.key = element_blank(),
                               legend.position='bottom'))

args <- commandArgs(TRUE)
stopifnot(length(args) >= 1)
control_paths <- args
#control_paths <- c("runs/10/0/LG08/uniform/control.json")

## control_paths <- c("runs/10/0/JTT92/gamma4-0.2/control.json")
## control_paths <- c("runs/10/0/JTT92/gamma4-0.2/control.json","runs/10/0/JTT92/uniform/control.json")

# load all the contents of control files into a list of named lists.
# why?!?   What is achieved by loading all this stuff into memory?
# why not process one file at a time?
controls <- llply(control_paths, fromJSON)

# calculate effective sample size
ess <- function(ll) { ll <- ll - max(ll); exp(-log(sum(exp(2 * ll))) + 2 * log(sum(exp(ll)))) }


# I think this is code may be wrong.   This guards against overflow when the
# log values are positive.   When the log values are negative (e.g. -2719) you want to
# scale by the SMALLEST value (i.e. the most negative value), not the largest.
# Fortunately if your range of values of not that great, then
# this bug doesn't matter, but it will eventually bite you.
logsumexp <- function(x) {
  log(sum(exp(x - max(x)))) + max(x)
}


kl_log <- function(log_freqs1, log_freqs2) {
  log_freqs1 <- log_freqs1 - logsumexp(log_freqs1)
  log_freqs2 <- log_freqs2 - logsumexp(log_freqs2)
  LR <- log_freqs1 - log_freqs2
  sum(exp(log_freqs1) * LR) / log(2)
}


# calculate summary statistics - the residual squared difference between
# the fitting procedures (lcfit and spline) and the actual likelihood.
compare_bls <- function(fit_ll, bpp_ll) {

  ## bpp_l <- exp(bpp_ll - max(bpp_ll))
  ## bpp_l <- bpp_l / sum(bpp_l)
  ## fit_l <- exp(fit_ll - max(fit_ll))
  ## fit_l <- fit_l / sum(fit_l)
  ## kl=KL.plugin(bpp_l, fit_l, unit = 'log2')
  
  w <- fit_weights(bpp_ll)  # weights for weighted RSS

  data.frame(
             wrss=sum((w*(fit_ll- bpp_ll))^2)
             )
}

calculate_lcfit <- function(bls, model, fit, weighted=F, keep=0) {
    sample.points <- fit[, c('branch_length', 'll')]
    names(sample.points) <- c('x', 'y')
    
    # Erick asks to keep only the four-highest likelihood points.
    if (keep > 0) {
        sample.points <- sample.points[order(sample.points$y, decreasing=T)[1:keep],]
    }
    model <- fit_model(model, sample.points, weighted)
    sapply(bls[['branch_length']], lcfit, model=model)
}

# sample extra points
samplePoints <- function(node, nextra=0) {
    node$fit$name <- 'Bio++'
    fit <- node$fit
    # sample some extra points
    if (nextra > 0) {
        # compute the categorical distribution of points
        bpp_ll <- node$bls[['bpp_ll']]
        bpp_l <- exp(bpp_ll - max(bpp_ll))
        bpp_l <- bpp_l / sum(bpp_l)

        # propose extra sample points.
        # reject the proposal if any of the new points duplicate existing points.
        #
        nrejected <- 0
        repeat {
            # sample w/o replacement according to the distribution calculated above
            i <- sample(1:nrow(node$bls), nextra, replace = FALSE, prob=bpp_l)
            # add the new points to the existing ones.
            extra.fit <- data.frame(node_id=0, branch_length=node$bls[i,"branch_length"], ll=node$bls[i,"bpp_ll"], name="extra")

            # Reject if any duplicates
            if (!any(extra.fit$branch_length %in% node$fit$branch_length))
                break
            nrejected <- nrejected + 1
            # too many rejections result in failure.

            if (nrejected > 200) {
                extra.fit <- NULL
                break
            }
        }

        if (is.null(extra.fit))
            # we didn't successfully generate a good set of extra points.
            # fail out of the routine.
            fit <- NULL
        else
            # succcess!
            fit <- rbind(node$fit, extra.fit)
    }

    if (!is.null(fit)) {
        if (sum(fit$branch_length == max(fit$branch_length)) != 1) {
            message("in samplePoints: Saw duplicated branch_lengths")
        }
    }
    return(fit)
}

##' Reads in simulation data for a particular tree and aggregates into a single R data structure
##'
##' This routine reads a JSON format control file to retrieve the names
##' of various simulation output files.  It reads the data from those output files and aggregates it
##' by node_id.
##' 
##' @param path path to a control file, e.g. "runs/10/0/JTT92/gamma4-0.2/control.json"
##' @return a list of lists, with each sublist containing the data for a single node
##' in the tree.  
##' @author chris
readSimulationData <- function(path) {
    p <- fromJSON(path)
    bls_file <- p$lcfit[1]
    maxima_file <- p$lcfit[2]
    fit_file <- p$lcfit[3]
    ml_est_file <- p$lcfit[4]
    
    bls <- read.csv(bls_file, as.is=TRUE)
    maxima <- read.csv(maxima_file, as.is=TRUE)
    fit <- read.csv(fit_file, as.is=TRUE)
    bls.list <- dlply(bls, .(node_id), function(d) d)
    fit.list <- dlply(fit, .(node_id), function(d) d)
    maxima.list <- dlply(maxima, .(node_id), function(d) d)
    m <- mapply(list, fit.list, bls.list, maxima.list, SIMPLIFY=FALSE)
    m <- lapply(m, setNames, c("fit", "bls", "model"))
}


##' plot estimates likelihood curves against actual likelihood.
##'
##' Given a path to a JSON control file, read in the data pointed to by the control file
##' and plot the various likelihood estimates on the actual likelihood as calulated by BPP.
##' @param path character path to control file, i.e. "runs/10/0/JTT92/uniform/control.json"
##' @param node_id identifier of node in tree, i.e. 4
##' @param nextra, number of extra points to fit, usually in the range 0-4
##' @return a data frame of error measures
##' @author chris
showplot <- function(path, node_id, nextra) {
    print(path)
    tree <- readSimulationData(path)
    node <- tree[[node_id]]

    fit <- samplePoints(node, nextra=nextra)
    if (is.null(fit)) {
        message(sprintf("cannot sample points for path %s, node %d", path, node$bls$node_id[[1]]))
    } else
        node$fit <- fit

    node$bls$fit_ll <- calculate_lcfit(node$bls, node$model, node$fit, weighted=F)
    node$bls$wfit_ll <- calculate_lcfit(node$bls, node$model, node$fit, weighted=T)
    node$bls$t4fit_ll <- calculate_lcfit(node$bls, node$model, node$fit, weighted=F, keep=4)

    bls <- melt(node$bls, id.vars=1:2)
    bl.cutoff <- 0.07

    p <- ggplot() +
         ylab("Log likelihood")


    mainplot <- p + geom_line(aes(x=branch_length, y=value, color=variable, linetype=variable),
                          data=subset(bls, branch_length < bl.cutoff)) +
                ggtitle(sprintf("Node #%s",  bls$node_id[1])) +
                geom_point(aes(x=branch_length, y=ll), data=subset(node$fit, branch_length < bl.cutoff)) +
                scale_color_discrete(name="", breaks=c('fit_ll', 'wfit_ll', 't4fit_ll'), labels=c('unweighted', 'weighted', 'top4')) +
                guides(color="legend", shape=FALSE, linetype=FALSE) +
                theme(legend.position="bottom") +
                ylab("Log likelihood") 
                #annotate(geom="text", label="foo", x=.5, y=.9)


    subplot <- p + 
        geom_line(aes(x=branch_length, y=value, color=variable, linetype=variable),
                  data=bls) +
        geom_point(aes(x=branch_length, y=ll), data=node$fit) +
        ggtitle(sprintf("Full plot",  bls$node_id[1])) +
        theme(plot.title = element_text(size = 10)) +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              # axis.ticks=element_blank(),
              plot.background=element_blank(),
              #panel.background=element_blank(),
              #panel.border=element_blank(),
              #panel.border = element_rect(linetype = "dashed", color = "blue")),
              legend.position="none"
              )

    
    vp <- viewport(width = 0.4, height = 0.4, x = 0.15, y = 0.17, just = c("left", "bottom"))

    print(mainplot)
    print(subplot, vp=vp)

    ldply(c(weighted="weighted", unweighted="unweighted", top4="top4"),
          function(analysis) {
              ll <- switch(analysis,
                           weighted = node$bls$wfit_ll,
                           unweighted = node$bls$fit_ll,
                           top4 = node$bls$t4fit_ll)
              df <- compare_bls(ll, node$bls$bpp_ll)
          })

}

showplot <- function(path, node_id, nextra) {
    print(path)
    tree <- readSimulationData(path)
    node <- tree[[node_id]]

    fit <- samplePoints(node, nextra=nextra)
    if (is.null(fit)) {
        message(sprintf("cannot sample points for path %s, node %d", path, node$bls$node_id[[1]]))
    } else
        node$fit <- fit

    node$bls$fit_ll <- calculate_lcfit(node$bls, node$model, node$fit, weighted=F)
    node$bls$wfit_ll <- calculate_lcfit(node$bls, node$model, node$fit, weighted=T)
    node$bls$t4fit_ll <- calculate_lcfit(node$bls, node$model, node$fit, weighted=F, keep=4)

    bls <- melt(node$bls, id.vars=1:2)

    p <- ggplot( ) +
         geom_line(aes(x=branch_length, y=value, color=variable, linetype=variable), data=bls) +
         ggtitle(sprintf("Node #%s",  bls$node_id[1])) +
         geom_point(aes(x=branch_length, y=ll), data=node$fit) +
         ylab("Log likelihood")

    p <- p + scale_color_discrete(name="", breaks=c('fit_ll', 'wfit_ll', 't4fit_ll'), labels=c('unweighted', 'weighted', 'top4')) 
    p <- p + guides(color="legend", shape=FALSE, linetype=FALSE) 
    p <- p + theme(legend.position="bottom")
    print(p)

    ldply(c(weighted="weighted", unweighted="unweighted", top4="top4"),
          function(analysis) {
              ll <- switch(analysis,
                           weighted = node$bls$wfit_ll,
                           unweighted = node$bls$fit_ll,
                           top4 = node$bls$t4fit_ll)
              df <- compare_bls(ll, node$bls$bpp_ll)
          })

}

# For each list in the list of control structures,
# syntax for naming the .id column taken from,
# http://stackoverflow.com/a/15559709/1135316

aggregate <- function(control_paths) {
    ldply(setNames(control_paths, control_paths), .id="tree", 
               function(path) {
                   message(path)
                   tree <- readSimulationData(path)
                   names(tree) <- sapply(tree, function(node) node$model$node_id)
                   # for each node
                   #     for N extra points
                   #         for unweighted, weighted, top4
                   #             calculate fit
                   #             compare to empirical
                   ldply(tree, .id="node_id", .parallel=TRUE,
                                 function(node) {

                                     # estimate log-likelihood with and without extra sample points
                                     tmp <- llply(0:4,
                                                  function(nextra) {
                                                      fit <- samplePoints(node, nextra=nextra)
                                                      if (is.null(fit)) {
                                                          message(sprintf("cannot sample points for path %s, node %d", path, node$bls$node_id[[1]]))
                                                      } else
                                                          node$fit <- fit
                                                      # estimate log-likelihood using a variety of methods,
                                                      # and return an error between the estimate and empircal log-likelihood
                                                      rss <- ldply(c(weighted="weighted", unweighted="unweighted", top4="top4"),
                                                                   function(analysis) {
                                                                       ll <- switch(analysis,
                                                                                    weighted = do.call(calculate_lcfit, c(node, weighted=T)),
                                                                                    unweighted = do.call(calculate_lcfit, c(node, weighted=F)),
                                                                                    top4 = do.call(calculate_lcfit, c(node, weighted=F, keep=4)))
                                                                       df <- compare_bls(ll, node$bls$bpp_ll)
                                                                   })
                                                      rownames(rss) <- paste0(rss$.id, "_", nextra)
                                                      rss <- rss[,-c(1), drop=FALSE]
                                                      # melting a matrix gives access to rownames
                                                      # http://stackoverflow.com/a/19826856/1135316
                                                      rss.m <- melt(t(rss))
                                                      rss <- dcast(rss.m, . ~ ...)
                                                      rss <- rss[,-c(1),drop=FALSE]
                                                      
                                                  })
                                     tmp <- do.call(cbind, tmp)
                                 })
                   # error <- ldply(bls.list, function(d) compare_bls(d))
               })
}

library(doMC)
registerDoMC(4)

control_paths <- args
#control_paths <- c("runs/10/0/JTT92/gamma4-0.2/control.json","runs/10/0/JTT92/uniform/control.json","runs/10/0/Binary-4.0/gamma4-0.2/control.json")
#if (is(try(load("aggdata.RData")), "try-error")) {
     aggdata <- aggregate(control_paths)
     save(aggdata, file="aggdata.RData")
#}


## at this point aggregate is a data frame that looks like this,

##                                      tree node_id kl_weighted_0 kl_unweighted_0
## 1 runs/10/0/JTT92/gamma4-0.2/control.json       0  7.594155e-06    6.757892e-06
## 2 runs/10/0/JTT92/gamma4-0.2/control.json       1  3.195024e-06    3.195024e-06
## 3 runs/10/0/JTT92/gamma4-0.2/control.json       2  6.245404e-06    6.670254e-06

# There are many more columns to the right representing additional error measures.
# a column name like " kl_weighted_0" is made up on the following components,
# <measure>_<fitting>_<samples>
# <measure> indicates which measure of error is used, usually just Kullback–Leibler divergence (kl) or weighted residual sum of squares (wrss)
# <fitting> indicates which error funciton was minimized when fitting the lcfit function to the sample points.
#           'weighted' for weighted rss, 'unweighted' for unweighted rss, top4 for unweighted using only the four highest-likelihood sample points.
# <samples> indicates how many additional sample points were included beyond those selected by the basic lcfit algorithm
#

m <- melt(aggdata, id=c("tree", "node_id"))
m <- cbind(m, colsplit(m$variable, names = c("measure", "fitting", "nextra"), pattern="_"))

tmp <- m[m$measure=="wrss" & m$fitting == 'weighted',]

p1 <- ggplot() +
      geom_boxplot(aes(x=factor(nextra), value), data=tmp) +
      xlab("number of extra points") +
      ggtitle("weighted rss error between weighted lcfit\n and empirical likelihoods") +
      theme(legend.position='bottom')
ggsave('wrss_weighted_vs_nextra.png', p1, width=14, height=9, dpi=72)

p1

tmp <- m[m$measure=="wrss" & m$fitting == 'weighted' & m$nextra > 2,]

p1 <- ggplot() +
      geom_boxplot(aes(x=factor(nextra), value), data=tmp) +
      xlab("number of extra points") +
      ggtitle("weighted rss error for 2 or more extra points\nbetween weighted lcfit and empirical likelihoods") +
      theme(legend.position='bottom')
ggsave('wrss_weighted_vs_nextra.gt.2.png', p1, width=14, height=9, dpi=72)

p1

tmp <- m[m$measure=="wrss" & m$fitting == 'weighted' & m$nextra == 4,]

p1 <- ggplot() +
      geom_boxplot(aes(x=factor(nextra), value), data=tmp) +
      xlab("number of extra points") +
      ggtitle("weighted rss error for 4 extra points\nbetween weighted lcfit and empirical likelihoods") +
      theme(legend.position='bottom')
ggsave('wrss_weighted_vs_nextra4.png', p1, width=14, height=9, dpi=72)

p1



tmp <- m[m$measure=="wrss" &  m$nextra == 2,]
p1 <- ggplot() +
      geom_boxplot(aes(x=factor(fitting), value), data=tmp) +
      xlab("minimization function") +
      ggtitle("weighted rss error for 2 extra points\n fitting function compared against empirical likelihoods") +
      theme(legend.position='bottom')
ggsave('wrss_nextra2_vs_measure.png', p1, width=14, height=9, dpi=72)
p1

tmp <- m[m$measure=="wrss" &  m$nextra == 3,]
p1 <- ggplot() +
      geom_boxplot(aes(x=factor(fitting), value), data=tmp) +
      xlab("minimization function") +
      ggtitle("weighted rss error for 3 extra points\n fitting function compared against empirical likelihoods") +
      theme(legend.position='bottom')
ggsave('wrss_nextra3_vs_measure.png', p1, width=14, height=9, dpi=72)
p1

showplot <- function(path, node_id, nextra) {

with(tmp[order(tmp$value, decreasing=T)[1], ],
     showplot(as.character(tree), node_id, nextra))
