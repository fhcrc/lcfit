library(lcfit)
library(grid)

# routine for creating our own signaling conditions
# http://adv-r.had.co.nz/Exceptions-Debugging.html

# R doesn’t come with a built-in constructor function for conditions,
# but we can easily add one. Conditions must contain message and call
# components, and may contain other useful components. When creating a
# new condition, it should always inherit from condition and one of
# error, warning, or message.
condition <- function(subclass, message, call = sys.call(-1), ...) {
  structure(
    class = c(subclass, "condition"),
    list(message = message, call = call),
    ...
  )
}
is.condition <- function(x) inherits(x, "condition")

# You can then use tryCatch() to take different actions for different
# types of errors. In this example we make a convenient custom_stop()
# function that allows us to signal error conditions with arbitrary
# classes. In a real application, it would be better to have
# individual S3 constructor functions that you could document,
# describing the error classes in more detail.
custom_stop <- function(subclass, message, call = sys.call(-1), 
                        ...) {
  c <- condition(c(subclass, "error"), message, call = call, ...)
  stop(c)
}


# copied from lcift_R/src/lcfit.h
# there doesn't seem to be a good way to import enums from C into R
LCFIT_SUCCESS	= 0
LCFIT_MAXITER	= 1
LCFIT_ERROR	= 2
LCFIT_NOPROG	= 27
LCFIT_ETOLF	= 29

lcfit_errmap <- data.frame(status=c(LCFIT_SUCCESS, LCFIT_MAXITER, LCFIT_NOPROG, LCFIT_ETOLF),
                           msg=c("success", "maximum iterations exceeded", "iteration is not making progress towards solution", "cannot reach the specified tolerance in F"),
                            stringsAsFactors=FALSE)

lcfit_strerror <- function(status) {
    i <- match(status, lcfit_errmap$status)
    if (!is.na(i))
        return(lcfit_errmap$msg[i])
    else
        return("non-specific lcfit error")
}

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

    model.fitted <- fit_model(model, sample.points, weighted)
        
    if (model.fitted["status"] == LCFIT_MAXITER) {
        message(sprintf("WARN: node_id=%d weighted=%s keep=%d %s", bls$node_id[[1]], weighted, keep, lcfit_strerror(model.fitted["status"])));
    } else if (model.fitted["status"] != LCFIT_SUCCESS) {
        message(sprintf("ERROR: node_id=%d weighted=%s keep=%d %s", bls$node_id[[1]], weighted, keep, lcfit_strerror(model.fitted["status"])));
    }
    ll <- list(status=model.fitted[["status"]], ll=sapply(bls[['branch_length']], lcfit, model=model.fitted))
    return(ll)
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
            extra.fit <- data.frame(node_id=unique(node$fit$node_id), branch_length=node$bls[i,"branch_length"], ll=node$bls[i,"bpp_ll"], name="extra")

            # Reject if any duplicates
            if (!all(is.na(extra.fit$branch_length %in% node$fit$branch_length)))
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

    # Expect all lists aggregated by node_id should have the same length!
    if (length(unique(sapply(list(bls.list,fit.list,maxima.list), length))) != 1) {
        message(sprintf("Error: inconsistent number of nodes in each data file length"))
        message(sprintf("%s: %d", fit_file, length(fit.list)))
        message(sprintf("%s: %d", maxima_file, length(maxima.list)))
        message(sprintf("%s: %d", bls_file, length(bls.list)))
        stopifnot(length(unique(sapply(list(bls.list,fit.list,maxima.list), length))) != 1)
    }
    
    m <- mapply(list, fit.list, bls.list, maxima.list, SIMPLIFY=FALSE)
    m <- lapply(m, setNames, c("fit", "bls", "model"))

    return(m)
}


##' plot estimated likelihood curves against actual likelihood.
##'
##' Given a path to a JSON control file, read in the data pointed to by the control file
##' and plot the various likelihood estimates against the actual likelihood as calulated by BPP.
##' @param path character path to control file, i.e. "runs/10/0/JTT92/uniform/control.json"
##' @param node_id identifier of node in tree, i.e. 4
##' @param nextra, number of extra points to fit, usually in the range 0-4
##' @return a data frame of error measures
##' @author chris
showplot <- function(path, node_id, nextra, model=NULL, zoom=F) {
    message(sprintf("showplot: path=%s", path))
    message(sprintf("showplot: cwd=%s", getwd()))
    message(sprintf("showplot: node_id=%s", node_id))
    tree <- readSimulationData(path)
    node_id <-  as.character(node_id)
    node <- tree[[as.character(node_id)]]

    fit <- samplePoints(node, nextra=nextra)
    if (is.null(fit)) {
        message(sprintf("cannot sample points for path %s, node %s", path, node_id))
        return(NULL)
    } 

    if (!is.null(model))
        node$model <-  model
    node$fit = fit

    title <-  sprintf("%s #%s",  path, node_id)

    return(plotnode(node, node_id, title=title, zoom=zoom))

}

##'
##' plot data extracted from the likelihood table for a particulat node.
##'
##' .. content for \details{} ..
##' @param path 
##' @param node_id 
##' @param nextra 
##' @param model 
##' @param zoom 
##' @return 
##' @author chris
plotnode <- function(node, node_id, zoom=F, title="") {

    fit <- node$fit
    node$bls$fit_ll <- calculate_lcfit(node$bls, node$model, fit, weighted=F)$ll
    node$bls$wfit_ll <- calculate_lcfit(node$bls, node$model, fit, weighted=T)$ll
    node$bls$t4fit_ll <- calculate_lcfit(node$bls, node$model, fit, weighted=F, keep=4)$ll
    bls <- node$bls

    # bls <- bls[,c('node_id','branch_length', 'bpp_ll', 'fit_ll')]
    bls <- melt(bls, id.vars=1:2)

    if (!zoom) {
        # full plot
        p <- ggplot() + ylab("Log likelihood")
        p <- p + geom_line(aes(x=branch_length, y=value, color=variable, linetype=variable),
                          data=bls) +
                ggtitle(title) +
                theme(plot.title=element_text(size=15)) +
                geom_point(aes(x=branch_length, y=ll), data=fit) +
                scale_color_discrete(name="", breaks=c('bpp_ll', 'fit_ll', 'wfit_ll', 't4fit_ll'), labels=c('bpp', 'unweighted', 'weighted', 'top4')) +
                guides(color="legend", shape=FALSE, linetype=FALSE) +
                theme(legend.position="bottom") +
                ylab("Log likelihood") 
                #annotate(geom="text", label="foo", x=.5, y=.9)
    } else {
        # zoomed-in partial plot
        # find the domain of the top 4 points
        tmp <- fit[order(fit$ll, decreasing=T)[1:4],]
        r <- extendrange(tmp$branch_length, f=2)
        p <- ggplot() + ylab("Log likelihood")
        p <- p + geom_line(aes(x=branch_length, y=value, color=variable, linetype=variable),
                          data=subset(bls, branch_length >= r[[1]] & branch_length <= r[[2]] )) +
                ggtitle(title) +
                theme(plot.title=element_text(size=15)) +
                    # geom_point(aes(x=branch_length, y=value), data=subset(bls, branch_length >= r[[1]] & branch_length <= r[[2]] )) +

                geom_point(aes(x=branch_length, y=ll), data=tmp) +
                scale_color_discrete(name="", breaks=c('bpp_ll', 'fit_ll', 'wfit_ll', 't4fit_ll'), labels=c('bpp', 'unweighted', 'weighted', 'top4')) +
                guides(color="legend", shape=FALSE, linetype=FALSE) +
                theme(legend.position="bottom") +
                ylab("Log likelihood") 
                #annotate(geom="text", label="foo", x=.5, y=.9)
    }
    return(p)
}
