#!/usr/bin/env Rscript

##' Aggregate error mesures for various types of fittig procedures across all the simulated trees.
##' 
##' usage: aggfit.R <control_file> ...
##' e.g.   aggfit.R runs/10/0/JTT92/gamma4-0.2/control.json runs/10/0/JTT92/uniform/control.json
##'
##' Usually you would run this 
##'	$ bin/aggfit.R `find runs -name 'control.json' | head -2`


library(ggplot2)
library(gridExtra)
library(plyr)
library(reshape2)
library(RJSONIO)
source("analysis/utils.R")
# source("analysis/aggutils.R")


theme_set(theme_bw(16) + theme(strip.background=element_blank(),
                               legend.key = element_blank(),
                               legend.position='bottom'))


##' Calculate error estimates various fitting techniques, and return status of numerical fitting procedure
##'
##' Returns a data frame that looks like this,
##'                   .id         wrss status
##'          1   weighted 1.318217e-05      0
##'          2 unweighted 1.607818e-02      0
##'          3       top4 7.852516e-05      1
##'
##' This has the weighted rss error estimate between three fitting procedures and the true liklihood,
##' along with the status from each of the fitting procedures.  A status of 0 indicates success, otherwise an error.
##' Error status values can be found in lcfit_R/src/lcfit.h
##' 
##' @param node named list containing bls, fit, and maxima parameters
##' @return dataframe of error mesures for each of the types of fitting techniques
##' @author chris
##' 
calculate.model.errors <- function(node) {

    calculate.helper <- function(analysis) {
        # prepare sample points
        pts <- node$fit[, c('branch_length', 'll')]
        names(pts) <- c('x', 'y')
        
        model <- switch(analysis,
                        weighted = fit_model(model=DEFAULT_MODEL, pts, weighted=T),
                        unweighted = fit_model(model=DEFAULT_MODEL, pts, weighted=F),
                        top4 = {
                            # Erick asks to keep only the four-highest likelihood points.
                            pts <- pts[order(pts$y, decreasing=T)[1:4],]
                            fit_model(model=DEFAULT_MODEL, pts, weighted=F)
                        }
                        )
        if (model['b'] < 0) {
            message("Warning: b < 0")
        }
        ll <- sapply(node$bls[['branch_length']], lcfit, model=model)
        if (any(is.na(ll))) {
            message(sprintf("NA found in %s likelhoods", analysis))
            message(sprintf("c=%f m=%f r=%f b=%f", model['c'], model['m'], model['r'], model['b']))
        }
        
        # See the definition of enum lcfit_status in lcfit_R/src/lcfit.h for error codes.
        #
        df <- compare_bls(ll, node$bls$bpp_ll)
        df$status <- model[['status']]
        return(df)
    }

    errs <- ldply(c(weighted="weighted",
                    unweighted="unweighted",
                    top4="top4"), calculate.helper)
    return(errs)

}
    
# for each node
#     for N extra points
#         for unweighted, weighted, top4
#             calculate fit
#             compare to empirical
#
# syntax for naming the .id column taken from,
# http://stackoverflow.com/a/15559709/1135316
aggregate <- function(control_paths, model) {
    if (!missing(model)) {
        message(paste0("Force default model ", paste0(model, collapse=",")))
    }
    # for each tree...
    ldply(setNames(nm=control_paths), .id="tree", 
          function(path) {
              # read control file, then read in various data files
              # pointed to by the control file.
              tree <- readSimulationData(path)
              names(tree) <- sapply(tree, function(node) node$model$node_id)
              message(path)
              # for each node in the tree...
              ldply(tree, .id="node_id", .parallel=FALSE,
                    function(node) {

                        # estimate log-likelihood with and without extra sample points
                        # for increasing numbers of additional sample points....
                        tmp <- llply(0:4,
                                     function(nextra) {
                                         fit <- samplePoints(node, nextra=nextra)
                                         if (is.null(fit)) {
                                             message(sprintf("cannot sample points for path %s, node %d", path, node$bls$node_id[[1]]))
                                         } else
                                             node$fit <- fit
                                         # estimate log-likelihood using a variety of methods,
                                         # and return an error between the estimate and empircal log-likelihood
                                         
                                         # use model as provided, otherwise us maxima model as determined from simulation
                                         if (!is.null(model)) node$model <- model

                                         rss <- calculate.model.errors(node)
                                         
                                         rownames(rss) <- paste0(rss$.id, "_", nextra)
                                         rss <- rss[,-match('.id', names(rss)), drop=FALSE]   # drop the '.id' column
                                         # melting a matrix gives access to rownames
                                         # http://stackoverflow.com/a/19826856/1135316
                                         rss.m <- melt(t(rss))
                                         rss <- dcast(rss.m, . ~ ...)
                                         rss <- rss[,-c(1),drop=FALSE]
                                         
                                     })
                        tmp <- do.call(cbind, c(tmp))
                    })
          })
}

library(doMC)
registerDoMC(8)
# DEFAULT_MODEL=c(c=1500, m=1000,r=2.0,b=0.5)  # not sure where from - probably from inside lcfit.c
DEFAULT_MODEL=c(c=1100, m=800,r=2.0,b=0.5)   # from lcfit/sims/Sconstruct:initial_values


args <- commandArgs(TRUE)
stopifnot(length(args) >= 1)
control_paths <- args
# control_paths <- c("runs/10/0/JTT92/gamma4-0.2/control.json","runs/10/0/JTT92/uniform/control.json")
# control_paths <- c("runs/100/9/JC/uniform/control.json")

aggdata <- aggregate(control_paths, model=DEFAULT_MODEL)

# The aggregated data should have no NA values!
stopifnot(!is.na(aggdata))

write.csv(aggdata, row.names=FALSE)

