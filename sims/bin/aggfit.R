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

DEFAULT_MODEL=c(c=1100,m=800,r=2.0,b=0.5)

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

source("analysis/aggutils.R")

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
    errs <- ldply(c(weighted="weighted", unweighted="unweighted", top4="top4"),
                 function(analysis) {
                     ll <- switch(analysis,
                                  weighted = calculate_lcfit(bls=node$bls, model=DEFAULT_MODEL, fit=node$fit, weighted=T),
                                  unweighted = calculate_lcfit(bls=node$bls, model=DEFAULT_MODEL, fit=node$fit, weighted=F),
                                  top4 = calculate_lcfit(bls=node$bls, model=DEFAULT_MODEL, fit=node$fit, weighted=F, keep=4)
                                  )

                     # If lcfit failed, set all error estimates to the negated value of the error status.
                     # negetive error estimates indicate a failure, and the particular value indicates
                     # the type of failure.
                     # See the definition of enum lcfit_status in lcfit_R/src/lcfit.h for error codes.
                     #
                     df <- compare_bls(ll$ll, node$bls$bpp_ll)
                     df$status <- ll$status
                     return(df)

                 })
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
aggregate <- function(control_paths, model=NULL) {
    if (!is.null(model)) {
        message(paste0("Force default model ", paste0(model, collapse=",")))
    }
    # for each tree...
    ldply(setNames(nm=control_paths), .id="tree", 
          function(path) {
              message(path)
              # read control file, then read in various data files
              # pointed to by the control file.
              tree <- readSimulationData(path)
              names(tree) <- sapply(tree, function(node) node$model$node_id)

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
                                         rss <- rss[,-c(1), drop=FALSE]
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


control_paths <- args
aggdata <- aggregate(control_paths, model=DEFAULT_MODEL)
write.csv(aggdata, row.names=FALSE)

