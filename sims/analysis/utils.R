library(lcfit)

# DEFAULT_MODEL=c(c=1500, m=1000,r=2.0,b=0.5)  # not sure where from - probably from inside lcfit.c
DEFAULT_MODEL=c(c=1100, m=800,r=2.0,b=0.5)   # from lcfit/sims/Sconstruct:initial_values
# DEFAULT_MODEL=c(c=1228, m=700,r=2.0,b=0.5)   # adjusted based on empirical observations.

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


ess <- function(ll) { ll <- ll - max(ll); exp(-log(sum(exp(2 * ll))) + 2 * log(sum(exp(ll)))) }

logsumexp <- function(x) {
  log(sum(exp(x - max(x)))) + max(x)
}

kl_log <- function(log_freqs1, log_freqs2) {
  log_freqs1 <- log_freqs1 - logsumexp(log_freqs1)
  log_freqs2 <- log_freqs2 - logsumexp(log_freqs2)
  LR <- log_freqs1 - log_freqs2
  sum(exp(log_freqs1) * LR) / log(2)
}


##' Calculates various error measurements between a fitted estimate of likelihood curve and the true likelihood curve.
##'
##' This is a poor modification of what Connor wrote.  I think his original only compared a particular known estimate with
##' the true likelihood.   I modified it to take the name of the column holding the estimate, but really this should
##' take just the true bpp_ll and the estimate.  I didn't change it yet b/c I would have to go find all th places that
##' depend on the current signature.
##'
##' There is another functionof the same name in aggutils.R that looks more like what I hope this would be!
##' 
##' @param d data frame with columns 'branch_length', 'bpp_ll', and one other column holding a likelihood estimate, usually 'fit_ll'
##' @param fit the ame of the column holding the likelihood estimate
##' @return a data frame holding various estimates of deviation between the tree likelihood in 'bpp_ll' and a likelihood estimate.
deprecated_compare_bls <- function(d, fit='fit_ll') {
    ## deprecated to find the places where this is called.  It is time to fix up those calls!

  branch_length <- d[['branch_length']]
  bpp_ll <- d[['bpp_ll']]
  fit_ll <- d[[fit]]

  bpp_l <- exp(bpp_ll - max(bpp_ll))
  bpp_l <- bpp_l / sum(bpp_l)
  fit_l <- exp(fit_ll - max(fit_ll))
  fit_l <- fit_l / sum(fit_l)

  bpp_d <- DiscreteDistribution(seq_along(branch_length), prob = bpp_l)
  fit_d <- DiscreteDistribution(seq_along(branch_length), prob = fit_l)
  data.frame(hellinger=HellingerDist(bpp_d, fit_d),
             kl=KL.plugin(bpp_l, fit_l, unit = 'log2'),
             my_kl=kl_log(bpp_ll, fit_ll),
             bpp_rel_ess=ess(bpp_ll) / length(bpp_ll),
             fit_rel_ess=ess(fit_ll) / length(bpp_ll))
}

## https://gist.github.com/doobwa/941125
log.sum.exp<- function(x) {
  # Computes log(sum(exp(x))
  # Uses offset trick to avoid numeric overflow: http://jblevins.org/notes/log-sum-exp
  if ( max(abs(x)) > max(x) )
    offset <- min(x)
  else
    offset <- max(x)
  log(sum(exp(x - offset))) + offset
}

##' Calculate exponential weights for sample points according to log-likelihood.
##'
##' These values are used to weight the minimization function when
##' fitting a nonlinear function to the sample points.  The weights are
##' distributed to favor points close to the maximum likelihood.
##'
##' This essentially turns the log-likelihood values into a discrete probability distribution.
##' Returns a numeric vector of weights corresponding to the vector of
##' likelihood values.
##' 
##' @param ll numeric vector of log-likelihood values
##' @return numeric vector of weights.
##' @author chris
fit_weights <- function(ll) {
    # 
    lambda <- 2
    weights <- ll * (1/lambda)
    weights <- weights - max(weights)
    sum = log.sum.exp(weights)
    weights <- weights - sum
    weights <- exp(weights)
    return(weights)
}

##' Default routine for calculating unitary weights across all sample points.
##'
##' This is the weight function to use if you want to do unweighted nonlinear approximation.
##' @param ll  nuemric vector of log-likelihood values.
##' @return numeric vector of 1.0, same length as ll
##' @author chris
fit_equal_weights <- function(ll) {
    weights <-  rep(1.0, length(ll))
    return(weights)
}


##' Fit lcfit model parameters to sample points
##'
##' Given base model parameters, adjust the parameters to include the sample points.
##' @title 
##' @param model	named numeric vector, { 'c', 'm', 'r', 'b' }
##' @param pts 		sampled points; dataframe { 'x', 'y' }
##' @return 		named numeric vector, { 'c', 'm', 'r', 'b' }
##' @author chris
fit_model <- function(model, pts, weighted=F, max.iter = 250) {
    stopifnot(is(pts, "data.frame"))
    stopifnot(all(c('x', 'y') %in% names(pts)))
    stopifnot(is.list(model) || is(model, "numeric"))
    stopifnot(all(c('c', 'm', 'r', 'b') %in% names(model)))

    # scale the model to the largest sample point.
    # What's going on here?
    #
    # lcfit.c::lcfit_bsm_scale_factor() calculates a scaling factor to keep the resulting log-likelihood values
    # within bounds.  It takes a single point as an argument, which represents the greatest log-likelihood value.
    # The easiest way to obtain this point is to take the maximum y value, but if the maximum point is repeated in the list,
    # then you will be get two points.  We warn against this and avoid it by using only a single max value.
    # 
    p <- pts[pts$y == min(pts$y),]
    p <- p[1,]	# make sure to use only a single point

    scale_factor <- lcfit_bsm_scale_factor(p$x, p$y, model);
    model['c'] <- model['c'] * scale_factor;
    model['m'] <- model['m'] * scale_factor;

    # unweighted
    if (weighted)
        weights <- fit_weights(pts$y)
    else
        weights <- fit_equal_weights(pts$y)

    model = lcfit_fit_bsm(pts$x, pts$y, weights, model, max.iter);
    return(model)
}


##' Apply the lcfit function to a set of branch lengths
##'
##' Give a set of branch lengths and model parameters, calculate the approximate likelihood values.
##' @param t  branch lengths
##' @param model named numeric vector with model parameters 'b', 'r', 'm', and 'c'
##' @return numeric vector (same length as t) of approximate likelihood values.
##' @author chris
lcfit <- function(t, model) {
    env <- list2env(lapply(model, function(x) x))
    with(env, {
        e <- exp(-r*(t+b))
        l <- c * log((1.0 + e)/2.0) + m * log((1.0 - e)/2.0)
    })
}


##' Sample additional points from a likelihood curve.
##'
##' .. content for \details{} ..
##' @param node_id 
##' @param nextra 
##' @param keep 
##' @return 
##' @author chris
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

            # too many rejections result in failure.
            nrejected <- nrejected + 1
            if (nrejected > 200) {
                extra.fit <- NULL
                break
            }
        }

        if (is.null(extra.fit))
            # we didn't successfully generate a good set of extra points.
            # fail out of the routine.
            # should raise an exception here.
            fit <- NULL
        else
            # succcess!
            fit <- rbind(node$fit, extra.fit)
    }

    if (!is.null(fit)) {
        # double-super-paranoid-extra check for duplicated branch lengths.
        # can almost certainly, probably, possibly, get rid of this.
        if (sum(fit$branch_length == max(fit$branch_length)) != 1) {
            message("in samplePoints: Saw duplicated branch_lengths")
        }
    }
    return(fit)
}

##' re-fit the lcfit model to the new sample points.
##'
##' .. content for \details{} ..
##' @param ndata 
##' @param weighted 
##' @param keep 
##' @return 
##' @author chris
calculate_lcfit <- function(ndata, weighted=F, keep=0, max.iter=250) {
    stop("stop calling deprecated calculate_lcfit() function in analysis/utils.R")
    # replaced with explicit call to fit_model() followed by sapply(lcfit) to apply the model to branch lengths.
    # model <- fit_model(model, sample.points, weighted, max.iter)
    # sapply(ndata$bls[['branch_length']], lcfit, model=model)
}


calculate_spline <- function(ndata) {
    # fit a spline to the sample points.  Interpolate at the branch_point values used for lcfit_ll.
    # ndata$bls$spline_ll <- 
    spline(ndata$fit$branch_length, y=ndata$fit$ll,
                                  xout=ndata$bls[['branch_length']], method = "natural")[[2]]

}

##' plot a spline approximation for the likelihood curve. 
##'
##' .. content for \details{} ..
##' @title 
##' @param node_id 	name of node to operaton on
##' @param nextra 	number of extra points to select
##' @return Returns nothing. Sends plot output to current graphics device.
##' @author chris
plot_node<- function(ndata) {

    bls <- melt(ndata$bls, id.vars=1:2)

    p <- ggplot( bls, aes(color=name, linetype=name)) +
         geom_line(aes(x=branch_length, y=value, color=variable, linetype=variable), data=bls) +
         ggtitle(sprintf("Node #%s",  bls$node_id[1])) +
         geom_point(aes(x=branch_length, y=ll, shape=name), data=ndata$fit) +
         ylab("Log likelihood")

    p <- p + scale_color_discrete(name="", breaks=c('fit_ll', 'wfit_ll', 't4fit_ll'), labels=c('unweighted', 'weighted', 'top4')) 
    p <- p + guides(color="legend", shape=FALSE, linetype=FALSE) 
    p <- p + theme(legend.position="bottom")

}



summaryTable <- function(ndata, compare.to='wfit_ll') {
    error.unweighted <- compare_bls(ndata$bls, fit='fit_ll')
    error.weighted <- compare_bls(ndata$bls, fit='wfit_ll')
    error.topfour <- compare_bls(ndata$bls, fit='t4fit_ll')

    w <- fit_weights(ndata$bls$bpp_ll)

    # calculate summary statistics - the residual squared difference between
    # the fitting procedures (lcfit and spline) and the actual likelihood.
    rss.unweighted <- sum((w*(ndata$bls$fit_ll- ndata$bls$bpp_ll))^2)
    rss.weighted <- sum((w*(ndata$bls$wfit_ll- ndata$bls$bpp_ll))^2)
    rss.topfour <- sum((w*(ndata$bls$t4fit_ll- ndata$bls$bpp_ll))^2)

    tbl <- data.frame(c(rss.unweighted, error.unweighted$kl),
                      c(rss.weighted, error.weighted$kl),
                      c(rss.topfour, error.topfour$kl)
                      )
    names(tbl) <- c('unweighted', "weighted", "topfour")
    rownames(tbl) <- c("RSS", "KL")
    return(tbl)
}




showSummary <- function(ndata) {
    tbl <- do.call(cbind, lapply(ndata, summaryTable, compare.to="wfit_ll"))
    rownames(tbl) <- c("RSS", "KL divergence")
    knitr::kable(tbl)
    }



showplots <- function(nodes, nextra) {
    set.seed(1234)
    ndata <- lapply(nodes,
                function(node) {
		    nd <- samplePoints(node, nextra=nextra)
                    nd$bls$fit_ll <- calculate_lcfit(nd, weighted=F)
                    nd$bls$wfit_ll <- calculate_lcfit(nd, weighted=T)
                    nd$bls$t4fit_ll <- calculate_lcfit(nd, weighted=F, keep=4)
                    nd })


    plots <- lapply(ndata, function(nd) {plot_node(nd) })
    multiplot(plotlist=plots, cols=3)
    return(ndata)
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
