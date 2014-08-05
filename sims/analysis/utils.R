library(lcfit)


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

# takes a data frame with columns 'branch_length', 'bpp_ll', and 'fit_ll'

compare_bls <- function(d, fit='fit_ll') {
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

plotPoly <- function(node_id, degree=3, nextra = 0) {
    ndata.bls <- data.bls[data.bls$node_id == node_id,]
    ndata.maxima <- data.maxima[data.maxima$node_id == node_id,]
    ndata.fit <- data.fit[data.fit$node_id == node_id,]
    ndata.fit$name <- 'Bio++'

    # sample some extra points
    if (nextra > 0) {
        i <- sample(1:nrow(ndata.bls), nextra, replace = FALSE)
        extra.fit <- data.frame(node_id=node_id, branch_length=ndata.bls[i,"branch_length"], ll=ndata.bls[i,"bpp_ll"], name="extra")
        ndata.fit <- rbind(ndata.fit, extra.fit)
    }
    
    # fit a polynomial to the sample points.  Interpolate over branch_length.
    fit <- lm(ll ~ poly(branch_length,degree,raw=TRUE), data=ndata.fit)
    
    ndata.bls$poly_ll <- predict(fit, data.frame(branch_length=ndata.bls[['branch_length']]))

    # calculate summary statistics
    error.lcfit <- compare_bls(ndata.bls)
    error.poly <- compare_bls(ndata.bls, fit='poly_ll')
    
    # calculate summary statistics - the residual squared difference between
    # the fitting procedures (lcfit and spline) and the actual likelihood.
    rss.poly <- with(ndata.bls, as.integer(sum((poly_ll - bpp_ll)^2)))
    rss.lcfit <- with(ndata.bls, as.integer(sum((fit_ll - bpp_ll)^2)))

    ndata.bls <- melt(ndata.bls, id.vars=1:2)

    p <- ggplot( ndata.bls) +
         geom_line(aes(x=branch_length, y=value, color=variable, linetype=variable), data=ndata.bls) +
         ggtitle(bquote(atop(.(sprintf("Node #%s",  ndata.bls$node_id[1])),
                             atop(scriptscriptstyle(italic(RSS[lcfit] ~ "=" ~ .(rss.lcfit))),
                                  scriptscriptstyle(italic(RSS[poly] ~ "=" ~ .(rss.poly))))
                             ~
                             atop(scriptscriptstyle(italic(KL[lcfit] ~ "=" ~ .(sprintf("%.3f",error.lcfit$kl)))),
                                  scriptscriptstyle(italic(KL[poly] ~ "=" ~ .(sprintf("%.3f",error.poly$kl)))))

                             ))) +
         geom_point(aes(x=branch_length, y=ll, shape=name), data=ndata.fit) +
         xlim(0, max(c(max(ndata.fit$branch_length), 1)))
    p <- p + scale_color_discrete(name="", breaks=c('bpp_ll', 'fit_ll', 'poly_ll'), labels=c('Bio++', 'lcfit', 'poly')) 
    p <- p + guides(color=FALSE, shape=FALSE, linetype=FALSE) 
    p <- p + guides(color="legend")
    p <- p + theme(legend.position="bottom")

}

##' Fit lcfit model parameters to sample points
##'
##' Given base model parameters, adjust the parameters to include the sample points.
##' @title 
##' @param model	named numeric vector, { 'c', 'm', 'r', 'b' }
##' @param pts 		sampled points; dataframe { 'x', 'y' }
##' @return 		named numeric vector, { 'c', 'm', 'r', 'b' }
##' @author chris
fit_model <- function(model, pts) {
    stopifnot(is(pts, "data.frame"))
    stopifnot(all(c("x", "y") %in% names(pts)))
    stopifnot(is.list(model) || is(model, "numeric"))
    
    # scale the model to the largest sample point
    # THIS IS WHERE THE ERROR OCCURS!
    p <- pts[pts$y == max(pts$y),]
    if (nrow(p) != 1)
        warning("Fitting duplicated points?  Saw multiple maximum values.")
    p <- p[1,]
    scale_factor <- lcfit_bsm_scale_factor(p$x, p$y, model);
    model['c'] <- model['c'] * scale_factor;
    model['m'] <- model['m'] * scale_factor;

    model = lcfit_fit_bsm(pts$x, pts$y, model);
    return(model)
}

lcfit <- function(t, model) {
    env <- list2env(lapply(model, function(x) x))
    with(env, {
        e <- exp(-r*(t+b))
        l <- c * log((1.0 + e)/2.0) + m * log((1.0 - e)/2.0)
    })
}

##' plot a spline approximation for the likelihood curve. 
##'
##' .. content for \details{} ..
##' @title 
##' @param node_id 	name of node to operaton on
##' @param nextra 	number of extra points to select
##' @return Returns nothing. Sends plot output to current graphics device.
##' @author chris
plotSpline <- function(node_id, nextra=0, keep=0) {
    print(sprintf("node_id = %d", node_id))
    ndata.bls <- data.bls[data.bls$node_id == node_id,]
    ndata.maxima <- data.maxima[data.maxima$node_id == node_id,]
    ndata.fit <- data.fit[data.fit$node_id == node_id,]
    ndata.fit$name <- 'Bio++'

    # sample some extra points
    if (nextra > 0) {
        # compute the categorical distribution of points
        bpp_ll <- ndata.bls[['bpp_ll']]
        bpp_l <- exp(bpp_ll - max(bpp_ll))
        bpp_l <- bpp_l / sum(bpp_l)

        # propose some more sample points, but reject the proposal if any of the new points duplicate existing points.
        repeat {
            # sample w/o replacement according to the distribution calculated above
            i <- sample(1:nrow(ndata.bls), nextra, replace = FALSE, prob=bpp_l)
            # add the new points to the existing ones.
            extra.fit <- data.frame(node_id=node_id, branch_length=ndata.bls[i,"branch_length"], ll=ndata.bls[i,"bpp_ll"], name="extra")
            if (!any(sapply(extra.fit$ll, function(x) x == ndata.fit$ll ))) {
                break
            }
            else {
                message("proposal rejected")
            }
        }

        ndata.fit <- rbind(ndata.fit, extra.fit)

        # Erick asks to keep only the four-highest likelihood points.
        if (keep > 0) {
            ndata.fit <- ndata.fit[order(ndata.fit$ll, decreasing=T)[1:keep],]
        }
    }

    # re-fit the lcfit model to the new sample points.
    model <-  ndata.maxima[, c('c', 'm', 'r', 'b')]
    sample.points <- ndata.fit[, c('branch_length', 'll')]
    names(sample.points) <- c('x', 'y')
    model <- fit_model(model, sample.points)
    ndata.bls$fit_ll <- sapply(ndata.bls[['branch_length']], lcfit, model=model)

    # fit a spline to the sample points.  Interpolate at the branch_point values used for lcfit_ll.
    ndata.bls$spline_ll <- spline(ndata.fit$branch_length, y=ndata.fit$ll,
                                  xout=ndata.bls[['branch_length']], method = "natural")[[2]]

    error.lcfit <- compare_bls(ndata.bls)
    error.spline <- compare_bls(ndata.bls, fit='spline_ll')
    

    # calculate summary statistics - the residual squared difference between
    # the fitting procedures (lcfit and spline) and the actual likelihood.
    rss.spline <- with(ndata.bls, as.integer(sum((spline_ll - bpp_ll)^2)))
    rss.lcfit <- with(ndata.bls, as.integer(sum((fit_ll - bpp_ll)^2)))

    ndata.bls <- melt(ndata.bls, id.vars=1:2)

    p <- ggplot( ndata.bls, aes(color=name, linetype=name)) +
         geom_line(aes(x=branch_length, y=value, color=variable, linetype=variable), data=ndata.bls) +
         ggtitle(bquote(atop(.(sprintf("Node #%s",  ndata.bls$node_id[1])),
                             atop(scriptscriptstyle(italic(RSS[lcfit] ~ "=" ~ .(rss.lcfit))),
                                  scriptscriptstyle(italic(RSS[spline] ~ "=" ~ .(rss.spline))))
                             ~
                             atop(scriptscriptstyle(italic(KL[lcfit] ~ "=" ~ .(sprintf("%.3f",error.lcfit$kl)))),
                                  scriptscriptstyle(italic(KL[spline] ~ "=" ~ .(sprintf("%.3f",error.spline$kl)))))

                             ))) +
         geom_point(aes(x=branch_length, y=ll, shape=name), data=ndata.fit) +
         xlim(0, max(c(max(ndata.fit$branch_length), 1))) 

    p <- p + scale_color_discrete(name="", breaks=c('bpp_ll', 'fit_ll', 'spline_ll'), labels=c('Bio++', 'lcfit', 'spline')) 
    p <- p + guides(color="legend", shape=FALSE, linetype=FALSE) 
    p <- p + theme(legend.position="bottom")
}
