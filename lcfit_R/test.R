
library(ggplot2)
library(plyr)
library(reshape2)
suppressPackageStartupMessages(library(distrEx))
suppressPackageStartupMessages(library(entropy))
source("../compare/utils.R")

theme_set(theme_bw(16))

library(lcfit)
## library(devtools)
## load_all(recompile=TRUE)
## dyn.load("src/lcfit.so")
## is.loaded('lcfit_rcpp_fit_bsm', PACKAGE = 'lcfit')

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
    stopifnot("x" %in% names(df))
    stopifnot("y" %in% names(df))
    stopifnot(is.list(model) || is(model, "numeric"))

    p <- pts[pts$y == max(pts$y),]
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


lcfit.root = "../"

sim.dir <- paste0(lcfit.root, "sims/runs/10/3/Binary-0.25/gamma4-0.2/")

name.bls <- paste0(sim.dir, "lcfit_bls.csv")
name.maxima <- paste0(sim.dir, "lcfit_maxima.csv")
name.fit <-  paste0(sim.dir, "lcfit_fit.csv")

data.bls <- read.csv(name.bls, as.is=TRUE)
data.maxima <- read.csv(name.maxima, as.is=TRUE)
data.fit <- read.csv(name.fit, as.is=TRUE)

nodes <- unique(data.fit$node_id)
node_id <- nodes[[1]]
print(sprintf("node_id = %d", node_id))
ndata.bls <- data.bls[data.bls$node_id == node_id,]
ndata.maxima <- data.maxima[data.maxima$node_id == node_id,]
ndata.fit <- data.fit[data.fit$node_id == node_id,]
ndata.fit$name <- 'Bio++'
# calculate summary statscomparing lcfit to likelihood
error.lcfit <- compare_bls(ndata.bls)

# calculate RSS between likelihood and lcfit
rss.lcfit <- with(ndata.bls, as.integer(sum((fit_ll - bpp_ll)^2)))


ndata.bls <- melt(ndata.bls, id.vars=1:2)

    p <- ggplot( ndata.bls) +
         geom_line(aes(x=branch_length, y=value, color=variable, linetype=variable), data=ndata.bls) +
         ggtitle(bquote(atop(.(sprintf("Node #%s",  ndata.bls$node_id[1])),
                             atop(scriptscriptstyle(italic(RSS[lcfit] ~ "=" ~ .(rss.lcfit))),
                                  phantom(0))
                             ~
                             atop(scriptscriptstyle(italic(KL[lcfit] ~ "=" ~ .(sprintf("%.3f",error.lcfit$kl)))),
                                  phantom(0))

                             ))) +
         geom_point(aes(x=branch_length, y=ll, shape=name), data=ndata.fit) +
         xlim(0, max(c(max(ndata.fit$branch_length), 1)))
    p <- p + scale_color_discrete(name="", breaks=c('bpp_ll', 'fit_ll'), labels=c('Bio++', 'lcfit')) 
    p <- p + guides(color="legend", shape=FALSE, linetype=FALSE) 
    p <- p + theme(legend.position="bottom")
print(p)

model <-  ndata.maxima[, c('c', 'm', 'r', 'b')]


DEFAULT_MODEL <- c(c=1800.0, m=400.0, r=1.0, b=0.5)
model <- DEFAULT_MODEL
df <- ndata.fit[, c('branch_length', 'll')]
names(df) <- c('x', 'y')
model <- fit_model(model, df)

y <- sapply(ndata.bls[['branch_length']], lcfit, model=model)

p <-  p + geom_line(aes(x=x,y=y), data=data.frame(x=ndata.bls[['branch_length']], y=y))
print(p)



#bsm <- c(c=4461.25, m=2764.04, r=7.94896, b=0.125468)


pts <- c(0.000151781,-2144.03,
         0.00151781,-2133.88,
         0.0151781,-2135.72,
         0.151781,-2232.7,
         0.20227,-2261.73,
         0.4648,-2361.22)
pts <- matrix(pts, ncol=2, byrow=T)
pts <- data.frame(x=pts[,1], y=pts[,1]) 


model.fit <- lcfit_fit_bsm(pts$x, pts$y, DEFAULT_MODEL)


lcfit <- function(model, t) {
    env <- list2env(lapply(model, function(x) x))
    with(env, {
        e <- exp(-r*(t+b))
        l <- c * log((1.0 + e)/2.0) + m * log((1.0 - e)/2.0)
    })
}



