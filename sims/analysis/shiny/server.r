library(shiny)
library(ggplot2)
library(grid)
library(plyr)
library(reshape2)
library(RJSONIO)
source("../utils.R")

theme_set(theme_bw(16))

source("../aggutils.R")

DEFAULT_MODEL=c(c=1500,m=1000,r=2.0,b=0.5) # taken from sim/bin/aggfit.R:DEFAULT_MODEL
DEFAULT_MODEL=c(c=1100,m=800,r=2.0,b=0.5)	# from lcfit/sims/Sconstruct:initial_values

m = NULL

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {

    # Implement the 'reset' button.
    # http://stackoverflow.com/a/24269691/1135316
    output$model_input <- renderUI({
        message("render model ui")
        times <- input$reset_input
        div(id=letters[(times %% length(letters)) + 1],
            sliderInput("model.c", "C", 800, 2500, DEFAULT_MODEL['c'], step=2),
            sliderInput("model.m", "M", 700, 1400, DEFAULT_MODEL['m']),
            sliderInput("model.r", "R", 0, 5, DEFAULT_MODEL['r'], step=0.1),
            sliderInput("model.b", "B", 0, 2, DEFAULT_MODEL['b'], step=0.1)
            )
    })


    
    
    ## load data for all the nodes in all the trees.
    loadAggData <- reactive({
        aggfile <- "../../aggfit.csv"
        aggdata <- read.csv(aggfile, stringsAsFactors=F)
        message(sprintf("loaded %d rows from %s", nrow(aggdata), aggfile))

        m.measures <- melt(aggdata, id.vars=c("tree", "node_id"), measure.vars=names(aggdata)[!grepl('tree|node_id|status_', names(aggdata))])
        m.measures <- cbind(m.measures, colsplit(m.measures$variable, names = c("measure", "fitting", "nextra"), pattern="_"))
        m.measures = m.measures[,-grep('variable',names(m.measures))]
        
        m.status <- melt(aggdata, id.vars=c("tree", "node_id"), measure.vars=names(aggdata)[grepl('status_', names(aggdata))], value.name='status')
        m.status <- cbind(m.status, colsplit(m.status$variable, names = c("measure", "fitting", "nextra"), pattern="_"))
        m.status = m.status[,-grep('variable|measure',names(m.status))]

        print(sprintf("len of measures table = %d", nrow(m.measures)))
        print(sprintf("len of status table = %d", nrow(m.status)))
        message(paste0("measures names = ", paste0(names(m.measures), collapse=',')))
        message(paste0("status names = ", paste0(names(m.status), collapse=',')))

        m <- join(m.measures, m.status, by=c("tree", "node_id", "fitting", "nextra"), type="left")
        print(sprintf("len of merged  table = %d", nrow(m)))

        # it seems there is no easy way to avoid having melt convert your character columns to factors.
        # convert them back to characters
        # wut?  When I check this it does not seem true.
        # > sapply(m, class)
        #         tree     node_id    variable       value 
        #  "character"   "integer"    "factor"   "numeric"
        # m$tree <- as.character(m$tree)
        m <- m[!is.na(m$value),]
        return(m)
    })


    output$caption <- renderText({
        m <- loadAggData()
        ntotal <- nrow(m)
        nfailures <- sum(m$status != 0)
        nselected <- nrow(selectedNodes())

        sprintf("%d selected out of %d failures, %d total samples.", nselected, nfailures, ntotal)
    })

    selectedNodes <- reactive({
        m <- loadAggData()
        
        ## let the user select which failure mode they are interested in.
        ## LCFIT_MAXITER = 1,	// exceeded maximum iterations without converging
        ## LCFIT_ERROR = 2,	// non-specific error
        ## LCFIT_ENOPROG = 27,	// iteration is not making progress towards solution
        ## LCFIT_ETOLF = 29	// cannot reach the specified tolerance in F
        errs = c(MAXITER=1, ERROR=2, ENOPROG=27, ETOLF=29)
        errs = errs[input$failures]
        message(sprintf("nrows of m = %d", nrow(m)))
        # select some nodes that meet our criteria
        tmp <- m[m$measure=="wrss" & m$nextra == input$nextra & m$fitting == input$fitting & m$status %in% errs,]
        if (nrow(tmp) == 0) 
            return()
        tmp <- tmp[order(tmp$value, decreasing=F),]
        return(tmp)

    })

        
    nodeData <- reactive({
        tmp <- selectedNodes()
            
        # choose one of the nodes at random
        n <- sample(1:nrow(tmp), 1)
        message(sprintf("selecting row %d out of %d", n, nrow(tmp)))
        print(tmp[n,])

        # read in the sampled points for that node in that particular tree...
        wd = getwd()
        tryCatch( {
            setwd("../../")
            tree <- readSimulationData(tmp[n,'tree'])
            message(sprintf("loaded tree data from %s", tmp[1,'tree']))

        }, finally={ setwd(wd) })
        node_id <- as.character(tmp[n, 'node_id'])
        node <- tree[[as.character(node_id)]]
        node$path <- tmp[n, 'tree']
        node$node_id <- node_id
        return(node)
    })

    output$distPlot <- renderPlot({
        message("In renderPlot")
        node <- nodeData()
        if (is.null(node))
             return()
        message(paste0(names(node), collapse=","))
        model <- c(c=input$model.c, m=input$model.m, r=input$model.r, b=input$model.b)

        message(sprintf("model = %s", paste0(model, collapse=", ")))

        if (is.null(model) || any(sapply(model, is.null)))
            return()

        if (is.null(input$fitting))
            return()
            
        bls <- node$bls
        fit = node$fit
        message(sprintf("fitting = %s", input$fitting))

        # plot the likelihood estimate using the model parameters as adjusted by the user.
        #
        switch(input$fitting,
               'unweighted'= {
                   lcfit.results <- calculate_lcfit(bls, model, fit, weighted=F)
                   bls$fit_ll <- lcfit.results$ll
               },
               'weighted' = {
                   lcfit.results <- calculate_lcfit(bls, model, fit, weighted=T)
                   bls$wfit_ll <- lcfit.results$ll},
               'top4' = {
                   lcfit.results <- calculate_lcfit(bls, model, fit, weighted=F, keep=4)
                   bls$t4fit_ll <- lcfit.results$ll})

        # If the model has been modified, plot the baseline plot
        # using the default model.   That way the user can always see where they are coming from.
        # 
        if (any(model != DEFAULT_MODEL)) {
            message(sprintf("model = %s", paste0(model, collapse=", ")))
                    
            switch(input$fitting,
                   'unweighted'= {
                       bls$baseline_ll <- calculate_lcfit(bls, DEFAULT_MODEL, fit, weighted=F)$ll
                   },
                   'weighted' = {
                       bls$baseline_ll <- calculate_lcfit(bls, DEFAULT_MODEL, fit, weighted=T)$ll
                   },
                   'top4' = {
                       bls$baseline_ll <- calculate_lcfit(bls, DEFAULT_MODEL, fit, weighted=F, keep=4)$ll
                   })
        } else if ('baseline_ll' %in% names(bls))
            bls <- bls[,-grep('baseline_ll', names(bls)),drop=FALSE]


        title <-  sprintf("%s #%s",  node$path, node$node_id)
        subtitle <- sprintf("c=%s m=%s r=%s b=%s", model['c'], model['m'], model['r'], model['b'] )
        print(sprintf("lcfit status = %s", lcfit.results$status))
        status <- switch(as.character(lcfit.results$status),
                         '0'="Branch length",
                         '1'="exceeded maximum iterations without converging",
                         '2'="non-specific error",
                         '27'="iteration is not making progress towards solution",
                         '29'="cannot reach the specified tolerance in F")
        p <- ggplot() + ylab("Log likelihood")
        p <- p + geom_line(aes(x=branch_length, y=value, color=variable, linetype=variable),
                           data=melt(bls, id.vars=1:2)) +
                ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
                theme(plot.title=element_text(size=15)) +
                geom_point(aes(x=branch_length, y=ll), data=fit) +
                scale_color_discrete(name="", breaks=c('bpp_ll', 'fit_ll', 'wfit_ll', 't4fit_ll', 'baseline_ll'), labels=c('bpp', 'unweighted', 'weighted', 'top4', 'baseline')) +
                guides(color="legend", shape=FALSE, linetype=FALSE) +
                theme(legend.position="bottom") +
                ylab("Log likelihood") +
                xlab(status) 

        print(p)
    })
})
