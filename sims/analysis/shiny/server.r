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
shinyServer(function(input, output, session) {

    # Implement the 'reset' button.
    # http://stackoverflow.com/a/24269691/1135316
    output$model_input <- renderUI({
        message("render model ui")
        times <- input$reset_input
        div(id=letters[(times %% length(letters)) + 1],
            sliderInput("model.c", "c: number of constant sites", 800, 2500, DEFAULT_MODEL['c'], step=2),
            sliderInput("model.m", "m: number of mutated sites", 650, 1400, DEFAULT_MODEL['m']),
            sliderInput("model.r", "r: mutation rate", 0, 5, DEFAULT_MODEL['r'], step=0.1),
            sliderInput("model.b", "b: a branch length offset", 0, 2, DEFAULT_MODEL['b'], step=0.1)
            )
    })

    ## load data for all the nodes in all the trees.
    loadAggData <- reactive({
        aggfile <- "../../aggfit.csv"
        aggdata <- read.csv(aggfile, stringsAsFactors=F)
        message(sprintf("loaded %d rows from %s", nrow(aggdata), aggfile))
        aggdaa <- aggdata[1:2000,]
        m.measures <- melt(aggdata, id.vars=c("tree", "node_id"), measure.vars=names(aggdata)[!grepl('tree|node_id|status_', names(aggdata))])
        m.measures <- cbind(m.measures, colsplit(m.measures$variable, names = c("measure", "fitting", "nextra"), pattern="_"))
        m.measures = m.measures[,-grep('variable',names(m.measures))]
        
        m.status <- melt(aggdata, id.vars=c("tree", "node_id"), measure.vars=names(aggdata)[grepl('status_', names(aggdata))], value.name='status')
        m.status <- cbind(m.status, colsplit(m.status$variable, names = c("measure", "fitting", "nextra"), pattern="_"))
        m.status = m.status[,-grep('variable|measure',names(m.status))]

        m <- join(m.measures, m.status, by=c("tree", "node_id", "fitting", "nextra"), type="left")

        m <- m[!is.na(m$value),]
        return(m)
    })


    output$caption <- renderUI({
        m <- loadAggData()
        ntotal <- nrow(m)
        nfailures <- sum(m$status != 0)
        nselected <- nrow(selectedNodes())

        tags$p(sprintf("%d selected out of %d failures, %d total samples.", nselected, nfailures, ntotal))
    })

    inode = 0

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

        # initialize the index to a random node from among those selected 
        inode <<- sample(1:nrow(tmp), 1)
        
        return(tmp)

    })

    iprev.save = 0
    inext.save = 0
    
    # React to the user clicking on the 'previous' or 'next' buttons
    # This updates the global variables `inode` keep track of
    # which node is currently selected.
    iNode <- reactive({
        tmp <- selectedNodes()
        if (is.null(tmp))
            return()

        iprev <- input$left
        inext <- input$right
        if (iprev != iprev.save) {
            iprev.save <<- iprev
            inode <<- inode - 1
            if (inode < 1) 
                inode <<- 1
        }
        if (inext != inext.save) {
            inext.save <<- inext
            inode <<- inode + 1
            if (inode > nrow(tmp))
                inode <<- nrow(tmp)
        }
        if (inode == 1) {
            session$sendCustomMessage("disableButton", "left")
        } else {
            session$sendCustomMessage("enableButton", "left")
        }
        if (inode == nrow(tmp)) {
            session$sendCustomMessage("disableButton", "right")
        } else {
            session$sendCustomMessage("enableButton", "right")
        }
        return(inode)

    })
    
    # Load the data for a single node in a tree.
    # We are just interested in looking at all sorts of failing nodes, so the particular tree has
    # significance beyond holding a node that failed lcfit.
    #
    nodeData <- reactive({
        tmp <- selectedNodes()
        if (is.null(tmp))
            return()
        n <- iNode()
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


    # Render a title consisting of the name of the file and the current model parameters.
    # This title changes when a different plot is selected, or when the model parameters change.
    output$title <- renderUI({
        node <- nodeData()
        if (is.null(node))
             return()
        model <- c(c=input$model.c, m=input$model.m, r=input$model.r, b=input$model.b)
        if (is.null(model) || any(sapply(model, is.null)))
            return()

        list(tags$div(style="text-align: center; font-style: bold; font-size: 120%;",
                      sprintf("%s #%s",  node$path, node$node_id)),
             tags$div(style="text-align: center; font-style: italic;",
                      sprintf("c=%s m=%s r=%s b=%s", model['c'], model['m'], model['r'], model['b'] )))
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
                geom_point(aes(x=branch_length, y=ll), data=fit) +
                scale_color_discrete(name="", breaks=c('bpp_ll', 'fit_ll', 'wfit_ll', 't4fit_ll', 'baseline_ll'), labels=c('bpp', 'unweighted', 'weighted', 'top4', 'baseline')) +
                guides(color="legend", shape=FALSE, linetype=FALSE) +
                theme(legend.position="bottom") +
                ylab("Log likelihood") +
                xlab(status) 

        print(p)
    })

    failures.action.value <- 0
    
    failingNodes <- reactive({
        tmp <- selectedNodes()
        if (is.null(tmp))
            return()
        
        if (input$apply != failures.action.value) {
            failures.action.value <<- input$apply
            print("applying model!")
            model <- isolate(c(c=input$model.c, m=input$model.m, r=input$model.r, b=input$model.b))
            fitting <- isolate(input$fitting)
            
            message(sprintf("model = %s", paste0(model, collapse=", ")))

            if (is.null(model) || any(sapply(model, is.null)))
                return()

            if (is.null(fitting))
                return()

            for (i in 1:nrow(tmp)) {
                # read in the sampled points for that node in that particular tree...
                wd = getwd()
                tryCatch( {
                    setwd("../../")
                    tree <- readSimulationData(tmp[i,'tree'])
                }, finally={ setwd(wd) })
                node_id <- as.character(tmp[i, 'node_id'])
                node <- tree[[as.character(node_id)]]

                bls <- node$bls
                fit = node$fit

                tmp[i,"status"] <- switch(fitting,
                                          'unweighted'= {
                                              calculate_lcfit(bls, model, fit, weighted=F)$status
                                          },
                                          'weighted' = {
                                              calculate_lcfit(bls, model, fit, weighted=T)$status
                                          },
                                          'top4' = {
                                              calculate_lcfit(bls, model, fit, weighted=F, keep=4)$status
                                          })
            }
        }

        return(tmp$status)

    })

    # http://stackoverflow.com/a/13261443/1135316
    output$failMap <- renderPlot({
        print("in renderPlot")
        tmp <- failingNodes()
        if (is.null(tmp))
            return()
        n <- ceiling(sqrt(length(tmp)))
        df <- expand.grid(x=seq(n), y=seq(n))
        df$status <- factor(c(tmp, rep(0, n**2 - length(tmp))))
        p <- ggplot(data=df[1:length(tmp),], aes(x=x, y=y))
        p <- p + geom_tile(aes(fill=status))
        p <- p + geom_rect(data=df[iNode(),], size=1, fill=NA, colour="black",
                           aes(xmin=x - 0.5, xmax=x + 0.5, ymin=y - 0.5, ymax=y + 0.5)) 
        p <- p +  theme(panel.grid = element_blank())
        p <- p +  theme(strip.background = element_blank())
        p <- p + scale_fill_hue(h=c(20, 60), l=80, c=150, breaks=c(1, 2, 27, 29, 0),
            labels=c("Iterations", "Other", "Progress", "Tolerance", "Success!"))
        p <- p + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
        p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
        p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
        p <- p + ggtitle("Failures in the current selection")
        p <- p + theme(panel.margin = unit(0.1, "lines"))
        p <- p + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"))
        print(p)
    })

})
