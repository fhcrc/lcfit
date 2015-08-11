library(shiny)
library(ggplot2)
library(grid)		# for unit() function
library(plyr)
library(reshape2)
library(RJSONIO)
source("../utils.R")

theme_set(theme_bw(16))

source("../aggutils.R")


m = NULL

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output, session) {

    # Implement the 'reset' button.
    # http://stackoverflow.com/a/24269691/1135316
    output$model_input <- renderUI({
        message("render model ui")
        times <- input$reset_input
        div(id=paste0("model.",letters[(times %% length(letters)) + 1]),
            sliderInput("model.c", "c: number of constant sites", 800, 2500, DEFAULT_MODEL['c'], step=2),
            sliderInput("model.m", "m: number of mutated sites", 650, 1400, DEFAULT_MODEL['m']),
            sliderInput("model.r", "r: mutation rate", 0, 5, DEFAULT_MODEL['r'], step=0.1),
            sliderInput("model.b", "b: a branch length offset", 0, 2, DEFAULT_MODEL['b'], step=0.1)
            )
    })

    output$max.iter <- renderUI({
        message("render max.iter ui")
        times <- input$reset_input
        div(id=paste0("max.iter.",letters[(times %% length(letters)) + 1]),
            sliderInput("max.iter", "Maximum iterations:", 1,
                        1000, 250)
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

        m <- join(m.measures, m.status, by=c("tree", "node_id", "fitting", "nextra"), type="left")

#        m <- m[!is.na(m$value),]
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



    # return a boolean vector that is TRUE for rows in loadAggData()
    # meeting the user-selected error criteria
    bSelectedNodes <- reactive({
        m <- loadAggData()
        
        ## let the user select which failure mode they are interested in.
        ## LCFIT_MAXITER = 1,	// exceeded maximum iterations without converging
        ## LCFIT_ERROR = 2,	// non-specific error
        ## LCFIT_ENOPROG = 27,	// iteration is not making progress towards solution
        ## LCFIT_ETOLF = 29	// cannot reach the specified tolerance in F
        errs = c(MAXITER=1, ERROR=2, ENOPROG=27, ETOLF=29)
        errs = errs[input$failures]
        # return a boolean vector (length=nrow(m)) that is TRUE for rows meeting the user-selected criteria
        bVector <- m$measure=="wrss" & m$nextra == input$nextra & m$fitting %in% input$fitting 
        if (any(is.na(errs))) {
            bVector <- bVector & (m$status %in% errs | is.na(m$value))
        } else {
            bVector <- bVector & m$status %in% errs
        }
        stopifnot(length(bVector) == nrow(m))
        message(sprintf("%d matching nodes", sum(bVector)))
        return(bVector)
    })


    # Sample nodes from the pool of non-failing nodes under the original model parameters.
    # This means we grab a subset of the nodes that did NOT fail when they were run with the original
    # defaul model parameters.   We will recalculate these nodes with the current modified model parameters to see if they still succeed.
    # We don't want our new model to make things worse!
    #
    # "input$nf.resample" resamples the nodes from the pool of
    # non-failing nodes.
    #
    # "input$nf.count" controls how many nodes to sample.  Constrained
    # by how long you are willing to wait for an interactive update
    # when applying the current model to these nodes.
    #
    nf.resample.value <- 0
    unselectedNodes <- reactive({
        
        if (input$nf.resample != nf.resample.value)
            nf.resample.value <<- input$nf.resample

        m <- loadAggData()
        # Calculate boolean index 
        b <- m$measure=="wrss" & m$nextra == input$nextra & m$fitting %in% input$fitting & m$status==0
        stopifnot(length(b) == nrow(m))
        i <- which( b )   		# convert to indicies
        i <- sample(i, min(input$nf.count, length(i)))	# sub-sample the indicies
        return(m[i,])
    })



    # Return a dataframe of nodes that are failing according to the user-selected criteria.
    #
    # We will select one of these nodes to display a detailed graph.
    # When directed to do so, we will apply a modified model to all
    # the nodes in this list to see how many failures are corrected.
    #
    selectedNodes <- reactive({
        i <- bSelectedNodes()
        m <- loadAggData()
        tmp <- m[i,]
        if (nrow(tmp) == 0) 
            return()
        tmp <- tmp[order(tmp$value, decreasing=F),]

        # initialize the index to a random node from among those selected 
        inode <<- sample(1:nrow(tmp), 1)
        
        return(tmp)

    })

    # save the value of the action buttons so we know when they are clicked.
    iprev.save = 0
    inext.save = 0
    
    # React to the user clicking on the 'previous' or 'next' buttons
    # This updates the global variable `inode` to keep track of
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

        # disable or enable the 'next' and 'previous' buttons as
        # appropriate depending on whether we are currently at the
        # first or last node in the list.
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
    
    # Load the data for a single node from selectedNodes()
    #
    # We are just interested in looking at all sorts of failing nodes,
    # so the particular tree has no significance beyond holding a node
    # that failed lcfit.
    #
    nodeData <- reactive({
        tmp <- selectedNodes()
        if (is.null(tmp))
            return()
        n <- iNode()
        message(sprintf("selecting row %d out of %d", n, nrow(tmp)))

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
        node$fitting <- tmp[n,'fitting']
        return(node)
    })


    # Render a title consisting of the name of the file and the current model parameters.
    # This title changes when a different plot is selected, or when the model parameters change.
    output$title <- renderUI({
        node <- nodeData()
        if (is.null(node))
             return()
        model <- c(c=input$model.c, m=input$model.m, r=input$model.r, b=input$model.b)
        iterations = input$max.iter
        if (is.null(model) || any(sapply(model, is.null)))
            return()

        list(tags$div(style="text-align: center; font-style: bold; font-size: 120%;",
                      sprintf("%s #%s",  node$path, node$node_id)),
             tags$div(style="text-align: center; font-style: italic;",
                      sprintf("c=%s m=%s r=%s b=%s  iterations=%d", model['c'], model['m'], model['r'], model['b'], iterations )))
    })



    ##' convenience function for fitting a model to data points according to style.
    ##'
    ##' @param type one of 'weighted', 'unweighted', or 'top4'
    ##' @param pts sampled points; dataframe { 'x', 'y' }
    ##' @param model named numeric vector, { 'c', 'm', 'r', 'b' }
    ##' @return fitted model numeric vector, { 'c', 'm', 'r', 'b', 'status' }
    ##' @author chris
    fit.model.bytype <- function(type, pts, model, max.iter) {
        names(pts) <- c('x', 'y')
        model <- switch(type,
                        weighted = fit_model(model, pts, weighted=T, max.iter),
                        unweighted = fit_model(model, pts, weighted=F, max.iter),
                        top4 = {
                                # Erick asks to keep only the four-highest likelihood points.
                            pts <- pts[order(pts$y, decreasing=T)[1:4],]
                            fit_model(model, pts, weighted=F, max.iter)
                        }
                        )
        return(model)
    }


    
    
    output$distPlot <- renderPlot({
        message("In output$distPlot")
        node <- nodeData()
        if (is.null(node))
             return()
        model.input <- c(c=input$model.c, m=input$model.m, r=input$model.r, b=input$model.b)
        max.iter <-  input$max.iter
        message(sprintf("starting model = %s", paste0(model.input, collapse=", ")))

        if (is.null( model.input) || any(sapply( model.input, is.null)))
            return()

        bls <- node$bls
        pts <- node$fit[,c('branch_length', 'll')]
        message(sprintf("fitting = %s", node$fitting))

        # If the model has been modified, plot the baseline plot
        # using the default model.   That way the user can always see where they are coming from.
        # 
        if (any(model.input != DEFAULT_MODEL) && !input$suppress.baseline) {
            model <- fit.model.bytype(node$fitting, pts, DEFAULT_MODEL, 250)
            bls$baseline_ll <- sapply(node$bls[['branch_length']], lcfit, model=model)
            
        } else if ('baseline_ll' %in% names(bls)) {
            # otherwise remove the baseline so it doesn't show up in the graph
            bls <- bls[,-grep('baseline_ll', names(bls)),drop=FALSE]
        }
        # plot the likelihood estimate using the model parameters as adjusted by the user.
        #
        model <- fit.model.bytype(node$fitting, pts, model.input, max.iter)
        message(sprintf("fitted model = %s", paste0(model, collapse=", ")))
        bls$fit_ll<- sapply(node$bls[['branch_length']], lcfit, model=model)
        
        print(sprintf("lcfit status = %s", model[['status']]))
        status <- switch(as.character(model[['status']]),
                         '0'="Branch length",
                         '1'="exceeded maximum iterations without converging",
                         '2'="non-specific error",
                         '27'="iteration is not making progress towards solution",
                         '29'="cannot reach the specified tolerance in F")
        p <- ggplot() + ylab("Log likelihood")
        p <- p + geom_line(aes(x=branch_length, y=value, color=variable, linetype=variable),
                           data=melt(bls, id.vars=1:2)) +
                geom_point(aes(x=branch_length, y=ll), data=pts) +
                scale_color_discrete(name="", breaks=c('bpp_ll', 'baseline_ll', 'fit_ll'), labels=c('bpp', 'baseline', 'fit')) +
                guides(color="legend", shape=FALSE, linetype=FALSE) +
                theme(legend.position="bottom") +
                ylab("Log likelihood") +
                xlab(status) 

        print(p)
    })


    recalculateNodes <- function(nodes, model, fitting, max.iter) {

        # helper function to apply to each row in `nodes`
        recalculate.helper <- function(row, model, fitting, max.iter) {
             # browser()
            tree <- readSimulationData(row['tree'])
            message(paste0("reading tree ", row['tree']))
            node_id <- as.character(as.numeric(row['node_id']))
            message(paste0("node_id ", node_id))
            node <- tree[[as.character(node_id)]]
            bls <- node$bls
            pts <- node$fit[,c('branch_length', 'll')]
            model <- fit.model.bytype(fitting, pts, model, max.iter)
            status = model[['status']]

            return(status)
        }
        message(paste0(names(nodes), sep=", "))

        wd = getwd()
        tryCatch( {
            setwd("../../")
            v = apply(nodes, 1, recalculate.helper, model, fitting, max.iter)
        }, finally={ message("TryCatch!!!")
                     setwd(wd) })

        message(class(v))
        return(v)
    }


    nf.apply.value <- 0
    failures.action.value <- 0

    # Get the success/failure status of a random sampling of the nodes NOT under the current selection criteria.
    #
    # Unless the actionbutton has been pressed, this status reflects
    # the sttaus calculated under the standard model.  After pressing
    # the action button, each node is recalculated with the current
    # adjusted model.  Since this can be an expensive operation, we
    # only do so when the tiobutton is pressed, and not when the model
    # is changed.
    #
    successNodes <- reactive({
        tmp <- unselectedNodes()
        if (is.null(tmp))
            return()

        # handle the 'Apply.." action button on the "Non-failures" panel
        if (input$nf.apply != nf.apply.value) {
            nf.apply.value <<- input$nf.apply
            print("applying model!")
            model <- isolate(c(c=input$model.c, m=input$model.m, r=input$model.r, b=input$model.b))
            max.iter <- isolate(input$max.iter)
            fitting <- isolate(input$fitting)
            
            message(sprintf("model = %s", paste0(model, collapse=", ")))

            if (is.null(model) || any(sapply(model, is.null)))
                return()

            if (is.null(fitting))
                return()

            return(recalculateNodes(tmp, model, fitting, max.iter))
        }
        
        return(tmp$status)
    })



    # Get the success/failure status of all the nodes under the current selection criteria.
    #
    # Unless the actionbutton has been pressed, this status reflects
    # the sttaus calculated under the standard model.  After pressing
    # the action button, each node is recalculated with the current
    # adjusted model.  Since this can be an expensive operation, we
    # only do so when the tiobutton is pressed, and not when the model
    # is changed.
    #
    failingNodes <- reactive({
        tmp <- selectedNodes()
        if (is.null(tmp))
            return()
        
        # handle the 'apply' action button on the "Failure Map" panel.
        if (input$apply != failures.action.value) {
            failures.action.value <<- input$apply
            print("applying model!")
            model <- isolate(c(c=input$model.c, m=input$model.m, r=input$model.r, b=input$model.b))
            fitting <- isolate(input$fitting)
            max.iter <- isolate(input$max.iter)
            
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
                pts <- node$fit[,c('branch_length', 'll')]
                fitted <- fit.model.bytype(fitting, pts, model, max.iter)
                tmp[i,"status"] <- fitted[['status']]

            }
        }

        return(tmp$status)
    })

    tilePlot <- function(status.vector, show.rect) {
        n <- ceiling(sqrt(length(status.vector)))
        df <- expand.grid(x=seq(n), y=seq(n))
        df$status <- factor(c(status.vector, rep(0, n**2 - length(status.vector))))
        p <- ggplot(data=df[1:length(status.vector),], aes(x=x, y=y))
        p <- p + geom_tile(aes(fill=status))

        # outline the tile corresponding to the currently selected tile
        if (!missing(show.rect)) 
            p <- p + geom_rect(data=df[show.rect,], size=1, fill=NA, colour="black",
                               aes(xmin=x - 0.5, xmax=x + 0.5, ymin=y - 0.5, ymax=y + 0.5)) 

        p <- p +  theme(panel.grid = element_blank())
        p <- p +  theme(strip.background = element_blank())
        p <- p + scale_fill_manual(values=c("0"="#1A9641", "1"="#FFB54F", "2"="#FFC400", "27"="#FFD300", "29"="#FFE000"),
                                   breaks=c(0, 1, 2, 27, 29),
                                   labels=c("Success!", "Iterations", "Other", "Progress", "Tolerance"))
        p <- p + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
        p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
        p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
        p <- p + theme(panel.margin = unit(0.1, "lines"))
        p <- p + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"))
        return(p)
    }
    
    # http://stackoverflow.com/a/13261443/1135316
    output$failMap <- renderPlot({
        print("in output$failMap")
        tmp <- failingNodes()
        if (is.null(tmp))
            return()
        p <- tilePlot(tmp, show.rect=iNode())
        p <- p + ggtitle("Failures in the current selection")

        print(p)
    })

    # http://stackoverflow.com/a/13261443/1135316
    output$successMap <- renderPlot({
        print("in output$successMap")
        tmp <- successNodes()
        if (is.null(tmp))
            return()
        p <- tilePlot(tmp)
        p <- p + ggtitle("Previously successful nodes")
        print(p)
    })
    

})
