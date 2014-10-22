library(shiny)

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(

    # Application title
    headerPanel("LCfit Failure Mode Analysis"),

    # Sidebar with a slider input for number of observations
    sidebarPanel(
        wellPanel(
            checkboxGroupInput("failures", "Failure modes:",
                               c("Progress" = "ENOPROG",
                                 "Tolerance" = "ETOLF",
                                 "Iterations" = "MAXITER",
                                 "Other" = "ERROR"), selected=c("ENOPROG", "ETOLF",  "ERROR"), inline=TRUE),
            code(textOutput("caption")),
            radioButtons("fitting", "Fitting technique:",
                         c("Unweighted" = "unweighted",
                           "Weighted" = "weighted",
                           "Top 4" = "top4"), selected=c("unweighted"), inline = TRUE),
            sliderInput("nextra", "Extra points:", 0, 4, 0)
            ),

        hr(),
        wellPanel(
            helpText(paste0(
                "Select initial model parameters.")),
            uiOutput('model_input'),
            hr(),
            sliderInput("niter", "Maximum iterations:", 100,
                        200, 120),
            hr(),
            actionButton("reset_input", "Reset")
            )
        ),

    # Show a plot of the generated distribution
    mainPanel(
        plotOutput("distPlot")
        )
    ))



