library(shiny)


shinyUI(fluidPage(theme = "bootstrap.css",
                  singleton(tags$head(HTML('
    <script type="text/javascript">
    $(document).ready(function() {
    
         // Disable an action button
         Shiny.addCustomMessageHandler("disableButton", function(message) {
         $("#"+message).attr("disabled", "true");
         });
    
         // Enable an action button
         Shiny.addCustomMessageHandler("enableButton", function(message) {
         $("#"+message).removeAttr("disabled");
         });

    })
    </script>')
                        )),
    titlePanel("LCfit Failure Mode Analysis"),
    fluidRow(
        column(4,
               wellPanel(
                   checkboxGroupInput("failures", "Failure modes:",
                                      c("Progress" = "ENOPROG",
                                        "Tolerance" = "ETOLF",
                                        "Iterations" = "MAXITER",
                                        "Other" = "ERROR"), selected=c("ENOPROG", "ETOLF",  "ERROR"), inline=TRUE),
                   htmlOutput("caption"),
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
                   uiOutput('max.iter'),

                   hr(),
                   actionButton("reset_input", "Reset")
                   )

               ),
        
        column(8,
               tabsetPanel(
                   tabPanel("Failure Plot",
                            fluidRow(
                                column(2, actionButton("left", "Prev")),
                                column(8,  htmlOutput("title")),
                                column(2, actionButton("right", "Next"))
                                ),
                            plotOutput("distPlot")),
                   tabPanel("Failure Map", 
                            actionButton("apply", "Apply current model to the failing nodes"),
                            plotOutput("failMap")),
                   tabPanel("Non-failures",
                            fluidRow(
                                span(style='float:left; ', actionButton("nf.apply", "Apply model")),
                                actionButton("nf.resample", "Resample"),
                                span(style='float:right; ',
                                     tags$label(style='vertical-align:top; display: inline-block;', class = "control-label", `for` = "nf.count", "# Cases: "),
                                     span(style='vertical-align:middle; display: inline-block;',
                                          selectInput("nf.count", label=NULL, choices=c(10, 50, 100, 300, 500, 1000))))
                                ),
                            plotOutput("successMap"))
                   )
               )
        )
    ))

