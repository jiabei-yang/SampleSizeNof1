library(shiny)
library(dplyr)
library(ggplot2)
library(forcats)     # for fct_explicit_na()
library(plotly)
library(stringr)     # for str_match() used in ggplotly legend
# library(shinycssloaders)     # withSpinner()
# library(shinyjs)             # hide()

source("gen.R")
source("Power.R")
source("se.R")

###################################################################################################
########################################## Define UI ##############################################
###################################################################################################
ui <- fluidPage(
  titlePanel("Sample Size Calculations for N-of-1 Trials"),
  
  # useShinyjs(),
  # extendShinyjs(text = "shinyjs.resetClick = function() { Shiny.onInputChange('plotly_click-A', 'null'); }",
  #               functions = "resetClick"),
  
  sidebarLayout(
    sidebarPanel(
      
      helpText("This app designs a series of n-of-1 trials using the following steps:", br(),
               "1) Find the optimized designs for estimating the population average treatment effect", br(),
               "2) Given # of measurements per participant or # of participants in optimized designs, finalize the design which satisfies the standard error requirement for individual-specific treatment effects."),
      
      # radioButtons("trt.eff.intst",
      #              label   = "Treatment effect of interest",
      #              choices = list("Population average treatment effect" = "pop.avg",
      #                             "Individual-specific treatment effect" = "indiv.spec")),
      
      # Appear if population average treatment effect is selected
      # conditionalPanel(condition = "input['trt.eff.intst'] == 'pop.avg'",
      
      # Possible sequences
      radioButtons("psbl.seq",
                   label   = "Possible sequences",
                   choices = list("Pairwise randomization" = "Pairwise Randomization",
                                  "Alternating sequences"  = "Alternating Sequences",
                                  "User-specified sequences" = "User-Specified Sequences")),
      
      conditionalPanel(condition = "input['psbl.seq'] == 'User-Specified Sequences'",
                       
                       helpText("Uploaded file should", br(),
                                "1) be in .csv format", br(), 
                                "2) have header (no restrictions on header names)", br(),
                                "3) have columns correspond to periods in sequences", br(),
                                "4) have rows correspond to sequences", br(),
                                "5) use 0 and 1 to indicate reference and intervention treatments respectively"), 
                       
                       fileInput("file.psbl.seq",
                                 label = "Upload the sequences in the design:")
                       
                       
      ),
      
      # fluidRow(
      #   column(6,
      #          radioButtons("model.pop.avg.intcpt",
      #                       label   = "Form of intercept",
      #                       choices = list("Fixed"  = "fixed",
      #                                      "Random" = "random"))),
      #   
      #   column(6,
      #          radioButtons("model.pop.avg.slp",
      #                       label   = "Form of slope",
      #                       choices = list("Common"  = "common",
      #                                      "Random"  = "random")))
      # ),
      
      checkboxGroupInput("model.popavg", 
                         label    = "Model form", 
                         choices  = list("Fixed intercepts-Common slope"  = "fixed-common", 
                                         "Random intercepts-Common slope"  = "random-common", 
                                         "Fixed intercepts-Random slopes"  = "fixed-random",
                                         "Random intercepts-Random slopes" = "random-random"),
                         selected = c("fixed-common", "random-common", "fixed-random", "random-random")),
      
      conditionalPanel(condition = "input['model.popavg'].includes('random-common') || input['model.popavg'].includes('random-random')",
                       div(style = "white-space: nowrap;", 
                           numericInput("var.intcpt", 
                                        label = "Variance of random intercepts",
                                        value = 4,
                                        width = 100)
                       )
      ),
      
      conditionalPanel(condition = "input['model.popavg'].includes('fixed-random') || input['model.popavg'].includes('random-random')",
                       div(style = "white-space: nowrap;", 
                           numericInput("var.slp", 
                                        label = "Variance of random slopes",
                                        value = 1,
                                        width = 100)
                       )
      ),
      
      conditionalPanel(condition = "input['model.popavg'].includes('random-random')",
                       div(style = "white-space: nowrap;", 
                           numericInput("cov.intcpt.slp", 
                                        label = "Covariance of random intercepts and random slopes",
                                        value = 1,
                                        width = 100)
                       )
      ),
      
      fluidRow(
        column(6,
               numericInput("type.1.err", 
                            label = "Type I error",
                            value = 0.05,
                            width = 100, 
                            min   = 0,
                            max   = 1,
                            step  = 0.1)),
        
        column(6,
               numericInput("power",
                            label = "Power",
                            value = 0.8,
                            width = 100,
                            min   = 0,
                            max   = 1,
                            step  = 0.1))
      ),
      # ),
      
      # # Appear if individual-specific treatment effect is selected
      # conditionalPanel(condition = "input['trt.eff.intst'] == 'indiv.spec'",
      #                  radioButtons("indiv.spec.estr",
      #                               label   = "Individual-specific treatment effect estimator",
      #                               choices = list("Naive estimate" = "naive",
      #                                              "Shrunken estimate" = "shrunken")),
      #                  
      #                  conditionalPanel(condition = "input['indiv.spec.estr'] == 'shrunken'",
      #                                   fluidRow(
      #                                     column(4, 
      #                                            radioButtons("shrunk.indiv.spec.intcpt",
      #                                                         label   = "Form of intercept",
      #                                                         choices = list("Fixed"  = "fixed",
      #                                                                        "Random" = "random"))),
      #                                     
      #                                     column(8, 
      #                                            strong("Form of slope"),
      #                                            helpText("The slope can only be random to estimate individual-specific treatment effect."))
      #                                   )
      #                  )
      # ),
      
      radioButtons("homo.resid.err",
                   label   = "Homogeneity of residual errors",
                   choices = list("Homogeneous"   = "homo",
                                  "Heterogeneous (Temporarily not implemented)" = "hetero")),
      
      conditionalPanel(condition = "input['homo.resid.err'] == 'homo'",
                       
                       div(style = "white-space: nowrap;", 
                           numericInput("sigma.resid.err",
                                        label = "Homogeneous residual standard error",
                                        value = 2,
                                        width = 100,
                                        min   = 0)
                       )
                       
      ),
      
      radioButtons("CorrStr",
                   label    = "Correlation structure for repeated measurements",
                   choices  = list("AR-1"         = "AR-1",
                                   "Exchangeable" = "Exchangeable",
                                   "Independent" = "Independent")),
      
      conditionalPanel(condition = "(input.CorrStr.indexOf('Exchangeable') != -1) || (input.CorrStr.indexOf('AR-1') != -1)",
                       
                       numericInput("corr", 
                                    label = "Correlation",
                                    value = 0.40,
                                    width = 100,
                                    min   = 0,
                                    max   = 1,
                                    step  = 0.1)
      ),
      
      # conditionalPanel(condition = "input['trt.eff.intst'] == 'pop.avg'", 
      
      div(style = "white-space: nowrap;", 
          numericInput("min.diff",
                       label = "Minimal clinically important treatment effect",
                       value = 1, 
                       width = 100,
                       min   = 0)
      ),
      
      radioButtons("pop.avg.plot",
                   label   = "Design option",
                   choices = list("Total # of measurements vs. # of measurements per participant" = "ttobs.vs.ppobs",
                                  "Total # of measurements vs. # of participants"                 = "ttobs.vs.pts")),
      
      conditionalPanel(condition = "input['pop.avg.plot'] == 'ttobs.vs.ppobs'",
                       sliderInput("range.ppobs",
                                   label = "Range for number of measurements per participant",
                                   min   = 2,
                                   max   = 500,
                                   value = c(2, 100)),
                       
                       div(style = "white-space: nowrap;", 
                           numericInput("max.IJ",
                                        label = "Maximum number of participants",
                                        value = 200, 
                                        width = 100,
                                        min   = 0)
                       )
      ),
      
      conditionalPanel(condition = "input['pop.avg.plot'] == 'ttobs.vs.pts'",
                       sliderInput("range.pts",
                                   label = "Range for number of participants",
                                   min   = 2,
                                   max   = 300,
                                   value = c(2, 100)),
                       
                       div(style = "white-space: nowrap;", 
                           numericInput("max.KL",
                                        label = "Maximum number of measurements per participant",
                                        value = 100, 
                                        width = 100,
                                        min   = 0)
                       )
                       
      ),
      
      # ),
      
      strong("Only optimize designs over y-axis"),
      checkboxInput("show.all.design", 
                    "Yes", value = FALSE),
      
      strong("Calculate standard error for individual-specific treatment effect estimates"),
      checkboxInput("calc.indiv.prec", 
                    "Yes", value = FALSE),
      
      actionButton("action.calc",
                   label = "Calculate"),
      
      actionButton("clear",
                   label = "Clear")
      
    ), # sidebarPanel(
    
    mainPanel( 
      
      tableOutput("params"),
      
      plotlyOutput("fig", width = "100%", height = 600),
      
      tags$br(),
      
      uiOutput("plot.note"),
      
      tags$br(),
      
      uiOutput("click.caption"),
      tableOutput("psbl.scnr"),
      
      tags$br(),
      uiOutput("click.caption.2"),
      plotlyOutput("fig.se", width = "100%", height = 600),
      
      tags$br(),
      uiOutput("click.caption.se.ovrl"),
      tableOutput("df.selected.data.fig.se")
      
      # tags$br(),
      # plotlyOutput("fig.se.ind", width = "100%", height = 600),
      # 
      # tags$br(),
      # tableOutput("df.se")
      
    )
  )
)



###################################################################################################
###################################### Define server logic ########################################
###################################################################################################
server <- function(input, output) {
  
  # for clear button to work
  v <- reactiveValues(calc.ct        = 0)
  
  observeEvent(input$clear, {
    v$data.fig                 <- NULL
    v$summ.data.fig.opt.forFig <- NULL
    v$summ.data.fig.opt        <- NULL
    v$selected.models          <- NULL
    v$selected.scnr            <- NULL
    v$df.psbl.scnr             <- NULL
    v$data.fig.se              <- NULL
    v$summ.data.fig.se         <- NULL
    v$selected.summ.data.fig.se <- NULL
    v$selected.data.fig.se      <- NULL
    v$selected.data.fig.se.ind  <- NULL
  })
  
  # Only update the figure/output when the "Calculate" button is clicked
  observeEvent(input$action.calc, {
    
    # js$resetClick()
    # If first time click calculate button, show selected scnr
    # if not, show selected scnr only on clicking the plot
    # count # of times clicking calculate button and set no.click to T when clicking the button
    v$calc.ct  <- v$calc.ct + 1
    v$no.click <- T
    
    # # Population average treatment effect
    # if (input$trt.eff.intst == "pop.avg") {
    
    # input <- list()
    # input$CorrStr = "AR-1"
    # # input$psbl.seq = "Alternating Sequences"
    # input$psbl.seq = "Pairwise Randomization"
    # # input$psbl.seq <- "User-Specified Sequences"
    # input$file.psbl.seq$datapath <- "TestDesgin.csv"
    # input$range.ppobs <- c(2, 100)
    # input$range.pts <- c(2, 100)
    # k <- 2
    # input$corr <- 0.4
    # input$sigma.resid.err <- 2
    # input$type.1.err      <- 0.05
    # input$power           <- 0.8
    # input$min.diff        <- 1
    # # input$model.pop.avg.intcpt <- "fixed"
    # # input$model.pop.avg.slp    <- "common"
    # input$var.intcpt <- 4
    # input$var.slp        <- 1
    # input$cov.intcpt.slp <- 1
    # input$max.IJ <- 100
    # input$max.KL <- 200
    # input$model.popavg <- c("fixed-common", "random-common", "fixed-random", "random-random")
    # input$show.all.design <- T
    
    data.fig <- NULL
    models   <- data.frame(matrix(unlist(strsplit(input$model.popavg, "-")), byrow = T, nrow = length(input$model.popavg)))
    colnames(models) <- c("intcpt", "slp")
    # print(models)
    n.models <- nrow(models)
    
    # When the output is "total # of measurements vs. total # of measurements per pts"
    if (input$pop.avg.plot == "ttobs.vs.ppobs") {
      
      withProgress(message = 'Calculating', value = 0, {
        
        if (input$psbl.seq == "User-Specified Sequences") {
          
          file.psbl.seq <- read.csv(input$file.psbl.seq$datapath)
          k <- ncol(file.psbl.seq)
          I <- nrow(file.psbl.seq)
          
          # Find possible l values
          start.l <- ceiling(input$range.ppobs[1] / k)
          end.l   <- floor(input$range.ppobs[2] / k)
          
          # if start.l greater than end.l, only equal to end.l
          if (start.l <= end.l) {
            psbl.l <- start.l:end.l
          } else if ((start.l * k) <= input$range.ppobs[2]) {
            psbl.l <- start.l
            # if (end.l != 1) {
            #   print(c("end.l", k))
            # }
          } else if ((end.l * k) >= input$range.ppobs[1]) {
            psbl.l <- end.l
          } else {
            print(c("startEnd", k, input$range.ppobs))
            stop("Number of measurements per participant out of allowed range!")
          }
          lth.psbl.l <- length(psbl.l)
          
          # Models
          for (models.i in 1:n.models) {
            
            incProgress(1/n.models, detail = input$model.popavg[models.i])
            print(models[models.i, ]) 
            
            # assign NA to new max.J's for different models
            max.J        <- rep(NA, length(psbl.l))
            names(max.J) <- as.character(psbl.l)
            
            for (l in psbl.l) {
              
              max.J.ind <- which(names(max.J) == as.character(l))
              
              if (!is.na(max.J[max.J.ind])) {
                
                # only prev.max.J = 1 will come here, will not from l loop
                if (max.J[max.J.ind] == 1) {
                  
                  TmpDsgn <- c(J = 1,
                               K = k, 
                               L = l,
                               power = NA)
                  
                }  else {
                  
                  TmpDsgn <- calc_j_popavg(K               = k, 
                                           L               = l, 
                                           max.J           = max.J[max.J.ind],
                                           psbl.seq        = input$psbl.seq,
                                           intcpt          = models$intcpt[models.i],
                                           slp             = models$slp[models.i],
                                           corstr          = input$CorrStr, 
                                           corr.resid.err  = input$corr,
                                           sigma.resid.err = input$sigma.resid.err,
                                           var.intcpt      = input$var.intcpt,
                                           var.slp         = input$var.slp,
                                           cov.intcpt.slp  = input$cov.intcpt.slp,
                                           type.1.err      = input$type.1.err, 
                                           power           = input$power, 
                                           min.diff        = input$min.diff, 
                                           file.psbl.seq   = file.psbl.seq)
                  
                }
                
              } else { # is.na(max.J[max.J.ind])
                
                TmpDsgn <- calc_j_popavg(K               = k, 
                                         L               = l, 
                                         max.J           = NA,
                                         psbl.seq        = input$psbl.seq,
                                         intcpt          = models$intcpt[models.i],
                                         slp             = models$slp[models.i],
                                         corstr          = input$CorrStr, 
                                         corr.resid.err  = input$corr,
                                         sigma.resid.err = input$sigma.resid.err,
                                         var.intcpt      = input$var.intcpt,
                                         var.slp         = input$var.slp,
                                         cov.intcpt.slp  = input$cov.intcpt.slp,
                                         type.1.err      = input$type.1.err, 
                                         power           = input$power, 
                                         min.diff        = input$min.diff, 
                                         file.psbl.seq   = file.psbl.seq)
                
                # K               = k 
                # L               = l 
                # max.J           = NA
                # psbl.seq        = input$psbl.seq
                # intcpt          = models$intcpt[models.i]
                # slp             = models$slp[models.i]
                # corstr          = input$CorrStr 
                # corr.resid.err  = input$corr
                # sigma.resid.err = input$sigma.resid.err
                # var.intcpt      = input$var.intcpt
                # var.slp         = input$var.slp
                # cov.intcpt.slp  = input$cov.intcpt.slp
                # type.1.err      = input$type.1.err 
                # power           = input$power 
                # min.diff        = input$min.diff
                # file.psbl.seq   = file.psbl.seq
                
              } # else is.na(max.J[max.J.ind])
              
              max.J[max.J.ind] <- TmpDsgn["J"]
              
              # If l the last one, then cannot assign to the next value
              if (l != psbl.l[lth.psbl.l]) {
                
                # need to see if the next value is NA
                if (is.na(max.J[max.J.ind + 1])) {
                  max.J[max.J.ind + 1] <- TmpDsgn["J"]
                } else {
                  
                  # if the next value is not NA, only assign if the value is greater than the current value
                  if (max.J[max.J.ind + 1] > TmpDsgn["J"]) {
                    max.J[max.J.ind + 1] <- TmpDsgn["J"]
                  }
                }
              }
              
              data.fig <- rbind(data.fig, 
                                c(intcpt = as.character(models$intcpt[models.i]),
                                  slp    = as.character(models$slp[models.i]), 
                                  I      = I,
                                  TmpDsgn))
              
              # if I == 1, no need to run more since for larger K and L, I will be 1 as well
              if (TmpDsgn["J"] == 1) {
                
                # can break directly since we will delete the later observations anyway
                max.J[max.J.ind:lth.psbl.l] <- 1
                break
                
              }
              
            } # for l loop
            
            l <- l + 1
            # fill in all l possibilities for later se calculation
            while (l <= psbl.l[lth.psbl.l]) {
              # j must be 1 because break from the l loop
              data.fig <- rbind(data.fig,
                                c(intcpt = as.character(models$intcpt[models.i]),
                                  slp    = as.character(models$slp[models.i]), 
                                  I      = I,
                                  J      = 1,
                                  K      = k,
                                  L      = l,
                                  power  = NA))
              l <- l + 1
            }
            
          } # for models.i loop
          
        } else { # if input$psbl.seq == "other"
          
          for (models.i in 1:n.models) {
            
            incProgress(1/n.models, detail = input$model.popavg[models.i])
            
            print(models[models.i, ])
            # Record the previous min.J for last k to save computation time
            prev.max.J <- NULL
            
            # Big k/i loop, go through l to calculate j
            # At least 2 periods per participant
            for (k in 2:(input$range.ppobs[2])) {
              # for (k in 2:6) {
              
              I <- calc_I_popavg(k, input$psbl.seq)
              if (I > input$max.IJ) {
                break
              }
              
              if ((k %% 10) == 0) {
                print(k) 
              }
              
              # Find earlier I values, new I should be equal or smaller since we increase k and l
              start.l <- ceiling(input$range.ppobs[1] / k)
              end.l   <- floor(input$range.ppobs[2] / k)
              
              # if start.l greater than end.l, only equal to end.l
              if (start.l <= end.l) {
                psbl.l <- start.l:end.l
              } else if ((start.l * k) <= input$range.ppobs[2]) {
                psbl.l <- start.l
                # if (end.l != 1) {
                #   print(c("end.l", k))
                # }
              } else if ((end.l * k) >= input$range.ppobs[1]) {
                psbl.l <- end.l
              } else {
                print(c("startEnd", k, input$range.ppobs))
                stop("Number of measurements per participant out of allowed range!")
              }
              max.J        <- rep(NA, length(psbl.l))
              names(max.J) <- as.character(psbl.l)
              
              if (!is.null(prev.max.J)) {
                max.J[names(max.J) %in% names(prev.max.J)] <- prev.max.J[names(prev.max.J) %in% names(max.J)]
              }
              
              # IF psbl.l contains only 1 and max.J is 1, then can skip the following loops
              lth.psbl.l <- length(psbl.l)
              
              for (l in psbl.l) {
                # for (l in 16:20) {
                
                max.J.ind <- which(names(max.J) == as.character(l))
                
                # tmp = max.J
                # K               = k
                # L               = l
                # max.J           = max.J[max.J.ind]
                # psbl.seq        = input$psbl.seq
                # intcpt          = input$model.pop.avg.intcpt
                # slp             = input$model.pop.avg.slp
                # corstr          = input$CorrStr
                # corr.resid.err  = input$corr
                # sigma.resid.err = input$sigma.resid.err
                # var.intcpt      = input$var.intcpt
                # var.slp         = input$var.slp
                # cov.intcpt.slp  = input$cov.intcpt.slp
                # type.1.err      = input$type.1.err
                # power           = input$power
                # min.diff        = input$min.diff
                # max.J = tmp
                
                # When equal to na and not equal to 1, use function to calculate;
                # otherwise assign to 1 directly
                if (!is.na(max.J[max.J.ind])) {
                  
                  # only prev.max.J = 1 will come here, will not from l loop
                  if (max.J[max.J.ind] == 1) {
                    
                    TmpDsgn <- c(J = 1,
                                 K = k, 
                                 L = l,
                                 power = NA)
                    
                  } else {
                    
                    TmpDsgn   <- calc_j_popavg(K               = k, 
                                               L               = l, 
                                               max.J           = max.J[max.J.ind],
                                               psbl.seq        = input$psbl.seq,
                                               intcpt          = models$intcpt[models.i],
                                               slp             = models$slp[models.i],
                                               corstr          = input$CorrStr, 
                                               corr.resid.err  = input$corr,
                                               sigma.resid.err = input$sigma.resid.err,
                                               var.intcpt      = input$var.intcpt,
                                               var.slp         = input$var.slp,
                                               cov.intcpt.slp  = input$cov.intcpt.slp,
                                               type.1.err      = input$type.1.err, 
                                               power           = input$power, 
                                               min.diff        = input$min.diff)
                  }
                } else {
                  
                  TmpDsgn   <- calc_j_popavg(K               = k, 
                                             L               = l, 
                                             max.J           = max.J[max.J.ind],
                                             psbl.seq        = input$psbl.seq,
                                             intcpt          = models$intcpt[models.i],
                                             slp             = models$slp[models.i],
                                             corstr          = input$CorrStr, 
                                             corr.resid.err  = input$corr,
                                             sigma.resid.err = input$sigma.resid.err,
                                             var.intcpt      = input$var.intcpt,
                                             var.slp         = input$var.slp,
                                             cov.intcpt.slp  = input$cov.intcpt.slp,
                                             type.1.err      = input$type.1.err, 
                                             power           = input$power, 
                                             min.diff        = input$min.diff)
                }
                
                max.J[max.J.ind] <- TmpDsgn["J"]
                
                # If l the last one, then cannot assign to the next value
                if (l != psbl.l[lth.psbl.l]) {
                  
                  # need to see if the next value is NA
                  if (is.na(max.J[max.J.ind + 1])) {
                    max.J[max.J.ind + 1] <- TmpDsgn["J"]
                  } else {
                    
                    # if the next value is not NA, only assign if the value is greater than the current value
                    if (max.J[max.J.ind + 1] > TmpDsgn["J"]) {
                      max.J[max.J.ind + 1] <- TmpDsgn["J"]
                    }
                  }
                }
                
                data.fig <- rbind(data.fig, 
                                  c(intcpt = as.character(models$intcpt[models.i]),
                                    slp    = as.character(models$slp[models.i]), 
                                    I      = I,
                                    TmpDsgn))
                
                # if I == 1, no need to run more since for larger K and L, I will be 1 as well
                if (TmpDsgn["J"] == 1) {
                  
                  # can break directly since we will delete the later observations anyway
                  max.J[max.J.ind:lth.psbl.l] <- 1
                  break
                  
                }
                
              } # for l loop
              
              l <- l + 1
              # fill in all l possibilities for later se calculation
              while (l <= psbl.l[lth.psbl.l]) {
                # j must be 1 because break from the l loop
                data.fig <- rbind(data.fig,
                                  c(intcpt = as.character(models$intcpt[models.i]),
                                    slp    = as.character(models$slp[models.i]), 
                                    I      = I,
                                    J      = 1,
                                    K      = k,
                                    L      = l,
                                    power  = NA))
                l <- l + 1
              }
              
              prev.max.J <- max.J
              
              # if both L and J are 1, greater K will lead to L and J being 1 as well
              # therefore break k loop
              # Fill in all k/i possibilities, with l combinations for later se calculation
              # nrow.data.fig <- nrow(data.fig)
              if ((TmpDsgn["J"] == 1) & (TmpDsgn["L"] == 1)) {
                
                k <- k + 1
                I <- calc_I_popavg(k, input$psbl.seq)
                
                while ((I <= input$max.IJ) & (k <= input$range.ppobs[2])) {
                  
                  # Find earlier I values, new I should be equal or smaller since we increase k and l
                  start.l <- ceiling(input$range.ppobs[1] / k)
                  end.l   <- floor(input$range.ppobs[2] / k)
                  
                  # if start.l greater than end.l, only equal to end.l
                  if (start.l <= end.l) {
                    psbl.l <- start.l:end.l
                  } else {
                    psbl.l <- end.l
                    if (end.l != 1) {
                      print(c("end.l", k))
                    }
                  }
                  
                  for (l in psbl.l) {
                    
                    data.fig <- rbind(data.fig,
                                      c(intcpt = as.character(models$intcpt[models.i]),
                                        slp    = as.character(models$slp[models.i]), 
                                        I      = I,
                                        J      = 1,
                                        K      = k,
                                        L      = l,
                                        power  = NA))
                  }
                  
                  k <- k + 1
                  I <- calc_I_popavg(k, input$psbl.seq)
                } # while k loop
                
                break
              } # if ((TmpDsgn["J"] == 1) & (TmpDsgn["L"] == 1)) {
              
            } # for k loop
            
          } # for model.i loop
          
        } # else if input$psbl.seq != "other"
        
      }) # withProgress
      
      data.fig           <- data.frame(data.fig)
      # colnames(data.fig) <- c("K", "L", "J", "I", "Power", "CorrStr")
      # data.fig <- data.fig %>%
      #   mutate(I = as.integer(levels(I))[I],
      #          J = as.integer(levels(J))[J],
      #          K = as.integer(levels(K))[K],
      #          L = as.integer(levels(L))[L],
      #          power = as.numeric(levels(power))[power])
      data.fig <- data.fig %>%
        mutate(I = as.integer(I),
               J = as.integer(J),
               K = as.integer(K),
               L = as.integer(L),
               power = as.numeric(power))
      
      data.fig <- data.fig %>%
        filter(((I * J) <= input$max.IJ) &
                 ((K * L) >= input$range.ppobs[1]) &
                 ((K * L) <= input$range.ppobs[2]))
      data.fig <- data.fig %>%
        mutate(power = round(power, 2))
      data.fig   <- data.fig %>%
        mutate(model = ifelse(slp == "random", "Random slope model", "Common slope model"))
      data.fig <- data.fig %>%
        mutate(model = factor(model))
      
      # Optimize the designs
      # NEED TO MODIFY THIS FUNCTION FOR DIFFERENT COLUMNS IN DATA.FIG
      data.fig.opt <- optimize_dsgn(data.fig)
      # Summarize by desired x-axis
      v$summ.data.fig.opt <- gen_data_fig_ppobs(data.fig.opt)
      
      data.fig <- gen_data_fig_ppobs(data.fig)
      
      # When the output is "total # of measurements vs. # of pts"
    } else { # input$pop.avg.plot == "ttobs.vs.pts"
      
      withProgress(message = 'Calculating', value = 0, {
        
        if (input$psbl.seq == "User-Specified Sequences") {
          
          file.psbl.seq <- read.csv(input$file.psbl.seq$datapath)
          k <- ncol(file.psbl.seq)
          I <- nrow(file.psbl.seq)
          
          for (models.i in 1:n.models) {
            
            # Big j loop, go through k & i combination within j loop to calculate l
            # all possible k and j combinations will be in the generated data frame for standard error calculation
            incProgress(1/n.models, detail = input$model.popavg[models.i])
            
            j <- ceiling(input$range.pts[1] / I)
            
            TmpDsgn   <- calc_l_popavg(J               = j,
                                       K               = k, 
                                       max.L           = NA,
                                       max.KL          = input$max.KL, 
                                       psbl.seq        = input$psbl.seq,
                                       intcpt          = models$intcpt[models.i],
                                       slp             = models$slp[models.i],
                                       corstr          = input$CorrStr, 
                                       corr.resid.err  = input$corr,
                                       sigma.resid.err = input$sigma.resid.err,
                                       var.intcpt      = input$var.intcpt,
                                       var.slp         = input$var.slp,
                                       cov.intcpt.slp  = input$cov.intcpt.slp,
                                       type.1.err      = input$type.1.err, 
                                       power           = input$power, 
                                       min.diff        = input$min.diff,
                                       file.psbl.seq   = file.psbl.seq)
            
            # J               = j
            # K               = k
            # max.L           = NA
            # max.KL          = input$max.KL
            # psbl.seq        = input$psbl.seq
            # intcpt          = models$intcpt[models.i]
            # slp             = models$slp[models.i]
            # corstr          = input$CorrStr
            # corr.resid.err  = input$corr
            # sigma.resid.err = input$sigma.resid.err
            # var.intcpt      = input$var.intcpt
            # var.slp         = input$var.slp
            # cov.intcpt.slp  = input$cov.intcpt.slp
            # type.1.err      = input$type.1.err
            # power           = input$power
            # min.diff        = input$min.diff
            # file.psbl.seq   = file.psbl.seq
            
            prev.max.L <- TmpDsgn["L"]
            names(prev.max.L) <- k
            
            data.fig <- rbind(data.fig, 
                              c(intcpt = as.character(models$intcpt[models.i]),
                                slp    = as.character(models$slp[models.i]), 
                                I      = I,
                                TmpDsgn))
            
            # iteration j + 1, start from k = 2 again
            j <- j + 1
            
            # The maximum number of J should be range.pts/2
            # since the minimum number of different sequences is 2
            # I = 2 anyway when k = 2
            while (j <= (input$range.pts[2] / I)) {
              # while (j <= 18) {
              if ((j %% 5) == 0) {
                print(j) 
              }
              
              # fill in the missing J's and break out of the j loop
              if (prev.max.L == 1) {
                
                # fill in all the k possibilities for later standard error calculations
                while ((I * j) <= input$range.pts[2]) {
                  
                  data.fig <- rbind(data.fig,
                                    c(intcpt = as.character(models$intcpt[models.i]),
                                      slp    = as.character(models$slp[models.i]), 
                                      I      = I,
                                      J      = j,
                                      K      = k,
                                      L      = 1, 
                                      power  = NA)) 
                  
                  j <- j + 1
                } 
                
                break
                
              } else { # else if (prev.max.L != 1)
                
                TmpDsgn   <- calc_l_popavg(J               = j,
                                           K               = k, 
                                           max.L           = prev.max.L, 
                                           max.KL          = input$max.KL, 
                                           psbl.seq        = input$psbl.seq,
                                           intcpt          = models$intcpt[models.i],
                                           slp             = models$slp[models.i],
                                           corstr          = input$CorrStr, 
                                           corr.resid.err  = input$corr,
                                           sigma.resid.err = input$sigma.resid.err,
                                           var.intcpt      = input$var.intcpt,
                                           var.slp         = input$var.slp,
                                           cov.intcpt.slp  = input$cov.intcpt.slp,
                                           type.1.err      = input$type.1.err, 
                                           power           = input$power, 
                                           min.diff        = input$min.diff,
                                           file.psbl.seq   = file.psbl.seq)
                prev.max.L <- TmpDsgn["L"]
                
                data.fig <- rbind(data.fig, 
                                  c(intcpt = as.character(models$intcpt[models.i]),
                                    slp    = as.character(models$slp[models.i]), 
                                    I      = I,
                                    TmpDsgn))
              }
              
              j <- j + 1
              
            } # while j loop            
            
          } # for models.i loop
          
        } else { # else if (input$psbl.seq != "other")
          
          for (models.i in 1:n.models) {
            
            # Big j loop, go through k & i combination within j loop to calculate l
            # all possible k and j combinations will be in the generated data frame for standard error calculation
            incProgress(1/n.models, detail = input$model.popavg[models.i])
            
            print(models[models.i, ])
            k <- 2
            # I = 2 anyway
            I <- calc_I_popavg(k, input$psbl.seq)
            
            j <- ceiling(input$range.pts[1] / I)
            
            TmpDsgn   <- calc_l_popavg(J               = j,
                                       K               = k, 
                                       max.L           = NA,
                                       max.KL          = input$max.KL, 
                                       psbl.seq        = input$psbl.seq,
                                       intcpt          = models$intcpt[models.i],
                                       slp             = models$slp[models.i],
                                       corstr          = input$CorrStr, 
                                       corr.resid.err  = input$corr,
                                       sigma.resid.err = input$sigma.resid.err,
                                       var.intcpt      = input$var.intcpt,
                                       var.slp         = input$var.slp,
                                       cov.intcpt.slp  = input$cov.intcpt.slp,
                                       type.1.err      = input$type.1.err, 
                                       power           = input$power, 
                                       min.diff        = input$min.diff)
            
            prev.max.L <- TmpDsgn["L"]
            names(prev.max.L) <- k
            
            data.fig <- rbind(data.fig, 
                              c(intcpt = as.character(models$intcpt[models.i]),
                                slp    = as.character(models$slp[models.i]), 
                                I      = I,
                                TmpDsgn))
            
            last.L <- TmpDsgn["L"]
            k <- k + 1
            I <- calc_I_popavg(k, input$psbl.seq)
            
            # continue if 1) L > 1 or last.L is NA, 2) k <= max.KL NOT k * l <= max.KL because last.L might be NA, 3) i * j < range.pts(max)
            while (((last.L > 1) | is.na(last.L)) & 
                   (k <= input$max.KL) & 
                   ((I * j) <= input$range.pts[2])) {
              
              TmpDsgn <- calc_l_popavg(J               = j,
                                       K               = k, 
                                       max.L           = last.L,
                                       max.KL          = input$max.KL, 
                                       psbl.seq        = input$psbl.seq,
                                       intcpt          = models$intcpt[models.i],
                                       slp             = models$slp[models.i],
                                       corstr          = input$CorrStr, 
                                       corr.resid.err  = input$corr,
                                       sigma.resid.err = input$sigma.resid.err,
                                       var.intcpt      = input$var.intcpt,
                                       var.slp         = input$var.slp,
                                       cov.intcpt.slp  = input$cov.intcpt.slp,
                                       type.1.err      = input$type.1.err, 
                                       power           = input$power, 
                                       min.diff        = input$min.diff)
              
              prev.max.L <- c(prev.max.L, TmpDsgn["L"])
              names(prev.max.L)[length(prev.max.L)] <- k
              
              last.L     <- TmpDsgn["L"]
              data.fig <- rbind(data.fig,
                                c(intcpt = as.character(models$intcpt[models.i]),
                                  slp    = as.character(models$slp[models.i]), 
                                  I      = I,
                                  TmpDsgn))
              
              k <- k + 1
              I <- calc_I_popavg(k, input$psbl.seq)
            }
            
            # fill in all the k possibilities for later standard error calculations
            while ((k <= input$max.KL) & ((I * j) <= input$range.pts[2])) {
              
              data.fig <- rbind(data.fig,
                                c(intcpt = as.character(models$intcpt[models.i]),
                                  slp    = as.character(models$slp[models.i]), 
                                  I      = I,
                                  J      = j,
                                  K      = k,
                                  L      = last.L, 
                                  power  = NA)) 
              
              k <- k + 1
              I <- calc_I_popavg(k, input$psbl.seq)
            }
            
            # iteration j + 1, start from k = 2 again
            j <- j + 1
            
            # The maximum number of J should be range.pts/2
            # since the minimum number of different sequences is 2
            # I = 2 anyway when k = 2
            while (j <= (input$range.pts[2] / 2)) {
              # while (j <= 18) {
              if ((j %% 5) == 0) {
                print(j) 
              }
              
              max.L     <- prev.max.L
              lth.max.L <- length(max.L)
              
              # prepare while loop
              k <- 2 
              I <- calc_I_popavg(k, input$psbl.seq)
              max.L.ind <- which(names(max.L) == k)
              
              # stop iterating if 1) last.L drops to 1, 2) k goes out of range and 3) # of participants goes out of range
              while (((max.L[max.L.ind] > 1) | is.na(max.L[max.L.ind])) & 
                     (k <= min(input$max.KL, as.numeric(names(max.L)[lth.max.L]))) &
                     ((I * j) <= input$range.pts[2])) {
                
                # last.L from outside while loop
                max.L.ind <- which(names(max.L) == k)
                last.L    <- max.L[max.L.ind]
                
                TmpDsgn   <- calc_l_popavg(J               = j,
                                           K               = k, 
                                           max.L           = last.L, 
                                           max.KL          = input$max.KL, 
                                           psbl.seq        = input$psbl.seq,
                                           intcpt          = models$intcpt[models.i],
                                           slp             = models$slp[models.i],
                                           corstr          = input$CorrStr, 
                                           corr.resid.err  = input$corr,
                                           sigma.resid.err = input$sigma.resid.err,
                                           var.intcpt      = input$var.intcpt,
                                           var.slp         = input$var.slp,
                                           cov.intcpt.slp  = input$cov.intcpt.slp,
                                           type.1.err      = input$type.1.err, 
                                           power           = input$power, 
                                           min.diff        = input$min.diff)
                
                # J               = j
                # K               = k
                # max.L           = last.L
                # max.KL          = input$max.KL
                # psbl.seq        = input$psbl.seq
                # intcpt          = models$intcpt[models.i]
                # slp             = models$slp[models.i]
                # corstr          = input$CorrStr
                # corr.resid.err  = input$corr
                # sigma.resid.err = input$sigma.resid.err
                # var.intcpt      = input$var.intcpt
                # var.slp         = input$var.slp
                # cov.intcpt.slp  = input$cov.intcpt.slp
                # type.1.err      = input$type.1.err
                # power           = input$power
                # min.diff        = input$min.diff
                
                max.L[max.L.ind] <- TmpDsgn["L"]
                
                # If k the last one, then cannot assign to the next value
                # if the current value is not NA
                if ((max.L.ind != lth.max.L) & (!is.na(max.L[max.L.ind]))) {
                  
                  # assign if the next value is NA
                  if (is.na(max.L[max.L.ind + 1])) {
                    
                    # if the L * (k+1) goes out of gen.max.KL, next L go back to NA
                    # in other words, only assign when L * (k+1) <= gen.max.KL
                    if ((TmpDsgn["L"] * as.numeric(names(max.L)[max.L.ind + 1])) <= input$max.KL) {
                      max.L[max.L.ind + 1] <- TmpDsgn["L"]
                    }
                    
                  } else {
                    # if the next value is not NA, only assign if the value is greater than the current value
                    if (max.L[max.L.ind + 1] > TmpDsgn["L"]) {
                      max.L[max.L.ind + 1] <- TmpDsgn["L"]
                    }
                  }
                  
                }
                
                data.fig <- rbind(data.fig, 
                                  c(intcpt = as.character(models$intcpt[models.i]),
                                    slp    = as.character(models$slp[models.i]), 
                                    I      = I,
                                    TmpDsgn))
                
                k <- k + 1
                I <- calc_I_popavg(k, input$psbl.seq)
                
              } # while k loop
              
              # fill in all the k possibilities for later standard error calculations
              while ((k <= input$max.KL) & ((I * j) <= input$range.pts[2])) {
                
                data.fig <- rbind(data.fig,
                                  c(intcpt = as.character(models$intcpt[models.i]),
                                    slp    = as.character(models$slp[models.i]), 
                                    I      = I,
                                    J      = j,
                                    K      = k,
                                    L      = TmpDsgn["L"], 
                                    power  = NA)) 
                
                k <- k + 1
                I <- calc_I_popavg(k, input$psbl.seq)
              }
              
              # if all NAs or min(max.L) != 1, just assign max.L to prev.max.L
              if (sum(is.na(max.L)) == lth.max.L) {
                prev.max.L <- max.L
              } else if (min(max.L, na.rm = T) != 1) {
                prev.max.L <- max.L
              } else {
                # assign values until the first L = 1
                prev.max.L <- max.L[1:(which(max.L == 1)[1])]
              }
              
              # DO WE NEED TO DECIDE WHETHER THE REMAINING ONE ELEMENT IS 1?
              # Yes, although always shrinking length, but there may be prev.max.L length of 1 originally
              if (length(prev.max.L) == 1) {
                
                # fill in the missing J's and break out of the j loop
                if (prev.max.L == 1) {
                  
                  j <- j + 1
                  # k and i will not change because possible number of k drop to 1
                  # L drop to 1
                  I <- as.numeric(data.fig[nrow(data.fig), "I"])
                  K <- as.numeric(data.fig[nrow(data.fig), "K"])
                  L <- as.numeric(data.fig[nrow(data.fig), "L"])
                  
                  # fill in all the k possibilities for later standard error calculations
                  while ((I * j) <= input$range.pts[2]) {
                    
                    data.fig <- rbind(data.fig,
                                      c(intcpt = as.character(models$intcpt[models.i]),
                                        slp    = as.character(models$slp[models.i]), 
                                        I      = I,
                                        J      = j,
                                        K      = K,
                                        L      = L, 
                                        power  = NA)) 
                    
                    j <- j + 1
                  } 
                  
                  break
                } # if (prev.max.L == 1) {
                
              } # if (length(prev.max.L) == 1) {
              
              j <- j + 1
              
            } # while j loop
            
          } # for model.i loop
          
        } # else if input$psbl.seq != "other"
        
      }) # withProgress
      
      data.fig <- data.frame(data.fig)
      data.fig <- data.fig %>%
        filter(!is.na(L))
      # colnames(data.fig) <- c("K", "L", "J", "I", "Power", "CorrStr")
      # data.fig <- data.fig %>%
      #   mutate(I = as.integer(levels(I))[I],
      #          J = as.integer(levels(J))[J],
      #          K = as.integer(levels(K))[K],
      #          L = as.integer(levels(L))[L],
      #          power = as.numeric(levels(power))[power])
      data.fig <- data.fig %>%
        mutate(I = as.integer(I),
               J = as.integer(J),
               K = as.integer(K),
               L = as.integer(L),
               power = as.numeric(power))
      
      data.fig <- data.fig %>%
        filter(((K * L) <= input$max.KL) & 
                 ((I * J) <= input$range.pts[2]) &
                 ((I * J) >= input$range.pts[1]))
      data.fig <- data.fig %>%
        mutate(power = round(power, 2))
      data.fig   <- data.fig %>%
        mutate(model = ifelse(slp == "random", "Random slope model", "Common slope model"))
      data.fig <- data.fig %>%
        mutate(model = factor(model))
      
      # Optimize the designs
      # NEED TO MODIFY THIS FUNCTION FOR DIFFERENT COLUMNS IN DATA.FIG
      data.fig.opt <- optimize_dsgn(data.fig)
      # Summarize by desired x-axis
      v$summ.data.fig.opt <- gen_data_fig_pts(data.fig.opt)
      
      data.fig <- gen_data_fig_pts(data.fig)
      
    } # else (input$pop.avg.plot == "ttobs.vs.pts")
    
    v$data.fig   <- data.fig
    v$summ.data.fig.opt <- v$summ.data.fig.opt %>%
      mutate(model = ifelse(slp == "random", "Random slope model", "Common slope model"))
    v$summ.data.fig.opt <- v$summ.data.fig.opt %>%
      mutate(model = factor(model))
    # lth.model <- length(unique(summ.data.fig.opt$model))
    
    # generate a dataframe with one form of intercept
    v$data.fig.forFig <- v$data.fig
    v$summ.data.fig.opt.forFig <- v$summ.data.fig.opt
    if (sum(grepl("-random", input$model.popavg)) == 2) {
      v$data.fig.forFig <- v$data.fig.forFig %>%
        filter(!((intcpt == "fixed") & (slp == "random")))
      v$summ.data.fig.opt.forFig <- v$summ.data.fig.opt.forFig %>%
        filter(!((intcpt == "fixed") & (slp == "random")))
    }
    
    if (sum(grepl("-common", input$model.popavg)) == 2) {
      v$data.fig.forFig <- v$data.fig.forFig %>%
        filter(!((intcpt == "fixed") & (slp == "common")))
      v$summ.data.fig.opt.forFig <- v$summ.data.fig.opt.forFig %>%
        filter(!((intcpt == "fixed") & (slp == "common")))
    }
    
    # }  # if (input$trt.eff.intst == "pop.avg")
    
  })  # observeEvent
  
  # These parameters should show no matter which plot is produced
  output$params <- renderTable({
    
    if (is.null(v$summ.data.fig.opt.forFig)) {
      return()
    }
    
    psbl.seq   <- isolate(input$psbl.seq)
    type.1.err <- isolate(input$type.1.err)
    power      <- isolate(input$power)
    min.diff   <- isolate(input$min.diff)
    sigma.resid.err <- isolate(input$sigma.resid.err)
    CorrStr         <- isolate(input$CorrStr)
    
    corr            <- isolate(input$corr)
    
    model.popavg <- isolate(input$model.popavg)
    var.intcpt   <- isolate(input$var.intcpt)
    var.slp      <- isolate(input$var.slp)
    cov.intcpt.slp <- isolate(input$cov.intcpt.slp)
    
    
    params <- data.frame(Parameters = "Randomization scheme",
                         Values     = psbl.seq)
    
    if (psbl.seq == "User-Specified Sequences") {
      file.psbl.seq <- read.csv(input$file.psbl.seq$datapath)
      
      tmp.seqs <- apply(file.psbl.seq, 1, paste, collapse = "")
      tmp.seqs <- paste0(tmp.seqs, collapse = ", ")
      
      params <- rbind(params,
                      data.frame(Parameters = "Sequences",
                                 Values     = tmp.seqs))
    }
    
    params <- rbind(params,
                    data.frame(Parameters = c("Type I error",
                                              "Power",
                                              "Minimal clinically important treatment effect",
                                              "Homogeneous residual standard error",
                                              "Correlation structure"),
                               Values     = c(type.1.err,
                                              power,
                                              min.diff,
                                              sigma.resid.err,
                                              CorrStr)))
    
    if (CorrStr %in% c("Exchangeable", "AR-1")) {
      params <- rbind(params,
                      data.frame(Parameters = "Correlation", 
                                 Values     = as.character(corr)))
    }
    
    if (("random-common" %in% model.popavg) | ("random-random" %in% model.popavg)) {
      params <- rbind(params,
                      data.frame(Parameters = "Variance of random intercepts", 
                                 Values     = as.character(var.intcpt)))
    }
    
    if (("fixed-random" %in% model.popavg) | ("random-random" %in% model.popavg)) {
      params <- rbind(params,
                      data.frame(Parameters = "Variance of random slopes", 
                                 Values     = as.character(var.slp)))
    }
    
    if ("random-random" %in% model.popavg) {
      params <- rbind(params,
                      data.frame(Parameters = "Covariance of random intercepts and random slopes", 
                                 Values     = as.character(cov.intcpt.slp)))
    }
    
    params
    
  })
  
  ######################################################################
  # Plot for optimized designs for population average treatment effect #
  ######################################################################
  output$fig <- renderPlotly({
    
    if (is.null(v$summ.data.fig.opt.forFig)) {
      return()
    }
    
    pop.avg.plot <- isolate(input$pop.avg.plot)
    psbl.seq     <- isolate(input$psbl.seq)
    show.all.design <- isolate(input$show.all.design)
    
    if (pop.avg.plot == "ttobs.vs.ppobs") {
      
      if (show.all.design) {
        p <- ggplot(v$data.fig.forFig, aes(x = KL, y = mean.IJKL, group = model))
      } else {
        p <- ggplot(v$summ.data.fig.opt.forFig, aes(x = KL, y = mean.IJKL, group = model))
      }
      
      p <- p +
        # aes(color = CorrStr, shape = CorrStr), 
        geom_point(aes(color = model, shape = model, 
                       text = paste("Number of measurements per participant: ", KL,
                                    "<br>Average number of measurements across trials: ", round(mean.IJKL, 2))), 
                   size = 1.5) + 
        xlab("Number of measurements per participant") 
      
    } else {
      
      if (show.all.design) {
        # print(dim(v$data.fig.forFig))
        p <- ggplot(v$data.fig.forFig, aes(x = IJ, y = mean.IJKL, group = model))
      } else {
        p <- ggplot(v$summ.data.fig.opt.forFig, aes(x = IJ, y = mean.IJKL, group = model))
      }
      
      p <- p +
        geom_point(aes(color = model, 
                       shape = model, 
                       text = paste("Number of participants: ", IJ,
                                    "<br>Average number of measurements across trials: ", round(mean.IJKL, 2))), 
                   size = 1.5) + 
        xlab("Number of participants")
      
    }
    
    p <- p + 
      geom_line(aes(color = model, linetype = model)) + 
      geom_ribbon(aes(ymin = min.IJKL, ymax = max.IJKL, fill = model), alpha = 0.3) +
      ylab("Total number of measurements across trials") +
      ggtitle(ifelse(show.all.design,
                     paste0("All possible optimized designs within range of x-axis when ", 
                            psbl.seq, 
                            ifelse(psbl.seq == "User-Specified Sequences",
                                   " are",
                                   " is"), 
                            " used"), 
                     paste0("Optimized designs when ", 
                            psbl.seq, 
                            ifelse(psbl.seq == "User-Specified Sequences",
                                   " are",
                                   " is"), 
                            " used"))) +
      labs(color    = "Model form",
           shape    = "Model form",
           fill     = "Model form",
           linetype = "Model form") +
      theme_classic() +
      theme(axis.title   = element_text(size = 14),
            axis.text    = element_text(size = 12), 
            legend.title = element_blank(),
            legend.text  = element_text(size = 12))
    
    # plotly_json(ggplotly(p1))
    tmp <- ggplotly(p, tooltip = "text", source = "A") %>%
      # style(hoverinfo = "none", traces = (lth.model+1):(2*lth.model)) %>%
      layout(legend=list(title=list(text='<b> Model form </b>')))
    
    # Modify the ggplotly label compared with ggplot
    for (i in 1:length(tmp$x$data)) {
      tmp$x$data[[i]]$name <- str_match(tmp$x$data[[i]]$name, "\\((.*?)\\,")[, 2]
    }
    
    tmp
    
  })
  
  # Show this note if there are 2 models behind the lines
  output$plot.note <- renderUI({
    
    if (is.null(v$summ.data.fig.opt.forFig)) {
      return()
    }
    
    # use isolate to stop reacting to the changes in input
    model.popavg <- isolate(input$model.popavg)
    
    if ((("fixed-common" %in% model.popavg) & ("random-common" %in% model.popavg)) | 
        (("fixed-random" %in% model.popavg) & ("random-random" %in% model.popavg))) {
      HTML(paste0("*Plot sample size for one of the intercept form because form of slope determines the sample size."))
    }
    
  })
  
  ############################################################################
  # Only calculate after observing click, find actual designs & calculate se #
  ############################################################################
  click.fig.pop.avg <- reactive({
    req(v$summ.data.fig.opt.forFig)
    event_data("plotly_click", source = "A")
  })
  
  # summ.data.fig.opt should be extracted from earlier result
  # 1) find optimized designs
  # 2) calculate standard error for the corresponding x value
  observeEvent(click.fig.pop.avg(), {
    
    s <- click.fig.pop.avg()
    s <- data.frame(s)
    
    print(s)
    # s <- data.frame(x = 12, curveNumber = 1)
    # selected.scnr <- summ.data.fig.opt %>%
    #   filter((KL == 36) & (model == "Common slope model"))
    
    if ((v$calc.ct > 1) & (v$no.click)) {
      v$selected.scnr  <- NULL
      v$no.click       <- F
      
      v$summ.data.fig.se          <- NULL
      v$selected.summ.data.fig.se <- NULL
      v$df.selected.data.fig.se   <- NULL
      v$selected.data.fig.se.ind  <- NULL
      return()
    } 
    
    pop.avg.plot    <- isolate(input$pop.avg.plot)
    show.all.design <- isolate(input$show.all.design)
    psbl.seq        <- isolate(input$psbl.seq)
    CorrStr         <- isolate(input$CorrStr)
    corr            <- isolate(input$corr)
    sigma.resid.err <- isolate(input$sigma.resid.err)
    var.intcpt      <- isolate(input$var.intcpt)
    var.slp         <- isolate(input$var.slp)
    cov.intcpt.slp  <- isolate(input$cov.intcpt.slp)
    type.1.err      <- isolate(input$type.1.err)
    min.diff        <- isolate(input$min.diff)
    
    if (psbl.seq == "User-Specified Sequences") {
      file.psbl.seq <- read.csv(input$file.psbl.seq$datapath)
    } else {
      file.psbl.seq <- NULL    
    }
    
    if (pop.avg.plot == "ttobs.vs.ppobs") {
      
      if (show.all.design) {
        v$selected.scnr <- v$data.fig %>%
          filter((KL == s$x) & (as.numeric(model) == (s$curveNumber + 1)))
      } else {
        v$selected.scnr <- v$summ.data.fig.opt %>%
          filter((KL == s$x) & (as.numeric(model) == (s$curveNumber + 1)))
      }
      
      # print(v$selected.scnr)
      
      # if ((v$calc.ct > 1) & (v$no.click)) {
      #   v$selected.scnr  <- NULL
      #   v$no.click       <- F
      #   return()
      # } 
      # check if last.selected.scnr is has the same x but not the same y as the selected scnr; if so, return nothing
      # if (!is.null(v$clear.selected.scnr)) {
      # if ((v$last.selected.scnr$x %in% v$selected.scnr$KL) & (!(v$last.selected.scnr$y %in% v$selected.scnr$mean.IJKL))) {
      
      # If first time click calculate button, show selected scnr
      # if not, show selected scnr only on clicking the plot
      
    } else { # if (pop.avg.plot == "ttobs.vs.ppobs") {
      
      if (show.all.design) {
        v$selected.scnr <- v$data.fig %>%
          filter((IJ == s$x) & (as.numeric(model) == (s$curveNumber + 1)))
      } else {
        v$selected.scnr <- v$summ.data.fig.opt %>%
          filter((IJ == s$x) & (as.numeric(model) == (s$curveNumber + 1)))
      }
      
      # print(v$selected.scnr)
      
      # if ((v$calc.ct > 1) & (v$no.click)) {
      #   v$selected.scnr  <- NULL
      #   v$no.click       <- F
      #   return()
      # } 
      
    } # else if (pop.avg.plot == "ttobs.vs.ppobs") {
    
    v$selected.models <- v$selected.scnr[!duplicated(v$selected.scnr[, c("intcpt", "slp")]), c("intcpt", "slp")]
    # print(v$selected.models)
    
    # organize the tables
    v$df.psbl.scnr <- NULL
    for (i in 1:nrow(v$selected.scnr)) {
      
      if (is.na(v$selected.scnr$power[i])) {
        
        v$selected.scnr$power[i] <- calc_power(K = v$selected.scnr$K[i], 
                                               L = v$selected.scnr$L[i], 
                                               J = v$selected.scnr$J[i], 
                                               psbl.seq = psbl.seq, 
                                               intcpt   = v$selected.scnr$intcpt[i], 
                                               slp      = v$selected.scnr$slp[i],
                                               corstr   = CorrStr, 
                                               corr.resid.err  = corr, 
                                               sigma.resid.err = sigma.resid.err,
                                               var.intcpt      = var.intcpt, 
                                               var.slp         = var.slp, 
                                               cov.intcpt.slp  = cov.intcpt.slp,
                                               type.1.err      = type.1.err, 
                                               min.diff        = min.diff, 
                                               file.psbl.seq   = file.psbl.seq)
      } 
      
      v$df.psbl.scnr <- rbind(v$df.psbl.scnr,
                              v$selected.scnr[i, ])
      
    }
    
    v$df.psbl.scnr <- v$df.psbl.scnr %>%
      select(intcpt, I, J, K, L, IJKL, power) %>%
      mutate(intcpt = ifelse(intcpt == "fixed", "Fixed", as.character(intcpt))) %>%
      mutate(intcpt = ifelse(intcpt == "random", "Random", as.character(intcpt)))
    
    colnames(v$df.psbl.scnr) <- c("Form of intercept",
                                  "# of possible sequences",
                                  "# of participants per sequence",
                                  "# of periods",
                                  "# of measurements per period",
                                  "Total # of measurements across trials",
                                  "Power")
    
    # Calculate the standard errors 
    calc.indiv.prec <- isolate(input$calc.indiv.prec)
    
    if (calc.indiv.prec) {
      
      # v$summ.data.fig.se          <- NULL
      v$selected.summ.data.fig.se <- NULL
      v$df.selected.data.fig.se   <- NULL
      v$selected.data.fig.se.ind <- NULL
      
      # Naive estimate can be used for any model from for the population average treatment effect
      # shrunken estimate should be from the random slope model
      if (pop.avg.plot == "ttobs.vs.ppobs") {
        
        withProgress(message = "Calculating", value = 0, detail = "Naive Estimates", {
          
          max.IJ <- isolate(input$max.IJ)
          
          selected.data.fig <- v$data.fig %>%
            filter(((K*L) == s$x) & (as.numeric(model) == (s$curveNumber + 1)))
          
          # naive estimates will apply to both common slope and random slope model
          # satisfy.niv should include all the k possibilities for I * J because the iterations goes through all j and all possible k within j
          satisfy.niv <- selected.data.fig %>%
            arrange(K, L) %>%
            select(K, L)
          satisfy.niv <- satisfy.niv[!duplicated(satisfy.niv), ]
          
          data.fig.se <- NULL
          
          # Calculate naive se
          tmp.data.fig.niv.se <- NULL
          for (i in 1:nrow(satisfy.niv)) {
            
            k <- satisfy.niv$K[i]
            l <- satisfy.niv$L[i]
            
            tmp.se <- calc_se_niv(K               = k,
                                  L               = l,
                                  psbl.seq        = psbl.seq,
                                  sigma.resid.err = sigma.resid.err,
                                  corstr          = CorrStr,
                                  corr.resid.err  = corr,
                                  file.psbl.seq   = file.psbl.seq)
            
            # K               = k
            # L               = l
            # psbl.seq        = psbl.seq
            # sigma.resid.err = sigma.resid.err
            # corstr          = CorrStr
            # corr.resid.err  = corr
            # file.psbl.seq   = file.psbl.seq
            
            tmp.data.fig.niv.se <- rbind(tmp.data.fig.niv.se,
                                         data.frame(K = k,
                                                    L = l,
                                                    tmp.se))
          }
          
          satisfy.niv.ijkl <- selected.data.fig %>%
            arrange(I, J) %>%
            select(I, J, K, L)
          satisfy.niv.ijkl <- satisfy.niv.ijkl[!duplicated(satisfy.niv.ijkl[, "I"]), ]
          
          for (row.i in 1:nrow(satisfy.niv.ijkl)) {
            
            # for fixed i, k will be the same; because KL is fixed, l will be the same
            i <- satisfy.niv.ijkl$I[row.i]
            k <- satisfy.niv.ijkl$K[row.i]
            l <- satisfy.niv.ijkl$L[row.i]
            
            tmp.se <- tmp.data.fig.niv.se[(tmp.data.fig.niv.se$K == k) & (tmp.data.fig.niv.se$L == l), ]
            
            for (j in (satisfy.niv.ijkl$J[row.i]):floor(max.IJ / i)) {
              
              data.fig.se <- rbind(data.fig.se,
                                   data.frame(Estimator = "Naive Estimates",
                                              I         = i,
                                              J         = j,
                                              tmp.se))
              
            } # for j loop
          } #   for (row.i in 1:nrow(satisfy.niv.ijkl)) {
          
          print("Progress in calculating:")
          
          if (unique(selected.data.fig$model) == "Random slope model") {
            
            for (row.i in 1:nrow(selected.data.fig)) {
              
              i <- selected.data.fig$I[row.i]
              k <- selected.data.fig$K[row.i]
              l <- selected.data.fig$L[row.i]
              
              incProgress(1/2 + (row.i - 1) / nrow(selected.data.fig) / 2,
                          detail = paste0(ifelse(selected.data.fig$intcpt[row.i] == "fixed",
                                                 "Shrunken Estimates-Fixed Intercepts",
                                                 "Shrunken Estimates-Random intercepts"), ": K = ", k))
              print(1/2 + (row.i - 1) / nrow(selected.data.fig) / 2)
              
              for (j in (selected.data.fig$J[row.i]):floor(max.IJ / i)) {
                
                # the first column in tmp.se is indexing which participant assigned to the specific sequence
                tmp.se <- calc_se_shrk(K = k,
                                       L = l,
                                       J = j,
                                       psbl.seq = psbl.seq,
                                       intcpt   = selected.data.fig$intcpt[row.i],
                                       sigma.resid.err = sigma.resid.err,
                                       corstr          = CorrStr,
                                       corr.resid.err  = corr,
                                       var.intcpt      = var.intcpt,
                                       var.slp         = var.slp,
                                       cov.intcpt.slp  = cov.intcpt.slp,
                                       var.rand.eff    = T,
                                       file.psbl.seq   = file.psbl.seq)
                
                # K = k
                # L = l
                # J = j
                # psbl.seq = psbl.seq
                # intcpt   = selected.data.fig$intcpt[row.i]
                # sigma.resid.err = sigma.resid.err
                # corstr          = CorrStr
                # corr.resid.err  = corr
                # var.intcpt      = var.intcpt
                # var.slp         = var.slp
                # cov.intcpt.slp  = cov.intcpt.slp
                # var.rand.eff    = T
                # file.psbl.seq   = file.psbl.seq
                
                # if (length(unique(tmp.se[, 1])) > 1) {
                #   print(selected.data.fig[row.i, ])
                #   stop("People assigned to the same sequence have different se")
                # }
                
                data.fig.se <- rbind(data.fig.se,
                                     data.frame(Estimator = ifelse(selected.data.fig$intcpt[row.i] == "fixed",
                                                                   "Shrunken Estimates-Fixed Intercepts",
                                                                   "Shrunken Estimates-Random intercepts"),
                                                I = i,
                                                J = j,
                                                K = k,
                                                L = l,
                                                tmp.se %>% select(seq, se)))
                
              } # for j loop
            } # for row.i loop
          } # if (unique(selected.data.fig$model) == "Random slope model") {
        }) # withProgress
        
        # data.fig.se <- data.frame(data.fig.se)
        # data.fig.se <- data.fig.se %>%
        #   mutate(I = as.integer(levels(I))[I],
        #          J = as.integer(levels(J))[J],
        #          K = as.integer(levels(K))[K],
        #          L = as.integer(levels(L))[L],
        #          se = as.numeric(levels(se))[se])
        # data.fig.se <- data.fig.se %>%
        #   mutate(I = as.integer(I),
        #          J = as.integer(J),
        #          K = as.integer(K),
        #          L = as.integer(L),
        #          se = as.numeric(se))
        data.fig.se <- data.fig.se %>%
          mutate(Estimator = factor(Estimator))
        
        summ.data.fig.se <- data.fig.se %>%
          mutate(IJ = I * J) %>%
          group_by(Estimator, IJ) %>%
          summarise(avg_se = mean(se),
                    max_se = max(se),
                    min_se = min(se),
                    n      = n())
        summ.data.fig.se <- data.frame(summ.data.fig.se)
        
      } else { #  if (pop.avg.plot == "ttobs.vs.ppobs") {
        
        withProgress(message = "Calculating", value = 0, detail = "Naive Estimates", {
          
          max.KL <- isolate(input$max.KL)
          
          selected.data.fig <- v$data.fig %>%
            filter(((I*J) == s$x) & (as.numeric(model) == (s$curveNumber + 1)))
          
          # naive estimates will apply to both common slope and random slope model
          # satisfy.niv should include all the k possibilities for I * J because the iterations goes through all j and all possible k within j
          satisfy.niv <- selected.data.fig %>%
            arrange(K, L) %>%
            select(I, J, K, L)
          satisfy.niv <- satisfy.niv[!duplicated(satisfy.niv[, "K"]), ]
          
          data.fig.se <- NULL
          # Calculate naive se
          for (row.i in 1:nrow(satisfy.niv)) {
            
            # for fixed k, i will be the same; because IJ is fixed, j will be the same
            k <- satisfy.niv$K[row.i]
            i <- satisfy.niv$I[row.i]
            j <- satisfy.niv$J[row.i]
            
            for (l in (satisfy.niv$L[row.i]):floor(max.KL / k)) {
              # for (l in 51:52) {
              tmp.se <- calc_se_niv(K               = k, 
                                    L               = l, 
                                    psbl.seq        = psbl.seq, 
                                    sigma.resid.err = sigma.resid.err,
                                    corstr          = CorrStr,
                                    corr.resid.err  = corr,
                                    file.psbl.seq   = file.psbl.seq)
              
              data.fig.se <- rbind(data.fig.se,
                                   cbind(Estimator = "Naive Estimates",
                                         I         = i,
                                         J         = j,
                                         K         = k, 
                                         L         = l, 
                                         tmp.se))
            } # for l loop
            
          } # for row.i loop satisfy.niv
          
          if (unique(selected.data.fig$model) == "Random slope model") {
            
            for (row.i in 1:nrow(selected.data.fig)) {
              
              k <- selected.data.fig$K[row.i]
              i <- selected.data.fig$I[row.i]
              j <- selected.data.fig$J[row.i]
              
              incProgress(1/2 + (row.i - 1) / nrow(selected.data.fig) / 2,
                          detail = paste0(ifelse(selected.data.fig$intcpt[row.i] == "fixed",
                                                 "Shrunken Estimates-Fixed Intercepts",
                                                 "Shrunken Estimates-Random intercepts"), ": I = ", i))
              print(1/2 + (row.i - 1) / nrow(selected.data.fig) / 2)
              
              for (l in (selected.data.fig$L[row.i]):floor(max.KL / k)) {
                # for (l in 51:52) {
                tmp.se <- calc_se_shrk(K = k, 
                                       L = l, 
                                       J = j, 
                                       psbl.seq = psbl.seq, 
                                       intcpt   = selected.data.fig$intcpt[row.i], 
                                       sigma.resid.err = sigma.resid.err,
                                       corstr          = CorrStr,
                                       corr.resid.err  = corr,
                                       var.intcpt      = var.intcpt, 
                                       var.slp         = var.slp, 
                                       cov.intcpt.slp  = cov.intcpt.slp, 
                                       var.rand.eff    = T,
                                       file.psbl.seq   = file.psbl.seq)
                if (length(unique(tmp.se[, 1])) > 1) {
                  print(selected.data.fig[i, ])
                  stop("People assigned to the same sequence have different se")
                }
                
                if (nrow(tmp.se) == 1) {
                  data.fig.se <- rbind(data.fig.se,
                                       c(Estimator = ifelse(selected.data.fig$intcpt[row.i] == "fixed", 
                                                            "Shrunken Estimates-Fixed Intercepts", 
                                                            "Shrunken Estimates-Random intercepts"), 
                                         I = i,
                                         J = j,
                                         K = k, 
                                         L = l, 
                                         tmp.se[2:3]))
                } else {
                  data.fig.se <- rbind(data.fig.se,
                                       cbind(Estimator = ifelse(selected.data.fig$intcpt[row.i] == "fixed", 
                                                                "Shrunken Estimates-Fixed Intercepts", 
                                                                "Shrunken Estimates-Random intercepts"), 
                                             I = i,
                                             J = j,
                                             K = k, 
                                             L = l, 
                                             tmp.se[, 2:3]))
                }
                
              } # for l loop
            } # for i loop selected.data.fig for shrunken estimates
          } # random slope model
        }) # withProgress
        
        data.fig.se <- data.frame(data.fig.se)
        # data.fig.se <- data.fig.se %>%
        #   mutate(I = as.integer(levels(I))[I],
        #          J = as.integer(levels(J))[J],
        #          K = as.integer(levels(K))[K],
        #          L = as.integer(levels(L))[L],
        #          se = as.numeric(levels(se))[se])
        data.fig.se <- data.fig.se %>%
          mutate(I = as.integer(I),
                 J = as.integer(J),
                 K = as.integer(K),
                 L = as.integer(L),
                 se = as.numeric(se))
        data.fig.se <- data.fig.se %>%
          mutate(Estimator = factor(Estimator))
        
        summ.data.fig.se <- data.fig.se %>%
          mutate(KL = K * L) %>%
          group_by(Estimator, KL) %>%
          summarise(avg_se = mean(se), 
                    max_se = max(se), 
                    min_se = min(se),
                    n      = n())
        summ.data.fig.se <- data.frame(summ.data.fig.se)
        
      } # if ttobs.vs.ppobs standard error calculation
      
      v$data.fig.se      <- data.fig.se
      v$summ.data.fig.se <- summ.data.fig.se
      
      # js$resetClick()
      
    } # if calc.indiv.prec
    
  }) # observeEvent for "plotly_click"
  
  ####################################
  # Show optimized design basic info #
  ####################################
  output$click.caption <- renderUI({
    
    # when there is nothing to show, return ""
    if ((is.null(v$summ.data.fig.opt)) | (is.null(v$selected.scnr))) {
      HTML("")
      return()
    }
    
    pop.avg.plot <- isolate(input$pop.avg.plot)
    
    if (pop.avg.plot == "ttobs.vs.ppobs") {
      
      text.caption <- paste0("<em>You have fixed the number of measurements per participant at <b>", unique(v$selected.scnr$KL), "</b>.</em><br/>")
      
      for (selected.models.i in 1:nrow(v$selected.models)) {
        text.caption <- paste0(text.caption,
                               "Under <b>", v$selected.models$intcpt[selected.models.i], "</b> intercept-<b>", v$selected.models$slp[selected.models.i], "</b> slope model, ",
                               "when the number of measurements per parcipant is fixed at <b>", unique(v$selected.scnr$KL), "</b>,",
                               "<ul><li>the minimum total number of measurements across trials for an optimized design is <b>", 
                               unique(v$selected.scnr$min.IJKL[v$selected.scnr$intcpt == v$selected.models$intcpt[selected.models.i]]), "</b>;",
                               "<li>the maximum total number of measurements across trials for an optimized design is <b>", unique(v$selected.scnr$max.IJKL[v$selected.scnr$intcpt == v$selected.models$intcpt[selected.models.i]]),
                               "</b>. </li></ul>")
      }
      
    } else { # if (input$pop.avg.plot == "ttobs.vs.ppobs") {
      
      text.caption <- paste0("<em>You have fixed the number of participants at <b>", unique(v$selected.scnr$IJ), "</b>.</em><br/>")
      for (selected.models.i in 1:nrow(v$selected.models)) {
        text.caption <- paste0(text.caption,
                               "Under <b>", v$selected.models$intcpt[selected.models.i], "</b> intercept-<b>", v$selected.models$slp[selected.models.i], "</b> slope model, ",
                               "when the number of parcipants is fixed at <b>", unique(v$selected.scnr$IJ), "</b>,",
                               "<ul><li>the minimum total number of measurements across trials for an optimized design is <b>", 
                               unique(v$selected.scnr$min.IJKL[v$selected.scnr$intcpt == v$selected.models$intcpt[selected.models.i]]), "</b>;",
                               "<li>the maximum total number of measurements across trials for an optimized design is <b>", unique(v$selected.scnr$max.IJKL[v$selected.scnr$intcpt == v$selected.models$intcpt[selected.models.i]]),
                               "</b>. </li></ul>")
      }
      
    } # if (input$pop.avg.plot == "ttobs.vs.ppobs") {
    
    
    text.caption <- paste0(text.caption,
                           "The possible designs are:")
    HTML(text.caption)
    
  })
  
  ####################################
  ## Show optimized design in table ##
  ####################################
  output$psbl.scnr <- renderTable({
    
    if (is.null(v$summ.data.fig.opt) | is.null(v$selected.scnr)) {
      return()
    }
    
    v$df.psbl.scnr
  }, digits = 2)
  
  output$click.caption.2 <- renderUI({
    
    if (is.null(v$summ.data.fig.se)) {
      return()
    }
    
    text.caption <- paste0("<em>You have selected to use <b>", ifelse(unique(v$selected.models$slp) == "random",  "Random slope model", "Common slope model"), 
                           "</b> to estimate the population average treatment effect. </em><br/>Individual-specific treatment effect can then be estimated using<ul>")
    
    lvl.indiv.est <- levels(v$summ.data.fig.se$Estimator)
    for (i in 1:length(lvl.indiv.est)) {
      text.caption <- paste0(text.caption,
                             "<li>", lvl.indiv.est[i])
    }
    
    text.caption <- paste0(text.caption, 
                           "</li></ul>")
    
    HTML(text.caption)
  })
  
  
  ######################################################################
  ############### Plot se for selected optimized design ################
  ######################################################################
  output$fig.se <- renderPlotly({
    
    if (is.null(v$summ.data.fig.se)) {
      return()
    }
    
    pop.avg.plot  <- isolate(input$pop.avg.plot)
    selected.scnr <- isolate(v$selected.scnr)
    
    
    if (pop.avg.plot == "ttobs.vs.ppobs") {
      p <- ggplot(v$summ.data.fig.se, aes(x = IJ, y = avg_se, group = Estimator)) +
        # aes(color = CorrStr, shape = CorrStr),
        geom_point(aes(color = Estimator, shape = Estimator,
                       text = paste("Number of participants: ", IJ,
                                    "<br>Average standard error: ", round(avg_se, 2))),
                   size = 1) +
        xlab("Number of participants") +
        ggtitle(paste0("When the number of measurements per participant is fixed at ", unique(selected.scnr$KL)))
      
    } else {
      
      p <- ggplot(v$summ.data.fig.se, aes(x = KL, y = avg_se, group = Estimator)) +
        geom_point(aes(color = Estimator, shape = Estimator,
                       text = paste("Number of measurements per participant: ", KL,
                                    "<br>Average standard error: ", round(avg_se, 2))),
                   size = 1) +
        xlab("Number of measurements per participant") +
        ggtitle(paste0("When the number of participants is fixed at ", unique(selected.scnr$IJ)))
      
    }
    
    p <- p +
      geom_line(aes(color = Estimator, linetype = Estimator)) +
      geom_ribbon(aes(ymin = min_se, ymax = max_se, fill = Estimator), alpha = 0.3) +
      ylab("Standard error of individual-specific treatment effect") +
      labs(color    = "Estimator",
           shape    = "Estimator",
           fill     = "Estimator",
           linetype = "Estimator") +
      theme_classic() +
      theme(axis.title   = element_text(size = 14),
            axis.text    = element_text(size = 12),
            legend.title = element_blank(),
            legend.text  = element_text(size = 12))
    
    # plotly_json(ggplotly(p1))
    tmp <- ggplotly(p, tooltip = "text", source = "B", type = "scattergl") %>%
      # style(hoverinfo = "none", traces = (lth.model+1):(2*lth.model)) %>%
      layout(legend=list(title=list(text='<b> Estimator </b>')))
    
    # Modify the ggplotly label compared with ggplot
    for (i in 1:length(tmp$x$data)) {
      tmp$x$data[[i]]$name <- str_match(tmp$x$data[[i]]$name, "\\((.*?)\\,")[, 2]
    }
    
    tmp
    
  })
  
  #######################################################################
  # Only calculate after observing click, find se at given x axis value #
  #######################################################################
  click.fig.se.ovrl <- reactive({
    # req(v$data.fig.se)
    event_data("plotly_click", source = "B")
  })
  
  # se plot click
  observeEvent(click.fig.se.ovrl(), {
    
    s2 <- click.fig.se.ovrl()
    # print(s)
    # if (is.null(s)) {
    #   v$selected.summ.data.fig.se <- NULL
    # }
    s2 <- data.frame(s2)
    
    print(s2)
    # s2 <- data.frame(x = 48, curveNumber = 1)
    # selected.scnr <- summ.data.fig.opt %>%
    #   filter((KL == 36) & (model == "Common slope model"))
    
    pop.avg.plot <- isolate(input$pop.avg.plot)
    psbl.seq     <- isolate(input$psbl.seq)
    
    if (pop.avg.plot == "ttobs.vs.ppobs") {
      selected.data.fig.se <- v$data.fig.se %>%
        filter(((I * J) == s2$x) & (as.numeric(Estimator) == (s2$curveNumber + 1)))
      
      selected.summ.data.fig.se <- selected.data.fig.se %>%
        group_by(Estimator, I) %>%
        summarise(avg_se = mean(se),
                  max_se = max(se),
                  min_se = min(se),
                  n      = n())
      
      selected.data.fig.se <- selected.data.fig.se %>%
        left_join(selected.summ.data.fig.se, by = c("Estimator", "I"))
      
    } else {
      selected.data.fig.se <- v$data.fig.se %>%
        filter(((K * L) == s2$x) & (as.numeric(Estimator) == (s2$curveNumber + 1)))
      
      selected.summ.data.fig.se <- selected.data.fig.se %>%
        group_by(Estimator, K) %>%
        summarise(avg_se = mean(se),
                  max_se = max(se),
                  min_se = min(se),
                  n      = n())
      
      selected.data.fig.se <- selected.data.fig.se %>%
        left_join(selected.summ.data.fig.se, by = c("Estimator", "K"))
      
    }
    
    # selected.data.fig.se <- selected.data.fig.se %>%
    #   mutate(seq = as.character(seq))
    # if (psbl.seq == "Pairwise Randomization") {
    #   for (i in 1:nrow(selected.data.fig.se)) {
    #     selected.data.fig.se$seq[i] <- gen.chr.seq(seq = as.character(selected.data.fig.se$seq[i]),
    #                                                psbl.seq = psbl.seq,
    #                                                K   = selected.data.fig.se$K[i])
    #   }
    # }
    
    df.selected.data.fig.se <- selected.data.fig.se %>%
      select(I, J, K, L, seq, se)
    
    colnames(df.selected.data.fig.se) <- 
      # c("Estimator for individual-specific treatment effect",
      c("# of possible sequences",
        "# of participants per sequence",
        "# of periods",
        "# of measurements per period",
        "Treatment order",
        "s.e.")
    
    v$selected.data.fig.se      <- data.frame(selected.data.fig.se)
    v$selected.summ.data.fig.se <- data.frame(selected.summ.data.fig.se)
    v$df.selected.data.fig.se   <- df.selected.data.fig.se
  })
  
  ###################################################
  # Show se for designs that satisfy both IJ and KL #
  ###################################################
  output$click.caption.se.ovrl <- renderUI({
    
    if (is.null(v$df.selected.data.fig.se)) {
      HTML("")
      return()
    }
    
    pop.avg.plot <- isolate(input$pop.avg.plot)
    
    if (pop.avg.plot == "ttobs.vs.ppobs") {
      
      text.caption <- paste0("<em>You have fixed the number of participants at <b>", unique(v$selected.data.fig.se$I * v$selected.data.fig.se$J),
                             "</b> and the number of measurements per participant at <b>", unique(v$selected.scnr$KL), "</b>.</em>",
                             "<ul><li>The total number of measurements across trials is ", unique(v$selected.data.fig.se$I * v$selected.data.fig.se$J) * unique(v$selected.scnr$KL), ".</li></ul>")
      
    } else {
      text.caption <- paste0("<em>You have fixed the number of participants at <b>", unique(v$selected.scnr$IJ),
                             "</b> and the number of measurements per participant at <b>", unique(v$selected.data.fig.se$K * v$selected.data.fig.se$L), "</b>.</em>",
                             "<ul><li>The total number of measurements across trials is ", unique(v$selected.data.fig.se$K * v$selected.data.fig.se$L) * unique(v$selected.scnr$IJ), ".</li></ul>")
    }
    
    text.caption <- paste0(text.caption,
                           "<em>You have selected to use <b>", ifelse(unique(v$selected.models$slp) == "random",  "Random slope model", "Common slope model"), 
                           "</b> to estimate the population average treatment effect and <b>", unique(v$selected.data.fig.se$Estimator), 
                           "</b> to estimate the individual-specific treamtent effects.</em><br/><br/>")
    
    text.caption <- paste0(text.caption,
                           "The standard errors for the individual-specific treatment effects for all the possible designs satisfying the power requirements are:")
    HTML(text.caption)
    
  })
  
  output$df.selected.data.fig.se <- renderTable({
    
    if (is.null(v$df.selected.data.fig.se)) {
      return()
    }
    
    v$df.selected.data.fig.se
    
  }, digits = 2)
  
  # ############################
  # # se individual level plot #
  # ############################
  # output$fig.se.ind <- renderPlotly({
  #   
  #   if (is.null(v$selected.summ.data.fig.se)) {
  #     return()
  #   } 
  #   
  #   pop.avg.plot         <- isolate(input$pop.avg.plot)
  #   # selected.scnr        <- isolate(v$selected.scnr)
  #   # selected.data.fig.se <- isolate(v$selected.data.fig.se)
  #   
  #   # print(selected.data.fig.se)
  #   
  #   if (pop.avg.plot == "ttobs.vs.ppobs") {
  #     
  #     if (length(unique(v$selected.summ.data.fig.se$I)) == 1) {
  #       return()
  #     }
  #     
  #     p <- ggplot(v$selected.summ.data.fig.se, aes(x = I, y = avg_se, group = Estimator)) +
  #       # aes(color = CorrStr, shape = CorrStr),
  #       geom_point(aes(color = Estimator, shape = Estimator,
  #                      text = paste("Number of sequences: ", I,
  #                                   "<br>Average standard error: ", round(avg_se, 2))),
  #                  size = 1) +
  #       xlab("Number of sequences") +
  #       ggtitle(paste0("When the number of participants is fixed at ", unique(v$selected.data.fig.se$I * v$selected.data.fig.se$J), 
  #                      " and the number of measurements per participant is fixed at ", unique(v$selected.scnr$KL)))
  #     
  #   } else {
  #     
  #     if (length(unique(v$selected.summ.data.fig.se$K)) == 1) {
  #       return()
  #     }
  #     
  #     p <- ggplot(v$selected.summ.data.fig.se, aes(x = K, y = avg_se, group = Estimator)) +
  #       geom_point(aes(color = Estimator, shape = Estimator,
  #                      text = paste("Number of periods per participant: ", K,
  #                                   "<br>Average standard error: ", round(avg_se, 2))),
  #                  size = 1) +
  #       xlab("Number of periods per participant") +
  #       ggtitle(paste0("When the number of participants is fixed at ", unique(v$selected.scnr$IJ), 
  #                      " and the number of measurements per participant is fixed at ", unique(v$selected.data.fig.se$K * v$selected.data.fig.se$L)))
  #     
  #   }
  #   
  #   p <- p +
  #     geom_line(aes(color = Estimator, linetype = Estimator)) +
  #     geom_ribbon(aes(ymin = min_se, ymax = max_se, fill = Estimator), alpha = 0.3) +
  #     ylab("Standard error of individual-specific treatment effect") +
  #     labs(color    = "Estimator",
  #          shape    = "Estimator",
  #          fill     = "Estimator",
  #          linetype = "Estimator") +
  #     theme_classic() +
  #     theme(axis.title   = element_text(size = 14),
  #           axis.text    = element_text(size = 12),
  #           legend.title = element_blank(),
  #           legend.text  = element_text(size = 12))
  #   
  #   # plotly_json(ggplotly(p1))
  #   tmp <- ggplotly(p, tooltip = "text", source = "C") %>%
  #     # style(hoverinfo = "none", traces = (lth.model+1):(2*lth.model)) %>%
  #     layout(legend=list(title=list(text='<b> Estimator </b>')))
  #   
  #   # Modify the ggplotly label compared with ggplot
  #   for (i in 1:length(tmp$x$data)) {
  #     tmp$x$data[[i]]$name <- str_match(tmp$x$data[[i]]$name, "\\((.*?)\\,")[, 2]
  #   }
  #   
  #   tmp
  #   
  # })
  # 
  # ##########################################################################################################
  # # Only calculate after observing click, find actual design that satisfy power and precision requirements #
  # ##########################################################################################################
  # click.fig.se.ind <- reactive({
  #   # req(v$selected.summ.data.fig.se)
  #   event_data("plotly_click", source = "C")
  # })
  # 
  # observeEvent(click.fig.se.ind(), {
  #   
  #   s3 <- click.fig.se.ind()
  #   s3 <- data.frame(s3)
  #   
  #   # print(s)
  #   # s <- data.frame(x = 16, curveNumber = 2)
  #   # selected.scnr <- summ.data.fig.opt %>%
  #   #   filter((KL == 36) & (model == "Common slope model"))
  #   pop.avg.plot <- isolate(input$pop.avg.plot)
  #   
  #   if (pop.avg.plot == "ttobs.vs.ppobs") {
  #     selected.data.fig.se.ind <- v$selected.data.fig.se %>%
  #       filter((I == s3$x) & (as.numeric(Estimator) == (s3$curveNumber + 1)))
  #     
  #   } else {
  #     selected.data.fig.se.ind <- v$selected.data.fig.se %>%
  #       filter((K == s3$x) & (as.numeric(Estimator) == (s3$curveNumber + 1)))
  #   }
  #   
  #   selected.data.fig.se.ind <- selected.data.fig.se.ind %>%
  #     select(Estimator, I, J, K, L, seq, se)
  #   
  #   colnames(selected.data.fig.se.ind) <- c("Estimator",
  #                                           "# of sequences",
  #                                           "# of participants per sequence",
  #                                           "# of periods per participant",
  #                                           "# of measurements per period",
  #                                           "Treatment order", 
  #                                           "s.e.")
  #   
  #   v$selected.data.fig.se.ind <- selected.data.fig.se.ind
  #   
  # })
  # 
  # output$df.se <- renderTable({
  #   
  #   if (is.null(v$selected.data.fig.se.ind)) {
  #     return()
  #   }
  #   
  #   pop.avg.plot         <- isolate(input$pop.avg.plot)
  #   
  #   if (pop.avg.plot == "ttobs.vs.ppobs") {
  #     
  #     if (length(unique(v$selected.summ.data.fig.se$I)) == 1) {
  #       return()
  #     }
  #     
  #   } else {
  #     
  #     if (length(unique(v$selected.summ.data.fig.se$K)) == 1) {
  #       return()
  #     }
  #     
  #   }
  #   
  #   v$selected.data.fig.se.ind
  # }, digits = 2)
  
  
}


################################### Run the app ###################################
shinyApp(ui = ui, server = server)