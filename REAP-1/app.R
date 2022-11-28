#

library(shiny)
library(shinyBS)
library(gridExtra)
library(msm)
library(stats)
library(ggplot2)
library(betareg)
library(boot)
library(reshape2)
library(xtable)
library(Metrics)
library(memisc)
library(DT)
library(kableExtra)
library(Rmisc)
library(plyr)
library(shinyjs)
library(fastDummies)
library(lmtest)
library(rpsychi)
library(dplyr)
library(ggprism)
library(ggnewscale)
library(scales)
# library(nloptr)
library(robustbase)
library(matlib)
library(BBmisc)
library(shinybusy)

pairwise_compare <- function(m, sd, n){
  k <- length(m)
  Xg <- sum(n * m)/sum(n)
  dfb <- k - 1
  dfw <- sum(n) - k
  MSb <- sum(n * (m - Xg)^2)/(k - 1)
  MSw <- sum((n - 1) * sd^2)/dfw
  SSb <- dfb * MSb
  SSw <- dfw * MSw
  SSt <- SSb + SSw
  f.value <- MSb/MSw
  anova.table <- data.frame(matrix(NA, ncol = 4, nrow = 3))
  rownames(anova.table) <- c("Between (A)", "Within", "Total")
  colnames(anova.table) <- c("SS", "df", "MS", "F")
  anova.table$SS <- c(SSb, SSw, SSt)
  anova.table$df <- c(dfb, dfw, dfb + dfw)
  anova.table$MS <- c(MSb, MSw, NA)
  anova.table$F <- c(f.value, NA, NA)
  class(anova.table) <- c("anova", "data.frame")
  return(anova.table)
}

displaytable = read.table("31780660_F1B_exampledata.csv", sep = ",", header = TRUE, check.names = FALSE)
displaytable = head(displaytable)

source("MDPDE2.R", local = TRUE)


# Define UI for application that draws a histogram
ui <- pageWithSidebar(
  
  
  # Application title
  headerPanel("Robust and Efficient Assessment of drug Potency (REAP)"),
  
  # Sidebar with a slider input for number of bins
  
  
  sidebarPanel(
    
    # Progress bar
    add_busy_spinner(spin = "fading-circle"),
    
    
    
    # fixedPanel(
    tags$style(type = "text/css", '#intro_panel {position : fixed; width: 20%}'),
    id = "intro_panel",
    conditionalPanel(condition="input.tabset==1",
                     h4("Outline"),
                     h5(a("Introduction", href = "#intro")),
                     h5(a("Example", href = "#exp")),
                     h5(a("Reference", href = "#ref"))
                     
    ),
    # a("TEST", href= "#tabset")),
    
    conditionalPanel(condition="input.tabset==2",
                     # File input
                     fileInput(
                       inputId = "filedata",
                       "Upload data. Choose csv file",
                       accept = c(".csv")
                     )
                     
                     
    ),
    
    conditionalPanel(condition="input.tabset==3",
                     h3("Model Feature"),
                     checkboxInput("checkbox_betareg",
                                   label = "Enable dose-dependent precision",
                                   value = FALSE),
                     numericInput(inputId = "effectpct",
                                  tags$span("Add effect estimation",
                                            tags$i(
                                              id = "effectpct",
                                              class = "glyphicon glyphicon-info-sign",
                                              style = "color:#0072B2; ",
                                              title = "Input percentage of effect, from 0 to 100; Effect estimates are shown with triangles in the dose-response curve plot"  
                                            )), 
                                  value=50, min=0, max=100),
                     hr(),
                     # h3("Model Comparisons"),
                     tags$h3(
                       tags$span(
                         "Model Comparisons",
                         tags$i(
                           id = "modelcompare",
                           class = "glyphicon glyphicon-info-sign",
                           style = "color:#0072B2; ",
                           title = "for effect estimations or slope comparison if only one box was checked, or the curve comparison if both boxes were checked."  
                         )
                       )
                     ),
                     
                     checkboxInput("checkbox4",
                                   label = tags$span("Effect estimations"),
                                   value = FALSE),
                     checkboxInput("checkbox5",
                                   label = tags$span("Slopes"),
                                   value = FALSE),
                     
                     
                     hr(),
                     h3("Plot Specifics"),
                     checkboxInput("checkbox1", 
                                   label = tags$span("Remove log-transformation of x-axis"),
                                   value = FALSE),
                     
                     checkboxInput("checkbox2", 
                                   label = tags$span("Add sample mean and sample S.D."),
                                   value = FALSE),
                     useShinyjs(),
                     shinyjs::hidden(
                       # numericInput(inputId = "num_width", "Width of CI", value = 0.02, min = 0, max = 0.1)
                       numericInput(inputId = "num_width",
                                    tags$span("Width of Error Bar",
                                              tags$i(
                                                id = "num_width1",
                                                class = "glyphicon glyphicon-info-sign",
                                                style = "color:#0072B2; ",
                                                title = "Input width of error bar, from 0 to 0.1"
                                              )),
                                    value=0.02, min=0, max=0.1)
                       
                     ),
                     
                     span(textOutput("CIinput"), style="color:red"),
                     
                     hr(),
                     # h4("Downloaded Plot"),
                     tags$h3(
                       tags$span(
                         # tags$h4("Download Plot"),
                         "Download Plot",
                         tags$i(
                           id = "icon",
                           class = "glyphicon glyphicon-info-sign",
                           style = "color:#0072B2; ",
                           title = "Specify the width and height of the downloaded plot"  
                         )
                       )
                     ),
                     fixedRow(
                       column(6, numericInput("graphwidth", "Width",NA, value = NULL, min = 0.000001)),
                       column(6, numericInput("graphheight", "Height",NA, value = NULL, min = 0.000001))
                     )
    ),
    
    
    # add hover message
    bsTooltip("icon", "Specify width and height of the downloaded plot", placement = "top", trigger = "hover",
              options = NULL),
    bsTooltip("effectpct", "Input percentage of effect, from 0 to 100; Effect estimates are shown triangles in the dose-response curve plot", placement = "top", trigger = "hover",
              options = NULL),
    bsTooltip("num_width1", "Input width of error bar, from 0 to 0.1", placement = "top", trigger = "hover",
              options = NULL),
    bsTooltip("modelcompare", "for effect estimations or slope comparison if only one box was checked, or the curve comparison if both boxes were checked.", placement = "top", trigger = "hover",
              options = NULL)
    
    
    
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tags$style(type = "text/css", '#main_panel {position : absolute; width: 75%; right: 3%}'),
    id = "main_panel",
    
    # Output: Tabset w/ dose-response curve plots and interaction index curve plot ----
    tabsetPanel(type = "tabs",
                id = "tabset",
                tabPanel("Introduction", 
                         value = 1,
                         # h2("Introduction of the Robust Dose Response Estimation"),
                         h3("Objective", id = "intro"),
                         p("The Robust and Efficient Assessment of drug Potency (REAP) is 
                             developed for convenient application of the robust dose-response 
                             estimation to real-world data analysis. It presents a straightforward 
                             analytic environment for robust estimation of dose-response curve 
                             and assessment of key statistics, including implementation of 
                             statistical comparisons and delivery of customized output for 
                             graphic presentation."),
                         
                         h3("Illustrative Example"),
                         h4("Dataset Input Requirements", id = "exp"),
                         tags$ul(
                           tags$li("The input dataset should be in a csv file"),
                           tags$li("The input dataset contains three columns: Concentration, Effect and Agent"),
                           tags$li("Columns in the input dataset should follow the order of Concentration, Effect and Agent")
                         ),
                         p("It is recommended that users normalize the response variable to the range of (0,1) by 
                               themselves. Otherwise, REAP will automatically truncate the values exceeding 
                               the boundaries to (0,1) using a truncation algorithm"),
                         p("Below is an example dataset illustrating the format of the input dataset"),
                         tags$a(href='31780660_F1B_exampledata.csv', 
                                target='blank', 'Sample Data', download = '31780660_F1B_exampledata.csv'),
                         tableOutput("table1"),
                         
                         h4("Output"),
                         p("We upload the example data and obtain the following results in the output tab:"),
                         tags$img(src='shinydemo.png', height="60%", width="60%", align="center"),
                         br(),
                         p("Under the output tab, we will see a dose-response curve plot, along with 
                         tabulation for effect and model estimations. We also enable hypothesis testing 
                         for comparisons of effect estimations, slopes and models (i.e., comparing both 
                         intercepts and slopes). By default, the x-axis of the dose-response plot is 
                         log-scaled. In the plot, users can choose to add mean values and confidence 
                         intervals for data points under the same agent and dose level. Triangles indicate 
                         effect estimations. Both plots and 
                         estimation tables are downloadable on REAP to plug in presentations and 
                         manuscripts for result dissemination."),
                         
                         h4("Reference", id= "ref"),
                         p("Zhou, S*, Liu, X*, Fang, X*, Chinchilli, VM, Wang, M, Wang, HG, Dokholyan, NV, Shen, C, Lee, JJ. (2021) Robust and Efficient Assessment of Potency (REAP): A Quantitative Tool for Dose-response Curve Estimation. doi:10.1101/2021.11.20.469388")
                ),
                tabPanel("Dataset", value = 2, 
                         fluidRow(
                           column(12, span(textOutput("text_warning1"), style="color:red")),
                           column(12, span(textOutput("text_warning2"), style="color:red")),
                           br(),
                           column(12, tableOutput("table"))
                         )),
                tabPanel("Output", value = 3,
                         plotOutput("plot123"),
                         span(textOutput("ci_warning"), style="color:red"),
                         br(),
                         fluidRow(
                           column(12, htmlOutput("modelsummary")),
                           column(12, htmlOutput("summary")),
                           # column(12, h5("|m|>1: P-values on hypothesis testing on H0: m>1 or H0: m<-1.")),
                           # column(12, h5("Pairwise comparison: P-values on pairwise comparisons.")),
                           column(12, htmlOutput("modelcomparison"))
                         ),
                         br(),
                         downloadButton("drplot", "Download Plot"),
                         downloadButton("dtable", "Download Table")
                )
    )
    
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  
  # Read in the input dataset
  dt <- reactive({
    req(input$filedata)
    read.csv(input$filedata$datapath)
  })
  
  # Truncation on dataset -----
  dt.truncated <- reactive({
    truncated <- 1e-9
    
    df <- dt()
    df <- df[complete.cases(df[,2]),]
    df <- df[complete.cases(df[,1]),]
    df <- df[df[,1] > 0, ]
    df[,2][df[,2] <= truncated] = truncated
    df[,2][df[,2] >= 1-truncated] = 1-truncated
    df
  })
  
  # Output warning message
  output$text_warning1 <- renderText({
    text1 = NULL
    if (any(dt()[,2] <=0)==TRUE | any(dt()[,2]>=1)==TRUE){
      text1 = "Observations about the data are unreliable. Some of your drug effects are out of (0,1).
            For better analysis, we truncate the extreme values to be within (0,1), which may increase
            the robustness of the dose response estimation. However, it is 
            suggested that you perform the truncation or transformation yourself 
            before uploading the dataset. "
    }
    
    text1
  })
  
  output$text_warning2 <- renderText({
    text2 = NULL
    if (nrow(dt()) > nrow(dt.truncated())){
      text2 = "Your data contains missing values or some dose levels are non-positive. We delete the relavant observation to 
      ensure the following analysis."
    }
    
    text2
  })
  
  # Output dataset table
  output$table <- renderTable({
    # validate(
    #     need(dt()[,2] >0 & dt()[,2] < 1,
    #          "Observations about the data is unreliable. Please check your dataset.")
    # )
    
    if (any(dt()[,2] <=0)==TRUE | any(dt()[,2]>=1)==TRUE){
      return(dt.truncated())
    } else {
      return(dt())
    }},
    digits = 3)
  
  output$table1 <- renderTable({
    return(displaytable)
  }, digits = 3)
  
  # CI width
  observe({
    if (input$checkbox2 == TRUE){
      shinyjs::show(id = "num_width")
    } else {
      shinyjs::hide(id = "num_width")
    }
  })
  
  # CI width input warning
  output$CIinput <- renderText({
    if (input$num_width > 0.1 | input$num_width < 0){
      "Input width of the error bar should be within (0, 0.1)."
    } else {
      NULL
    }
  })
  
  
  drplots <- reactiveValues()
  
  
  # Save results of MDPDE fitting
  beta.fit <- reactive({
    nms <- colnames(dt.truncated())
    
    n = length(unique(dt.truncated()[,3]))
    beta.fit = c()
    
    for (i in 1:n){
      d.dt <- subset(dt.truncated(), dt.truncated()[,3] == as.character(unique(dt.truncated()[,3]))[i])
      nms.dt <- colnames(d.dt)
      y <- d.dt[,2] # response value
      X <- matrix(c(rep(1,nrow(d.dt)),log(d.dt[,1])),ncol=2,byrow=F) #regressor matrix for the mean submodel
      
      if (input$checkbox_betareg==TRUE){
        Z <- matrix(c(rep(1,nrow(d.dt)), log(d.dt[,1])),ncol=2,byrow=F) #regressor matrix for the precision submodel
      } else {
        Z <- matrix(c(rep(1,nrow(d.dt))),ncol=1,byrow=F) 
      }
      
      # d.betareg <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=TRUE, q0=0.5, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
      #                         linkmu="logit", linkphi="identity", weights=FALSE)
      d.betareg <- MDPDE_BETA2(y=y, X=X, Z=Z)
      
      beta.fit[[i]] = d.betareg
    }
    beta.fit
  })
  
  
  # Generate plots of dose-response curves  -----
  output$plot123 <- renderPlot({
    
    nms <- colnames(dt.truncated())
    
    n = length(unique(dt.truncated()[,3]))
    # beta.fit = c()
    # 
    # for (i in 1:n){
    #     d.dt <- subset(dt.truncated(), dt.truncated()[,3] == as.character(unique(dt.truncated()[,3]))[i])
    #     nms.dt <- colnames(d.dt)
    #     y <- d.dt[,2] # response value
    #     X <- matrix(c(rep(1,nrow(d.dt)),log(d.dt[,1])),ncol=2,byrow=F) #regressor matrix for the mean submodel
    #     
    #     if (input$checkbox_betareg==TRUE){
    #       Z <- matrix(c(rep(1,nrow(d.dt)), log(d.dt[,1])),ncol=2,byrow=F) #regressor matrix for the precision submodel
    #     } else {
    #       Z <- matrix(c(rep(1,nrow(d.dt))),ncol=1,byrow=F) 
    #     }
    # 
    #     d.betareg <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=TRUE, q0=0.5, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
    #                             linkmu="logit", linkphi="identity", weights=FALSE)
    #   
    #     beta.fit[[i]] = d.betareg
    # }
    
    dose.dt <- as.data.frame(seq(min(dt.truncated()[,1]), max(dt.truncated()[,1]), length.out=1000))
    colnames(dose.dt) <- nms[1]
    
    beta.pred <- c()
    for (i in 1:n){
      beta0 <- beta.fit()[[i]]$beta[1]
      beta1 <- beta.fit()[[i]]$beta[2]
      temp <- inv.logit(beta0 + beta1*log(dose.dt[,1]))
      beta.pred <- c(beta.pred, temp)
    }
    dt.pred <- as.data.frame(rep(dose.dt[,1],n))
    colnames(dt.pred) <- nms[1]
    dt.pred[,nms[2]] <- beta.pred
    agentnames <- as.character(unique(dt.truncated()[,3]))
    dt.pred[,nms[3]] <- rep(agentnames, each=1000) # Prediction dataset is set up
    # First 1000 data predictions are from the first agent in the dataset; second 1000 are from the second agent ...
    
    dt.ci <- summarySE(dt.truncated(), measurevar=nms[2], groupvars=c(nms[1],nms[3]))
    
    dt.ci$LCL = dt.ci[,nms[2]]-dt.ci$sd
    dt.ci$UCL = dt.ci[,nms[2]]+dt.ci$sd
    dt.ci$UCL_l = ifelse(dt.ci$UCL > 1, 1 , NA)
    dt.ci$LCL_l = ifelse(dt.ci$LCL < 0, 0 , NA)
    UCL_index = which(dt.ci$UCL > 1)
    LCL_index = which(dt.ci$LCL < 0)
    
    # Add (0,1) data points or (1,1) data points
    tmpdt <- subset(dt.truncated(), dt.truncated()[,3] == as.character(unique(dt.truncated()[,3]))[1])
    tmpdt <- tmpdt[order(tmpdt[,2]),]
    if (tmpdt[,1][1] <= tmpdt[,1][length(tmpdt[,1])]){
      df01 <- data.frame(rep(0,n),rep(0,n),unique(dt.pred[,3]))
      colnames(df01) <- nms
      dt.pred1 <- rbind(dt.pred,df01)
    } else {
      df01 <- data.frame(rep(0,n),rep(1,n),unique(dt.pred[,3]))
      colnames(df01) <- nms
      dt.pred1 <- rbind(dt.pred,df01)
    }
    
    ## Add IC50/%effect estimation points on x axis
    agentnms = unique(dt.truncated()[,3]) 
    if (is.na(input$effectpct)==FALSE){
      effect_num = 0.01*input$effectpct
      
      ic50dt = data.frame(agent = agentnms,
                          xval = NA,
                          yval = NA)
      for (i in 1:n){
        dt.pred.select = subset(dt.pred, dt.pred[,3]==agentnms[i])
        index = which.min(abs(dt.pred.select[,2]-effect_num))
        ic50dt[i,2] = dt.pred.select[index,1]
        ic50dt[i,3] = -0.03
        # print(dt.pred.select[index,])
        if (min(dt.pred.select[,2]) > effect_num | max(dt.pred.select[,2]) < effect_num){
          ic50dt[i,2] = NA
          ic50dt[i,3] = NA
        }
      }
      
      for (i in 1:nrow(ic50dt)){
        ind = which(abs(ic50dt$xval - ic50dt$xval[i])==0)
        if (length(ind) >=2){
          k=-0.03
          for (j in 2:length(ind)){
            k = k+0.015
            ic50dt$yval[ind[j]] = k
          }
        }
        
      }
    }
    
    
    ciwidth = input$num_width
    if (input$num_width > 0.1 | input$num_width < 0){ ciwidth = 0.1 }
    
    shapenum = 4
    sizenum=10* (ciwidth+0.2)
    linesize = 1
    xlabel = nms[1]
    xlabel = gsub("."," ", xlabel,fixed=TRUE)
    ylabel = nms[2]
    ylabel = gsub("."," ", ylabel,fixed=TRUE)
    ## plotci_scale -----
    if(any(is.na(dt.ci$LCL_l)==FALSE) & any(is.na(dt.ci$UCL_l)==FALSE)){ # both ucl and lcl have out-of-limit values
      plotci_scale = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum,stroke=linesize) +
        geom_line(aes(x = dt.pred[,1], y = dt.pred[,2], color=dt.pred[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + 
        # ylim(-0.05,1) +
        labs(color=nms[3], shape=nms[3]) + scale_x_continuous(trans = 'log10') + 
        geom_errorbar(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci$UCL[LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]][LCL_index]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci[,nms[2]][LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]][LCL_index]), size=linesize) +
        geom_errorbar(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci$LCL[UCL_index],color=dt.ci[,nms[3]][UCL_index]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci[,nms[2]][UCL_index],color=dt.ci[,nms[3]][UCL_index]), size=linesize) +
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$UCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)+ 
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$LCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)
      
      
    } else if (any(is.na(dt.ci$LCL_l)==FALSE)){
      plotci_scale = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum,stroke=linesize) +
        geom_line(aes(x = dt.pred[,1], y = dt.pred[,2], color=dt.pred[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + 
        labs(color=nms[3], shape=nms[3]) + scale_x_continuous(trans = 'log10') + 
        geom_errorbar(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci$UCL[LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci[,nms[2]][LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]]), size=linesize) +
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$LCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)
      
    } else if (any(is.na(dt.ci$UCL_l)==FALSE)){
      plotci_scale = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum,stroke=linesize) +
        geom_line(aes(x = dt.pred[,1], y = dt.pred[,2], color=dt.pred[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + 
        labs(color=nms[3], shape=nms[3]) + scale_x_continuous(trans = 'log10') + 
        geom_errorbar(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci$LCL[UCL_index],color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci[,nms[2]][UCL_index],color=dt.ci[,nms[3]]), size=linesize) +
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$UCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)
      
    } else { # neither of them have out-of-limit values
      plotci_scale = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum,stroke=linesize) +
        geom_line(aes(x = dt.pred[,1], y = dt.pred[,2], color=dt.pred[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + 
        labs(color=nms[3], shape=nms[3]) + scale_x_continuous(trans = 'log10')  
    }
    
    ## plotci_scale_ic50 -----
    if(any(is.na(dt.ci$LCL_l)==FALSE) & any(is.na(dt.ci$UCL_l)==FALSE)){ # both ucl and lcl have out-of-limit values
      plotci_scale_ic50 = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum,stroke=linesize) +
        geom_line(aes(x = dt.pred[,1], y = dt.pred[,2], color=dt.pred[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + 
        labs(color=nms[3], shape=nms[3]) + scale_x_continuous(trans = 'log10') + 
        geom_errorbar(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci$UCL[LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]][LCL_index]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci[,nms[2]][LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]][LCL_index]), size=linesize) +
        geom_errorbar(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci$LCL[UCL_index],color=dt.ci[,nms[3]][UCL_index]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci[,nms[2]][UCL_index],color=dt.ci[,nms[3]][UCL_index]), size=linesize) +
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$UCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)+ 
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$LCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)+
        geom_point(aes(x=ic50dt$xval, y=ic50dt$yval,color=ic50dt$agent),shape=17,size=2) 
      
      
    } else if (any(is.na(dt.ci$LCL_l)==FALSE)){
      plotci_scale_ic50 = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum,stroke=linesize) +
        geom_line(aes(x = dt.pred[,1], y = dt.pred[,2], color=dt.pred[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(nms[1]) + ylab(paste(ylabel,"(in %)")) +
        labs(color=nms[3], shape=nms[3]) + scale_x_continuous(trans = 'log10') + 
        geom_errorbar(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci$UCL[LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci[,nms[2]][LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]]), size=linesize) +
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$LCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)+
        geom_point(aes(x=ic50dt$xval, y=ic50dt$yval,color=ic50dt$agent),shape=17,size=2)
      
    } else if (any(is.na(dt.ci$UCL_l)==FALSE)){
      plotci_scale_ic50 = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum, stroke=linesize) +
        geom_line(aes(x = dt.pred[,1], y = dt.pred[,2], color=dt.pred[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) +
        labs(color=nms[3], shape=nms[3]) + scale_x_continuous(trans = 'log10') + 
        geom_errorbar(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci$LCL[UCL_index],color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci[,nms[2]][UCL_index],color=dt.ci[,nms[3]]), size=linesize) +
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$UCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)+
        geom_point(aes(x=ic50dt$xval, y=ic50dt$yval,color=ic50dt$agent),shape=17,size=2)
      
    } else { # neither of them have out-of-limit values
      plotci_scale_ic50 = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum, stroke=linesize) +
        geom_line(aes(x = dt.pred[,1], y = dt.pred[,2], color=dt.pred[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) +
        labs(color=nms[3], shape=nms[3]) + scale_x_continuous(trans = 'log10')  +
        geom_point(aes(x=ic50dt$xval, y=ic50dt$yval,color=ic50dt$agent),shape=17,size=2)
    }
    
    ## plotci_ic50 -----
    if(any(is.na(dt.ci$LCL_l)==FALSE) & any(is.na(dt.ci$UCL_l)==FALSE)){ # both ucl and lcl have out-of-limit values
      plotci_ic50 = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum, stroke=linesize) +
        geom_line(aes(x = dt.pred1[,1], y = dt.pred1[,2], color=dt.pred1[,3]),size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + 
        # ylim(-0.05,1)+
        labs(color=nms[3], shape=nms[3]) + 
        geom_errorbar(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci$UCL[LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]][LCL_index]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci[,nms[2]][LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]][LCL_index]), size=linesize) +
        geom_errorbar(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci$LCL[UCL_index],color=dt.ci[,nms[3]][UCL_index]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci[,nms[2]][UCL_index],color=dt.ci[,nms[3]][UCL_index]), size=linesize) +
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$UCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)+ 
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$LCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)+
        geom_point(aes(x=ic50dt$xval, y=ic50dt$yval,color=ic50dt$agent),shape=17,size=2)
      
      
    } else if (any(is.na(dt.ci$LCL_l)==FALSE)){
      plotci_ic50 = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum, stroke=linesize) +
        geom_line(aes(x = dt.pred1[,1], y = dt.pred1[,2], color=dt.pred1[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + 
        labs(color=nms[3], shape=nms[3]) + 
        geom_errorbar(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci$UCL[LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]][LCL_index]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci[,nms[2]][LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]][LCL_index]), size=linesize) +
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$LCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)+
        geom_point(aes(x=ic50dt$xval, y=ic50dt$yval,color=ic50dt$agent),shape=17,size=2)
      
    } else if (any(is.na(dt.ci$UCL_l)==FALSE)){
      plotci_ic50 = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum, stroke=linesize) +
        geom_line(aes(x = dt.pred1[,1], y = dt.pred1[,2], color=dt.pred1[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) +
        labs(color=nms[3], shape=nms[3]) + 
        geom_errorbar(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci$LCL[UCL_index],color=dt.ci[,nms[3]][UCL_index]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci[,nms[2]][UCL_index],color=dt.ci[,nms[3]][UCL_index]), size=linesize) +
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$UCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)+
        geom_point(aes(x=ic50dt$xval, y=ic50dt$yval,color=ic50dt$agent),shape=17,size=2)
      
    } else { # neither of them have out-of-limit values
      plotci_ic50 = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum, stroke=linesize) +
        geom_line(aes(x = dt.pred1[,1], y = dt.pred1[,2], color=dt.pred1[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + coord_cartesian(ylim=c(-0.05, 1)) +
        labs(color=nms[3], shape=nms[3]) +
        geom_point(aes(x=ic50dt$xval, y=ic50dt$yval,color=ic50dt$agent),shape=17,size=2)
      
    }
    
    
    ## plotci -----
    if(any(is.na(dt.ci$LCL_l)==FALSE) & any(is.na(dt.ci$UCL_l)==FALSE)){ # both ucl and lcl have out-of-limit values
      plotci = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum, stroke=linesize) +
        geom_line(aes(x = dt.pred1[,1], y = dt.pred1[,2], color=dt.pred1[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + 
        labs(color=nms[3], shape=nms[3]) + 
        geom_errorbar(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci$UCL[LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]][LCL_index]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci[,nms[2]][LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]][LCL_index]), size=linesize) +
        geom_errorbar(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci$LCL[UCL_index],color=dt.ci[,nms[3]][UCL_index]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci[,nms[2]][UCL_index],color=dt.ci[,nms[3]][UCL_index]), size=linesize) +
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$UCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)+ 
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$LCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)
      
    } else if (any(is.na(dt.ci$LCL_l)==FALSE)){
      plotci = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum, stroke=linesize) +
        geom_line(aes(x = dt.pred1[,1], y = dt.pred1[,2], color=dt.pred1[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + 
        labs(color=nms[3], shape=nms[3]) + 
        geom_errorbar(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci$UCL[LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]][LCL_index]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][LCL_index],ymin = dt.ci[,nms[2]][LCL_index], ymax = dt.ci$UCL[LCL_index],color=dt.ci[,nms[3]][LCL_index]), size=linesize) +
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$LCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize)
      
    } else if (any(is.na(dt.ci$UCL_l)==FALSE)){
      plotci = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum, stroke=linesize) +
        geom_line(aes(x = dt.pred1[,1], y = dt.pred1[,2], color=dt.pred1[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + 
        labs(color=nms[3], shape=nms[3]) + 
        geom_errorbar(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci$LCL[UCL_index],color=dt.ci[,nms[3]][UCL_index]), width=ciwidth, size=linesize) + 
        geom_linerange(aes(x = dt.ci[,nms[1]][UCL_index],ymin = dt.ci$LCL[UCL_index], ymax = dt.ci[,nms[2]][UCL_index],color=dt.ci[,nms[3]][UCL_index]), size=linesize) +
        geom_segment(aes(x = dt.ci[,nms[1]], xend = dt.ci[,nms[1]], y = dt.ci[,nms[2]], yend = dt.ci$UCL_l, color=dt.ci[,nms[3]]),  arrow = arrow(length = unit(ciwidth+0.2, "cm")), size=linesize) 
      
    } else { # neither of them have out-of-limit values
      plotci = ggplot() +
        geom_point(aes(x=dt.ci[,nms[1]], y=dt.ci[,nms[2]], color=dt.ci[,nms[3]]),shape=shapenum,size=sizenum, stroke=linesize) +
        geom_line(aes(x = dt.pred1[,1], y = dt.pred1[,2], color=dt.pred1[,3]), size=linesize) +
        geom_errorbar(aes(x=dt.ci[,1],ymin=dt.ci$LCL, ymax=dt.ci$UCL, 
                          color=dt.ci[,nms[3]]), width=ciwidth, size=linesize) + 
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + 
        labs(color=nms[3], shape=nms[3]) 
      
    }
    
    
    
    
    if (input$checkbox1 == TRUE){
      drplots$plottotal = ggplot() +
        geom_line(aes(x = dt.pred1[,1], y = dt.pred1[,2], color=dt.pred1[,3]), size=linesize) +
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) + 
        labs(color=nms[3], shape=nms[3]) +
        scale_color_prism("colors") + 
        scale_fill_prism("colors") + 
        theme_prism(palette = "colors", base_size = 16) +
        scale_y_continuous(
          limits = c(-0.05, 1), 
          # breaks = seq(-100, 500, 100),
          guide = "prism_offset",
          labels = function(x) x*100
        )
      
      if (input$checkbox2 == TRUE){
        drplots$plottotal = plotci +
          scale_color_prism("colors") + 
          scale_fill_prism("colors") + 
          theme_prism(palette = "colors", base_size = 16) +
          scale_y_continuous(
            limits = c(-0.05, 1), 
            # breaks = seq(-100, 500, 100),
            guide = "prism_offset",
            labels = function(x) x*100
          )
      }
      if (is.na(input$effectpct)==FALSE){
        drplots$plottotal = ggplot() +
          geom_line(aes(x = dt.pred1[,1], y = dt.pred1[,2], color=dt.pred1[,3]), size=linesize) +
          xlab(xlabel) + ylab(paste(ylabel,"(in %)")) +
          labs(color=nms[3], shape=nms[3]) +
          geom_point(aes(x=ic50dt$xval, y=ic50dt$yval,color=ic50dt$agent),shape=17,size=2)+
          scale_color_prism("colors") + 
          scale_fill_prism("colors") + 
          theme_prism(palette = "colors", base_size = 16) +
          scale_y_continuous(
            limits = c(-0.05, 1), 
            guide = "prism_offset",
            labels = function(x) x*100
          )
      }
      
    } else {
      drplots$plottotal = ggplot() +
        geom_line(aes(x = dt.pred[,1], y = dt.pred[,2], color=dt.pred[,3]), size=linesize) +
        xlab(xlabel) + ylab(paste(ylabel,"(in %)")) +
        scale_x_continuous(trans = 'log10') + 
        labs(color=nms[3], shape=nms[3]) +
        scale_color_prism("colors") + 
        scale_fill_prism("colors") + 
        theme_prism(palette = "colors", base_size = 16) +
        scale_y_continuous(
          limits = c(-0.05, 1), 
          # breaks = seq(-100, 500, 100),
          guide = "prism_offset",
          labels = function(x) x*100
        )
      
      if (input$checkbox2 == TRUE){
        drplots$plottotal = plotci_scale +
          scale_color_prism("colors") + 
          scale_fill_prism("colors") + 
          theme_prism(palette = "colors", base_size = 16) +
          scale_y_continuous(
            limits = c(-0.05, 1), 
            # breaks = seq(-100, 500, 100),
            guide = "prism_offset",
            labels = function(x) x*100
          )
      }
      
      if (is.na(input$effectpct)==FALSE){
        drplots$plottotal = ggplot() +
          geom_line(aes(x = dt.pred[,1], y = dt.pred[,2], color=dt.pred[,3]), size=linesize) +
          xlab(xlabel) + ylab(paste(ylabel,"(in %)")) +
          # xlab(expression(paste("Assay Concentration"," (",mu,"M)",sep=""))) + ylab(paste("Cell viability","(in %)")) +
          scale_x_continuous(trans = 'log10') + 
          labs(color=nms[3], shape=nms[3]) +
          geom_point(aes(x=ic50dt$xval, y=ic50dt$yval,color=ic50dt$agent),shape=17,size=2) +
          scale_color_prism("colors") + 
          scale_fill_prism("colors") + 
          theme_prism(palette = "colors", base_size = 16) +
          scale_y_continuous(
            limits = c(-0.05, 1), 
            # breaks = seq(-100, 500, 100),
            guide = "prism_offset",
            labels = function(x) x*100
          )
      }
    }
    
    if (input$checkbox2 == TRUE & is.na(input$effectpct)==FALSE){
      drplots$plottotal = plotci_scale_ic50 +
        scale_color_prism("colors") + 
        scale_fill_prism("colors") + 
        theme_prism(palette = "colors", base_size = 16) +
        scale_y_continuous(
          limits = c(-0.05, 1), 
          # breaks = seq(-100, 500, 100),
          guide = "prism_offset",
          labels = function(x) x*100
        )
    }
    if (input$checkbox1 == TRUE & input$checkbox2 == TRUE & is.na(input$effectpct)==FALSE){
      drplots$plottotal = plotci_ic50 +
        scale_color_prism("colors") + 
        scale_fill_prism("colors") + 
        theme_prism(palette = "colors", base_size = 16) +
        scale_y_continuous(
          limits = c(-0.05, 1), 
          # breaks = seq(-100, 500, 100),
          guide = "prism_offset",
          labels = function(x) x*100
        )
    }
    
    # temp <- summarySE(dt.truncated(), measurevar=nms[2], groupvars=c(nms[1],nms[3]))
    # if (input$checkbox2 == TRUE & is.na(temp) == FALSE){
    #   output$ci_warning = "Confidence interval cannot be calculated because of not enough data points."
    # } else {
    #   output$ci_warning =NULL
    # }
    
    drplots$plottotal 
  })
  
  
  model.dt <- reactive({
    n = length(unique(dt.truncated()[,3]))
    
    beta.intercept = c()
    beta.slope = c()
    # beta.fit = c()
    beta.slope.se = c()
    beta.slope.pval = c()
    # logic50.se = c()
    logic10.se = c() # can be any % effect
    for (i in 1:n){
      # d.dt <- subset(dt.truncated(), dt.truncated()[,3] == as.character(unique(dt.truncated()[,3]))[i])
      # nms.dt <- colnames(d.dt)
      # y <- d.dt[,2] # response value
      # X <- matrix(c(rep(1,nrow(d.dt)),log(d.dt[,1])),ncol=2,byrow=F) #regressor matrix for the mean submodel
      # 
      # if (input$checkbox_betareg==TRUE){
      #   Z <- matrix(c(rep(1,nrow(d.dt)), log(d.dt[,1])),ncol=2,byrow=F) #regressor matrix for the precision submodel
      # } else {
      #   Z <- matrix(c(rep(1,nrow(d.dt))),ncol=1,byrow=F) 
      # }
      # # fcn <- as.formula(paste(nms.dt[2], "~log(",nms.dt[1],")|log(",nms.dt[1],")"))
      # d.betareg <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=TRUE, q0=0.5, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
      #                         linkmu="logit", linkphi="identity", weights=FALSE)
      # beta.fit[[i]] = d.betareg
      d.betareg <- beta.fit()[[i]]
      beta.intercept = c(beta.intercept, d.betareg$beta[1])
      beta.slope = c(beta.slope, d.betareg$beta[2])
      # sum <- summary(d.betareg)
      beta.slope.se = c(beta.slope.se, d.betareg$se_beta[2])
      beta.slope.pval = c(beta.slope.pval, d.betareg$pvalue_beta[2])
      # logic50.se = c(logic50.se, deltamethod(~-x1/x2, coef(d.betareg)[1:2],vcov(d.betareg)[1:2,1:2]))
      form <- sprintf("~ (%f-x1)/x2",logit(0.01*input$effectpct))
      logic10.se = c(logic10.se, deltamethod(as.formula(form), d.betareg$beta,d.betareg[[7]][1:2,1:2]))
    }
    
    # z-test for slope on |m|>1
    beta.slope.z.pval <- c()
    for (i in 1:length(beta.slope)){
      slope = beta.slope[i]
      se = beta.slope.se[i]
      if (slope > 0){
        z = (slope-1)/se
        pval = pnorm(z, lower.tail = FALSE)
      } else {
        z = (slope+1)/se
        pval = pnorm(z, lower.tail = TRUE)
      }
      beta.slope.z.pval = c(beta.slope.z.pval, pval)
    }
    beta.slope.z.pval = ifelse(beta.slope.z.pval<0.0001,"<.0001",round(beta.slope.z.pval,4))
    
    
    
    # EC10 estimation panel -----
    ic10 = exp((logit(0.01*input$effectpct)-beta.intercept)/beta.slope)
    ic10.se = sqrt((logic10.se)^2 * ic10^2)
    ic10.z = ic10/ic10.se
    
    # Compare ic10
    nms <- colnames(dt.truncated())
    ic10.count = dt.truncated() %>%
      group_by(dt.truncated()[,3]) %>%
      dplyr::summarise(
        count_agent = n()
      )
    
    ic10.sd<-ic10.se*sqrt(as.list(ic10.count[,2])$count_agent)
    ic10.n<-as.list(ic10.count[,2])$count_agent
    # ic10.ind.anova <- ind.oneway.second(ic10,ic10.sd,ic10.n)
    ic10.ind.anova <- pairwise_compare(m=ic10,sd=ic10.sd,n=ic10.n)
    ic10.anova.pval <- pf(ic10.ind.anova$F[1], 
                          ic10.ind.anova$df[1], 
                          ic10.ind.anova$df[2], lower.tail = FALSE)
    
    
    ic10.pval = 2*pnorm(abs(ic10.z), lower.tail = FALSE)
    ic10.pval = ifelse(ic10.pval<1e-4, "<.0001", round(ic10.pval,4))
    beta.slope.pval = ifelse(beta.slope.pval<1e-4, "<.0001", round(beta.slope.pval,4))
    models <- as.character(unique(dt.truncated()[,3]))
    
    # ic50 piecewise comparison
    # first rank ic50 from low to high
    dt.ic50 <- data.frame(
      # ic50=ic50,
      # se = ic50.se,
      # sd = ic50.sd,
      # n = ic50.n,
      ic10=ic10,
      se10 = ic10.se,
      sd10 = ic10.sd,
      n10 = ic10.n,
      model = models,
      beta.slope = beta.slope,
      beta.slope.se = beta.slope.se,
      beta.slope.pval = beta.slope.pval,
      beta.slope.z.pval = beta.slope.z.pval,
      beta.intercept = beta.intercept
    )
    dt.ic50 <- dt.ic50[order(ic10),]
    ic10.pwise.pval <- c()
    
    if (nrow(dt.ic50)==1) {
      ic10.pwise.pval = "-"
    } else {
      for (i in 1:(nrow(dt.ic50)-1)){
        dt <- dt.ic50[c(i,i+1),]
        # ic10.pwise.anova <- ind.oneway.second(dt[,"ic10"],dt[,"sd10"],dt[,"n10"])
        ic10.pwise.anova <- pairwise_compare(m=dt[,"ic10"],sd=dt[,"sd10"],n=dt[,"n10"])
        pval <- pf(ic10.pwise.anova$F[1], 
                   ic10.pwise.anova$df[1], 
                   ic10.pwise.anova$df[2], lower.tail = FALSE)
        ic10.pwise.pval = c(ic10.pwise.pval, pval)
      }
      ic10.pwise.pval <- ifelse(ic10.pwise.pval<0.0001, "<0.0001", round(ic10.pwise.pval,4))
      ic10.pwise.pval = c("-", ic10.pwise.pval)
    }
  
    model <- data.frame(Model = dt.ic50$model,
                        Intercept = round(dt.ic50$beta.intercept,3),
                        Slope = round(dt.ic50$beta.slope,3),
                        Slope.Std.Err = round(dt.ic50$beta.slope.se,3),
                        Slope.Pvalue = dt.ic50$beta.slope.pval,
                        Slope.z.Pvalue = dt.ic50$beta.slope.z.pval,
                        IC10 = round(dt.ic50$ic10,3),
                        IC10.Std.Err = round(dt.ic50$se10,3),
                        IC10.Pvalue = ic10.pwise.pval
                        )
    rownames(model) = NULL
    sum.res = c()
    sum.res$model <- model
    # sum.res$ic50.anova.pval <- ifelse(ic50.anova.pval<0.0001, "<.0001",round(ic50.anova.pval,4))
    sum.res$ic10.anova.pval <- ifelse(ic10.anova.pval<0.0001, "<.0001",round(ic10.anova.pval,4))
    sum.res
  })
  
  # Compare models and slopes -----
  compare.pval <- reactive({
    nms <- colnames(dt.truncated())
    dt.dummy <- dummy_cols(dt.truncated(), select_columns = nms[3])
    
    if (length(unique(dt.truncated()[,3])) == 1){
      compare.pval = 0
    } else {
      colnames(dt.dummy)[4:ncol(dt.dummy)] <- paste("dummy",1:length(unique(dt.truncated()[,3])),sep="")
      colnms <- colnames(dt.dummy)
      
      interactpart <- dt.dummy[,4:(ncol(dt.dummy)-1)]*log(dt.dummy[,1])
      colnames(interactpart) <- paste("drug_dummy",1:(length(unique(dt.truncated()[,3]))-1),sep="")
      dt.dummy <- cbind(dt.dummy, interactpart)
      
      y <- dt.dummy[,2] # response value
      
      if (input$checkbox_betareg == TRUE){
        X0 <- matrix(c(rep(1,nrow(dt.dummy)),log(dt.dummy[,1])),ncol=2,byrow=F) #regressor matrix for the mean submodel
        Z <- matrix(c(rep(1,nrow(dt.dummy)), log(dt.dummy[,1])),ncol=2,byrow=F) #regressor matrix for the precision submodel
        # model0 <- MDPDE_BETA(y=y, X=X0, Z=Z, qoptimal=TRUE, q0=0.5, m=3, L=0.02, qmin=0.2, spac=0.02, method="BFGS", startV="CP",
        #                      linkmu="logit", linkphi="identity", weights=FALSE)
        model0 <- MDPDE_BETA2(y=y, X=X0, Z=Z)
        df0 <- ncol(X0)+ncol(Z)
        
        # Multiple intercepts and multiple slopes
        X1 <- as.matrix(cbind(rep(1,nrow(dt.dummy)),
                              dt.dummy[,4:(4+length(unique(dt.truncated()[,3]))-2)],
                              log(dt.dummy[,1]),interactpart))
        
        # model1 <- MDPDE_BETA(y=y, X=X1, Z=Z, qoptimal=TRUE, q0=0.5, m=3, L=0.02, qmin=0.2, spac=0.02, method="BFGS", startV="CP",
        #                      linkmu="logit", linkphi="identity", weights=FALSE)
        model1 <- MDPDE_BETA2(y=y, X=X1, Z=Z)
        df1 <- ncol(X1)+ncol(Z)
        
        # Multiple intercepts and one slope
        X2 <- as.matrix(cbind(rep(1,nrow(dt.dummy)),
                              dt.dummy[,4:(4+length(unique(dt.truncated()[,3]))-2)],
                              log(dt.dummy[,1])))
        # model2 <- MDPDE_BETA(y=y, X=X2, Z=Z, qoptimal=TRUE, q0=0.5, m=3, L=0.02, qmin=0.2, spac=0.02, method="BFGS", startV="CP",
        #                      linkmu="logit", linkphi="identity", weights=FALSE)
        model2 <- MDPDE_BETA2(y=y, X=X2, Z=Z)
        df2 <- ncol(X2)+ncol(Z)
        
      } else {
        X0 <- matrix(c(rep(1,nrow(dt.dummy)),log(dt.dummy[,1])),ncol=2,byrow=F) #regressor matrix for the mean submodel
        Z <- matrix(c(rep(1,nrow(dt.dummy))),ncol=1,byrow=F)  #regressor matrix for the precision submodel
        # model0 <- MDPDE_BETA(y=y, X=X0, Z=Z, qoptimal=TRUE, q0=0.5, m=3, L=0.02, qmin=0.2, spac=0.02, method="BFGS", startV="CP",
        #                      linkmu="logit", linkphi="identity", weights=FALSE)
        model0 <- MDPDE_BETA2(y=y, X=X0, Z=Z)
        df0 <- ncol(X0)+ncol(Z)
        
        # Multiple intercepts and multiple slopes
        X1 <- as.matrix(cbind(rep(1,nrow(dt.dummy)),
                              dt.dummy[,4:(4+length(unique(dt.truncated()[,3]))-2)],
                              log(dt.dummy[,1]),interactpart))
        
        # model1 <- MDPDE_BETA(y=y, X=X1, Z=Z, qoptimal=TRUE, q0=0.5, m=3, L=0.02, qmin=0.2, spac=0.02, method="BFGS", startV="CP",
        #                      linkmu="logit", linkphi="identity", weights=FALSE)
        model1 <- MDPDE_BETA2(y=y, X=X1, Z=Z)
        df1 <- ncol(X1)+ncol(Z)
        
        # Multiple intercepts and one slope
        X2 <- as.matrix(cbind(rep(1,nrow(dt.dummy)),
                              dt.dummy[,4:(4+length(unique(dt.truncated()[,3]))-2)],
                              log(dt.dummy[,1])))
        # model2 <- MDPDE_BETA(y=y, X=X2, Z=Z, qoptimal=TRUE, q0=0.5, m=3, L=0.02, qmin=0.2, spac=0.02, method="BFGS", startV="CP",
        #                      linkmu="logit", linkphi="identity", weights=FALSE)
        model2 <- MDPDE_BETA2(y=y, X=X2, Z=Z)
        df2 <- ncol(X2)+ncol(Z)
      }
      
      # betaregmodel <- betareg(value ~ dummy1+dummy2+dummy3+dummy4+dummy5+dummy6+
      #                           dummy7+log(Auranofin), data=dt.dummy)
      
      compare01 = 1 - pchisq(2*model1$log_lik - 2*model0$log_lik, df1-df0)
      compare12 = 1 - pchisq(2*model1$log_lik - 2*model2$log_lik, df1-df2)
      compare01 = ifelse(compare01<0.0001, "<.0001",round(compare01,4))
      compare12 = ifelse(compare12<0.0001, "<.0001",round(compare12,4))
      compare.pval = c(compare01,compare12)
    }
    
    compare.pval
  })
  
  # Model Formula table -----
  output$modelsummary <- renderText({
    
    n = length(model.dt()$model$Model)
    nms = colnames(dt())
    left = paste("logit(E) =", sep="")
    right = paste("* log(",nms[1],")", sep="")
    modelsum = c()
    
    for (i in 1:n){
      if (model.dt()$model$Slope[i] >= 0){
        sgn = "+"
      } else {
        sgn = "-"
      }
      temp = paste(left, model.dt()$model$Intercept[i], sgn,
                   abs(model.dt()$model$Slope[i]), right)
      modelsum = c(modelsum, temp)
    }
    
    modelsum.dt <- data.frame(Model = model.dt()$model$Model,
                              Formula = modelsum)
    
    modelsum.dt %>% 
      kable(align = "l") %>%
      kable_styling("striped") 
  })
  
  # Model comparison -----
  output$modelcomparison <- renderText({
    
    if (length(compare.pval())==1){
      comp_models <- data.frame(col1 = c("Null hypothesis", "Alternative hypothesis",
                                         "P-value"),
                                col2 = c("Fitted models same for all the agents",
                                         "At least one fitted model is different from other agents",
                                         "Only one model for the dataset. No comparison is conducted."))
      
      comp_ic50 <- data.frame(col1 = c("Null hypothesis", "Alternative hypothesis",
                                       "P-value"),
                              col2 = c("Effect estimation same for all the agents",
                                       "At least one effect estimation is different from other agents",
                                       "Only one model for the dataset. No comparison is conducted."))
      
      comp_slope <- data.frame(col1 = c("Null hypothesis", "Alternative hypothesis",
                                        "P-value"),
                               col2 = c("Slope estimation same for all the agents",
                                        "At lease one slope estimation is different from other agents",
                                        "Only one model for the dataset. No comparison is conducted."))
    } else {
      comp_models <- data.frame(col1 = c("Null hypothesis", "Alternative hypothesis",
                                         "P-value"),
                                col2 = c("Fitted models same for all the agents",
                                         "At least one fitted model is different from other agents",
                                         compare.pval()[1]))
      
      comp_ic50 <- data.frame(col1 = c("Null hypothesis", "Alternative hypothesis",
                                       "P-value"),
                              col2 = c("Effect estimation same for all the agents",
                                       "At least one effect estimation is different from other agents",
                                       model.dt()$ic10.anova.pval))
      
      comp_slope <- data.frame(col1 = c("Null hypothesis", "Alternative hypothesis",
                                        "P-value"),
                               col2 = c("Slope estimation same for all the agents",
                                        "At lease one slope estimation is different from other agents",
                                        compare.pval()[2]))
    }
    
    
    comp=NULL
    if (input$checkbox4 == TRUE){
      comp = comp_ic50
    }
    if (input$checkbox5 == TRUE){
      comp = comp_slope
    }
    if (input$checkbox4 == TRUE & input$checkbox5 == TRUE){
      comp = comp_models
    }
    if (input$checkbox4 == FALSE & input$checkbox5 == FALSE){
      comp
    } else {
      names_spaced <- c("Model Comparisons"," ")
      comp %>% 
        kable(align = "l",col.names = names_spaced) %>%
        kable_styling("striped") 
    }
    
    
  })
  
  # Estimation summary -----
  output$summary <- renderText({
    
    names_spaced <- c(
      ' ', 'Estimate (m)', 'Std.Err.', 
      'm > 1',' ',                 
      'Estimate', 'Std.Err.','Pairwise comparison')
    
    hd = paste("IC/EC",input$effectpct,"Estimation")
    
    dt.footnote <- model.dt()$model
    # names(dt.footnote)[4] <- paste0(names(dt.footnote)[4], 
    #                                 footnote_marker_symbol(1))
    # names(dt.footnote)[4] <- paste0(names(dt.footnote)[7], 
    #                                 footnote_marker_symbol(2))
    
    options(knitr.kable.NA = '')
    dt.footnote %>%
      dplyr::select(Model,Slope,Slope.Std.Err,Slope.z.Pvalue,IC10,IC10.Std.Err,IC10.Pvalue)%>%
      mutate(Slope = abs(Slope)) %>%
      tibble::add_column(new_col = NA, .after = c("Slope.z.Pvalue"))%>%
      kable(align = "c",col.names=names_spaced,format = "html") %>%
      kable_styling("striped") %>%
      column_spec(c(2,4,6,8), width = "6em") %>%
      column_spec(c(5), width = "4em") %>%
      add_header_above(c(" " = 1, "Hill Coefficient" = 3, " "=1, "Effect Estimation" = 3))%>%
      footnote(
        number = c("m > 1: p-value based on one-sided t-test for hypothesis testing on hill coefficient > 1", 
                   "Pairwise comparison: p-value based on ANOVA test (Cohen, 2000). Concentrations that give specified effect (default at 50%) by group were sorted from low to high. Hypothesis testings on equal potency (i.e., concentration for ED50/IC50 by default) were conducted pairwise with the group right above (one rank lower).",
                   "95% confidence intervals can be approximated by Estimate +/- t-value(0.975, df=n-1)*Std.Err.",
                   "Effect estimate is indicated by triangles in the dose-response curve plot.")
        # alphabet = c("Footnote A; ", "Footnote B; "),
        # symbol = c("Footnote Symbol 1; ", "Footnote Symbol 2"),
        # number_title = "Type I: ",
        # alphabet_title = "Type II: ", symbol_title = "Type III: ",
        # footnote_as_chunk = T, title_format = c("italic", "underline")
      )
    
  })
  
  
  output$dtable <- downloadHandler(
    filename = function() {"Dose-response curves results.csv"},
    content = function(file) {
      output = model.dt()$model
      output = output[,-5]
      colnames(output) <- c("Model","Intercept","Slope (m)","Std. Err for m",
                            "P-value for m>1",
                            "Effect estimation","Std. Err for effect estimation",
                            "Pairwise comparison")
      write.csv(output, file, row.names = FALSE)
    }
  )
  output$drplot <- downloadHandler(
    filename = function() {"Dose-response curves plot.pdf"},
    content = function(file) {
      ggsave(file,drplots$plottotal, width = input$graphwidth,
             height = input$graphheight)
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)

# rsconnect::deployApp("~/Downloads/Github/robust-dose-response")
