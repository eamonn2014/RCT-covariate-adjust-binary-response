#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some notes on the topic of randomisation
# https://twitter.com/ildiazm/status/1303002930723913728
# Rshiny ideas from on https://gallery.shinyapps.io/multi_regression/
# It follows that covariate imbalance, contrary # to what has been claimed by Altman, is just as much of a problem for large Studies as for  small ones
# https://twitter.com/f2harrell/status/1299755896319475712
# https://twitter.com/f2harrell/status/1298640944405807105
# but adjusted estimation does not have to be robust to be a major improvement over unadjusted analysis.  
# Using observed imbalances to find covariates to adjust for is arbitrary and reduces power by maximizing co-linearity with treatment

# https://discourse.datamethods.org/t/should-we-ignore-covariate-imbalance-and-stop-presenting-a-stratified-table-one-for-randomized-trials/547/32
# 'randomisation entitles us to ignore covariates we have not measured.'
# To me the goal of a parallel-group randomized clinical trial is to answer this question: do two patients starting out at the same point 
# (same age, severity of disease, etc.), one on treatment A and one on treatment B, end up with the same expected outcomes? This is fundamentally a completely conditional model.

# https://www.linkedin.com/pulse/stop-obsessing-balance-stephen-senn/
# https://discourse.datamethods.org/t/guidelines-for-covariate-adjustment-in-rcts/2814/2
#' Can you reconcile these two points?
#'   
#' @f2harrell' : One covariate imbalance is likely to be counterbalanced by another in opposite direction.
#' @stephensenn : One covariate imbalance likely coincides with other imbalances in same direction (thus, adjusting for one adjusts for them all)
#' Stephen John Senn #' @stephensenn
#' 22 Aug 2019
#' 1/2 So the first is true(ish) of unobserved covariates. The second  covers the fact that given observed  imbalance/balance in one covariate there may be imbalance/balance in another. So what?!
#

# https://twitter.com/f2harrell/status/1303002649080532995
# Unadjusted estimates are biased in comparison with adjusted estimates from nonlinear models, a fact well studied for decades.  Some do not call this 'bias' but the damage is the same. Literature is summarized in my ANCOVA chapter in http://hbiostat.org/doc/bbr.pdf
# https://discourse.datamethods.org/t/should-we-ignore-covariate-imbalance-and-stop-presenting-a-stratified-table-one-for-randomized-trials/547/2
# scenarios investigated in this app
# y	prognostic		      adj
# y	prognostic		      not adj
# y	unrelated		        adj
# y	unrelated		        not adj
# y	mix prog		        adj
# y	mix prog	         	not adj
# y	correlated prog		  adj
# y	correlated prog		  not adjusted
# n	correlated not prog	adj                <----no
# n	correlated not prog	not adjusted       <----no
# y	imbalances prog		  adj
# y	imbalances prog		  not adjusted
# y	imbalances not prog	adj
# y	imbalances not prog	not adjusted

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    rm(list=ls()) 
    set.seed(333) # reproducible
    library(directlabels)
    library(shiny) 
    library(shinyjs)  #refresh'
    library(shinyWidgets)
    library(shinythemes)  # more funky looking apps
    library(shinyalert)
    library(Hmisc)
    library(reshape)
    library(rms)
    library(ggplot2)
    library(tidyverse)
    library(Matrix)
    library(shinycssloaders)
    #library(googleVis)
    library(xtable)

    options(max.print=1000000)    
    
    fig.width <- 1200
    fig.height <- 500
    fig.width1 <- 1380
    fig.width8 <- 1380
    fig.height1 <- 700
    fig.width7 <- 700
    fig.height7 <- 500
    fig.width6 <- 680
    
    ## convenience functions
    p0 <- function(x) {formatC(x, format="f", digits=0)}
    p1 <- function(x) {formatC(x, format="f", digits=1)}
    p2 <- function(x) {formatC(x, format="f", digits=2)}
    p3 <- function(x) {formatC(x, format="f", digits=3)}
    p4 <- function(x) {formatC(x, format="f", digits=4)}
    p5 <- function(x) {formatC(x, format="f", digits=5)}
    
    logit <- function(p) log(1/(1/p-1))
    expit <- function(x) 1/(1/exp(x) + 1)
    inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
    is.even <- function(x){ x %% 2 == 0 } # function to identify odd maybe useful
    
    options(width=200)
    options(scipen=999)
    w=4  # line type
    ww=3 # line thickness
    wz=1 

# not used,  MSE for linear regression
    calc.mse <- function(obs, pred, rsq = FALSE){
        if(is.vector(obs)) obs <- as.matrix(obs)
        if(is.vector(pred)) pred <- as.matrix(pred)
        
        n <- nrow(obs)
        rss <- colSums((obs - pred)^2, na.rm = TRUE)
        if(rsq == FALSE) rss/n else {
            tss <- diag(var(obs, na.rm = TRUE)) * (n - 1)
            1 - rss/tss
        }
    }

    RR=.37  # used to limit correlations between variables
    sd1= 3  # for X covariates sd
    
    # links to Rdata objects uploaded to Git, these are pre run simulations see the save function
    # see line 1224 for the save function
    pp<- "https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/A%205000%20default%20settings%20theta%20log1.5%20-1.00%20-0.67%20-0.43.Rdata" # 5000 default log1.5 -1 -.67 -.43
    pp2<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/B%205000%20default%20settings%20theta%20log0.5%20-1.68%20-1.39%20%200.71.Rdata"
    pp3<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/C%205000%20default%20settings%20theta%20log2%20-3.46%20-1.05%20%201.15%20p1=.75.Rdata"
    pp4<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/D_10Ksims_5covariates_p1_0.12_theta_log1.3_covariates_-1.02%20_0.42_0.43_0.61%20_1.01_3_prog.Rdata"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ui <-  fluidPage(theme = shinytheme("journal"), #https://www.rdocumentation.org/packages/shinythemes/versions/1.1.2 , paper another option to try
                 # paper
                 useShinyalert(),  # Set up shinyalert
                 setBackgroundColor(
                     color = c( "#2171B5", "#F7FBFF"), 
                     gradient = "linear",
                     direction = "bottom"
                 ),
                 
                 
                 h2("Covariate adjustment in randomised controlled trials (RCTs) with a binary response"), 
                 
                 h4("The main value of randomization in RCTs is that it eliminates selection bias, treatment groups are on average comparable in terms of known
                and unknown patient characteristics [1]. But as Stephen Senn has stated '1) randomised controlled trials don't deliver balance *even* if they are very large and 
                2) valid inference does *not* depend on having balanced groups', facts that do not seem 
                to be common knowledge [2]. As Senn says elsewhere, 'Balance is valuable as a contribution to efficiency. It has nothing to do with validity' [3]. 
                This app looks into these points and investigates a related common misconception concerning RCTs; 
                it is mistakenly thought there is no need to include baseline covariates in the analysis.
                
                Many RCTs are analysed in a simple manner using only the randomised treatment as the independent variable. 
                When the response outcome is continuous, 
                precision of the treatment effect estimate is improved when adjusting for baseline covariates. 
                
                Due to randomisation we do not expect covariates to be related to the treatment assignment, 
                but they may be related to the outcome and so are not considered confounding. 
                Differences between the outcome which can be attributed to differences in the covariates can be removed, 
                resulting in a more precise estimate of treatment effect. This should be considered more often as sample sizes can be reduced.
                
                As Frank Harrell has said, 'unadjusted analysis makes the most severe assumptions of all (that risk factors do not exist)' [4].

                Note: The total effect of covariates has to be bounded. For example the range of human fasting blood glucose levels is approx. 70 to 130 mg/dL 
                and if we were simulating this response and adding similar covariates into a model this will result in a response the variance of which keeps on increasing 
                and soon implausible fasting blood glucose levels will result. 
                In fact a single continuous covariate could be used as a linear predictor or risk score that summarizes the multivariable contribution of a set of predictor variables [5,6]."), 
                 
                 h4("With this app Monte Carlo simulation is used to generate an RCT with patients randomised 1:1 to treatment and control with a binary response, estimating treatment effects whilst examining adjusting and not adjusting for covariates related to the outcome, 
                covariates not related to the outcome, collinear or correlated covariates related to the outcome and imbalanced covariates both of prognostic value and unrelated to the outcome.
                As the variance of the response increases with more covariates in the simulation, it is advisable to limit the number of covariates, the default is 3. 
                As the number of simulations to get smooth curves is high, the application may time out before simulations complete. Therefore take the code and run on your own machine. 
                There is also a tab that will allow the user to load in pre run simulations that can be examined graphically using the radio buttons.
                Note, the prognostic strength of treatment may be small compared with patient characteristics,
                such as age as in the GUSTO-1 trial. This phenomenon is observed in many prognostic evaluations of RCTs:
                treatment has a 'statistically significant' impact on outcome, but its relevance is small
                compared to other prognostic factors [7,8].
                The limited simulations I have done support adjusting over not adjusting, the AIC being smaller when adjusting. 
                Note we have ideal data conforming to statistical distributions, no missing data and all covariates are truly linear and continuous.
                "), 
                 
                h3("  "), 
                 sidebarLayout(
                     
                     sidebarPanel( width=3 ,
                                   
                                   tags$style(type="text/css", ".span8 .well { background-color: #00FFFF; }"),
                                   
                                   actionButton(inputId='ab1', label="R Shiny ",   icon = icon("th"),   
                                                onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/RCT-covariate-adjust-binary-response/master/cov-adj-binary-response/app.R', '_blank')"), 
                                   actionButton(inputId='ab1', label="R code",   icon = icon("th"),   
                                                onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/RCT-covariate-adjust-binary-response/master/cov-adj-binary-response/R%20code.R', '_blank')"),  
                                   actionButton("resample", "Simulate a new sample"),
                                   br(),  
                                   tags$style(".well {background-color:#b6aebd ;}"), 
                                   
                                   # h4("User inputs"),
                                   div(
                                       # font colours for button font
                                       tags$head(
                                           tags$style(HTML('#upload{color:black}'))    
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#upload2{color:black}'))    
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#upload3{color:black}'))    
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#upload4{color:black}'))    
                                       ),
                                       
                                       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
                                       
                                       # colours for button background
                                     
                                       tags$head(
                                           tags$style(HTML('#upload{background-color:chartreuse}'))
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#upload2{background-color:chartreuse}'))
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#upload3{background-color:chartreuse}'))
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#upload4{background-color:chartreuse}'))
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#ab1{background-color:orange}'))
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#resample{background-color:orange}'))
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#resample2{background-color:orange}'))
                                       ),
                                       
                      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
                                       splitLayout(
                                           
                                           textInput('K', 
                                                     div(h5(tags$span(style="color:blue", "No of covariates (min 2)"))), "3"),
                                           
                                           textInput('Kp', 
                                                     div(h5(tags$span(style="color:blue", "Make covariates 1:n prognostic"))), "2")
                                       ),
                      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
                      
                                       tags$hr(),
                                       splitLayout(
                                           textInput('p1', #
                                                     div(h5(tags$span(style="color:blue", "Prob in placebo"))), ".35"),  
                                           
                                           textInput('theta', 
                                                     div(h5(tags$span(style="color:blue", "Treatment effect log odds ratio"))), "log(1.5)")  #log(1.5)
                                           
                                       ),
                      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
                      
                                       splitLayout(
                                           
                                           textInput('pow', 
                                                     div(h5(tags$span(style="color:blue", "Power (%)"))), "90"),
                                           textInput('alpha', 
                                                     div(h5(tags$span(style="color:blue", "Alpha level two sided (%)"))), "5")
                                       ),
                                       tags$hr(),
                                       textInput('Fact', 
                                                 div(h5(tags$span(style="color:blue", "Covariate coefficients. Here a multiplicative factor is selected so that betas are randomly chosen between (-X*treatment effect) and (X*treatment effect)"))), "3"),
                                       
                                       
                                       textInput('simuls', 
                                                 div(h5(tags$span(style="color:blue", "Number of simulations (simulation tab only)"))), "50"),
                                       tags$hr(), 
                                       
                                       textInput('covar', 
                                                 div(h5(tags$span(style="color:blue", "Covariate distribution 1: uniform(-1,1), 2: normal(0,1)"))), "2"),
                                       
                                       ###https://stackoverflow.com/questions/49616376/r-shiny-radiobuttons-how-to-change-the-colors-of-some-of-the-choices
                                       
                                       radioButtons(
                                           inputId = "dist",
                                           label =  div(h5(tags$span(style="color:blue","Plot choices for all figures, select to present :"))),
                                           choiceNames = list(
                                               HTML("<font color='blue'>All scenarios</font>"), 
                                               tags$span(style = "color:blue", "Covariates all of prognostic value only"), 
                                               tags$span(style = "color:blue", "Covariates all of no prognostic value only"), 
                                               tags$span(style = "color:blue", "Mix of prognostic and non prognostic covariates only"),
                                               tags$span(style = "color:blue", "Correlated covariates all of prognostic value only"),
                                               tags$span(style = "color:blue", "Imbalanced covariates all of prognostic value only"),
                                               tags$span(style = "color:blue", "Imbalanced covariates of no prognostic value only")
                                               
                                           ),
                                           choiceValues = c("All", "d1", "d3", "d5",  "d7", "d9", "d11")
                                       )
                                       
                                   ),
                                   
                     ),
                     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tab panels
                     mainPanel(width=9, #eight=4,
                               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               navbarPage(       
                                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
                                   tags$style(HTML("
                            .navbar-default .navbar-brand {color: orange;}
                            .navbar-default .navbar-brand:hover {color: blue;}
                            .navbar { background-color: #b6aebd;}
                            .navbar-default .navbar-nav > li > a {color:black;}
                            .navbar-default .navbar-nav > .active > a,
                            .navbar-default .navbar-nav > .active > a:focus,
                            .navbar-default .navbar-nav > .active > a:hover {color: pink;background-color: purple;}
                            .navbar-default .navbar-nav > li > a:hover {color: black;background-color:yellow;text-decoration:underline;}
                            .navbar-default .navbar-nav > li > a[data-value='t1'] {color: red;background-color: pink;}
                            .navbar-default .navbar-nav > li > a[data-value='t2'] {color: blue;background-color: lightblue;}
                            .navbar-default .navbar-nav > li > a[data-value='t3'] {color: green;background-color: lightgreen;}
                   ")),
                                   
                                   tabPanel( "1 Simulation - to adjust or not to adjust",
                                             
                                             h4(htmlOutput("textWithNumber1a") ),
                                             fluidRow(
                                                 column(width = 6, offset = 0, style='padding:1px;',
                                                        shinycssloaders::withSpinner(
                                                            div(plotOutput("reg.plotx",  width=fig.width8, height=fig.height7)),
                                                        ),
                                                        div(plotOutput("reg.ploty",  width=fig.width8, height=fig.height7)),
                                                 ) ,
                                                 
                                                 
                                                 fluidRow(
                                                     column(width = 6, offset = 0, style='padding:1px;',
                                                            
                                                     ))),#
                                             
                                             h4(paste("Here we perform simulations investigating the treatment effect estimate and associated standard error when there are
                                           covariates that are prognostic, covariates unrelated to the outcome, a mix of prognostic and covariates unrelated to the outcome, correlated covariates
                                           and imbalance prognostic covariates and imbalanced covariates of no prognostic value. For each scenario we adjust and also do not adjust for the covariates.
                                           The default number of simulations
                                           is set at a lowly 50 so that results appear quickly. It is advisable to increase this number.
                                                    The top panel shows the distribution of the treatment effect estimates, the lower panel the associated standard error estimates. The true value is shown
                                                    by the grey vertical lines. The same covariates are used for investigations of covariates with prognostic value,
                                                    covariates unrelated to the outcome, a mix of prognostic and covariates unrelated to the outcome.
                                                    For imbalanced and correlated investigations covariates will be unique. Correlations are capped at +/- 0.37.
                                                    In the case of the imbalanced scenario, an imbalance is induced for all covariates by way of the treatment arm being derived from a Normal(0.3, 1) 
                                                    and the control arm
                                                    from a Normal(0, 1) distribution. We can also investigate scenarios of imbalanced covariates derived from a uniform distribution Uniform(-1,1) in control
                                                    and Uniform(-0.8,1.2) in the treatment arm.

                                                 ")),
                                             
                                             
                                             
                                             
                                             
                                             h4(paste("Table 2 Summary, sorted by smallest AIC estimate")),
                                             
                                             div( verbatimTextOutput("zz") )  ,
                                             # h4(htmlOutput("textWithNumber99",) ),
                                             # div( verbatimTextOutput("mse.target") )  ,
                                             h4(paste("Here are the true log odds ratio coefficients of the covariates used in the simulation: ")),
                                             div( verbatimTextOutput("betas") )  ,
                                             width = 30 )     ,
                                   
                                   # this was a n intital thought to allow user to load Rdata by browsing, It works but not very practical
                                   # better to supply links on web which is the next tab
                                   
                                   # tabPanel( "99 Load",
                                   #           
                                   #           fileInput(inputId="file1", "Upload a pre-run simulation",
                                   #                     multiple = FALSE,
                                   #                     accept = c(".Rdata" )),
                                   #           
                                   #           div(plotOutput("reg.plotL",  width=fig.width8, height=fig.height7)),
                                   #           div(plotOutput("reg.plotM",  width=fig.width8, height=fig.height7)),
                                   #           tableOutput('contents3'),
                                   # ) ,
                                   
                                   # here have action buttons clicking will load a pre run simulation that can be examined
                                   tabPanel( "2 Load in pre-run simulations",
                                             
                                             # fancier action buttons shinywidgets
                                             # actionBttn(
                                             #     inputId = "upload",
                                             #     label = "Hit to load, default settings, the covariate coeficients used were -1, -.67, -.43",
                                             #     color = "royal",
                                             #     style = "float",
                                             #     icon = icon("sliders"),
                                             #     block = TRUE,
                                             #     no_outline=TRUE
                                             # ),
                                          #   h4(""),
                                             # actionBttn(
                                             #     inputId = "upload2",
                                             #     label = "Hit to load, default settings except that treatment effect is log(0.5). The covariate coeficients used were -1.68, -1.39, 0.71.",
                                             #     color = "royal",
                                             #     style = "float",
                                             #     icon = icon("sliders"),
                                             #     block = TRUE,
                                             #     no_outline=FALSE
                                             # ),
                                             #  tags$hr(),
                                             #br(),
                                            # h4(""),
                                             # actionBttn(
                                             #     inputId = "upload3",
                                             #     label = "Hit to load, default settings except that treatment effect is log(2), intercept probability 0.7. The covariate coeficients used were -3.46, -1.05, 1.15",
                                             #     color = "royal",
                                             #     style = "float",
                                             #     icon = icon("sliders"),
                                             #     block = TRUE
                                             # ),
                                             
                                          #   h4(""),
 
                                             # actionBttn(
                                             #     inputId = "upload4",
                                             #     label = "Hit to load, 10K simulations, 5 covariates (3 prognostic), treatment effect is log(1.3),  intercept probability 0.12 The covariate coeficients used were -1.02  0.42  0.43  0.61  1.01",  
                                             #     color = "royal",
                                             #     style = "float",
                                             #     icon = icon("sliders"),
                                             #     block = TRUE
                                             # ),
 
                                             # actionBttn(
                                             #     inputId = "upload4",
                                             #     label = "Hit to load, 10K simulations, 5 covariates (3 prognostic), treatment effect is log(1.3),  intercept probability 0.12 The covariate coeficients used were -1.02  0.42  0.43  0.61  1.01",  
                                             #     color = "royal",
                                             #     style = "float",
                                             #     icon = icon("sliders"),
                                             #     block = TRUE
                                             # ),
 
                                             # 
                                             # label1 <- c("Hit to load in a pre-run simulation, 5000 sims, default settings except that treatment effect is log(1.5). The covariate coefficients are -1, -.67, -.43"),
                                             # label2 <- c("Hit to load in a pre-run simulation, 5000 sims, default settings except that treatment effect is log(1.5). The covariate coefficients are -1, -.67, -.43"),
                                             # label3 <- c("Hit to load in a pre-run simulation, 5000 sims, default settings except that treatment effect is log(1.5). The covariate coefficients are -1, -.67, -.43"),
                                             # 
                                             # label1 <-c("xxxx"),
                                             # label2<-c("xxxx"),
                                             # label3<-c("xxxx"),  
 
                                             # label1 <-c("xxxx"),
                                             # label2<-c("xxxx"),
                                             # label3<-c("xxxx"),
 
                                             
                                             
                                             
                                             #     actionButton(inputId='upload', label=label1,    icon = icon("th") ),
                                             #    actionButton(inputId='upload2', label=label2,   icon = icon("th")),
                                             #   actionButton("upload3", label3,                 icon = icon("th")),
                                             
                                             
                                             
                                             fluidRow(
                                                 
                                                 column(1,
                                                        actionBttn(
                                                            inputId = "upload",
                                                            label = "",
                                                            color = "royal",
                                                            style = "float",
                                                            icon = icon("sliders"),
                                                            block = TRUE,
                                                            no_outline=TRUE
                                                        ),
                                                        
                                                 ),
                                                 
                                                 h4("Hit to load, default settings, the covariate coefficients used were -1 -0.67, -0.43"),
                                             ),
                                             
                                             
                                             fluidRow(
                                                 
                                                 column(1,
                                                        actionBttn(
                                                            inputId = "upload2",
                                                            label = "",
                                                            color = "royal",
                                                            style = "float",
                                                            icon = icon("sliders"),
                                                            block = TRUE,
                                                            no_outline=FALSE
                                                        ),
                                                        
                                                 ),
                                                 
                                                 h4("Hit to load, default settings except that treatment effect is log(0.5). The covariate coefficients used were -1.68 -1.39 0.71."),
                                             ),
                                             
                                             fluidRow(
                                                 
                                                 column(1,
                                                        actionBttn(
                                                            inputId = "upload3",
                                                            label = "",
                                                            color = "royal",
                                                            style = "float",
                                                            icon = icon("sliders"),
                                                            block = TRUE
                                                        ), 
                                                        
                                                 ),
                                                 
                                                 h4("Hit to load, default settings except that treatment effect is log(2), intercept probability 0.7. The covariate coefficients used were -3.46 -1.05 1.15"),
                                             ),
                                             
                                             
 
                                             fluidRow(
                                                 
                                                 column(1,
                                                      
                                                        
                                                        actionBttn(
                                                            inputId = "upload4",
                                                            label = "",  
                                                            color = "royal",
                                                            style = "float",
                                                            icon = icon("sliders"),
                                                            block = TRUE
                                                        ),
                                                        
                                                            
                                                        
                                                        
                                                 ),
                                                 h4("Hit to load, 10K simulations, 5 covariates (3 prognostic), treatment effect is log(1.3), intercept prob. 0.12. The covariate coefficients used were -1.02  0.42  0.43  0.61  1.01"),
                                                 
                                             ),
 
                                          shinycssloaders::withSpinner(plotOutput("plot1",  width=fig.width8, height=fig.height7),5),
                                          
                                          shinycssloaders::withSpinner(plotOutput("plot2",  width=fig.width8, height=fig.height7),5),
                                          
                                          shinycssloaders::withSpinner(verbatimTextOutput("content1"),type = 5),
                                          
                                          
                                          
                              
                                             # # this spinner indicating something is loading does not seem to work
                                             # shinycssloaders::withSpinner(
                                             #     div(plotOutput("plot1",  width=fig.width8, height=fig.height7)),  #trt est plot
                                             # ) ,
                                             # # this spinner indicating something is loading does not seem to work
                                             # shinycssloaders::withSpinner(
                                             #     div(plotOutput("plot2",  width=fig.width8, height=fig.height7)),  #se est plot
                                             # ) ,
                                             # # this spinner indicating something is loading does work
                                             # shinycssloaders::withSpinner( 
                                             #     verbatimTextOutput('content1'),  #summary table
                                             #     
                                             # ),
                                   ),            
                                   
                                   tabPanel("3 Understanding bivariate Log Regresion", value=3, 
                                            
                                            actionButton("resample2", "Simulate a new sample"),
                                            
                                            h4("The relationship between 2by2 table and logistic regression is explained") ,
                                            
                                            splitLayout(
                                                
                                                textInput('NN', 
                                                          div(h5(tags$span(style="color:blue", "Sample size"))), "300"),
                                                
                                                textInput('pp1', 
                                                          div(h5(tags$span(style="color:blue", "Expected baseline proportion"))), ".25"),
                                                
                                                textInput('pp2', 
                                                          div(h5(tags$span(style="color:blue", "Proportion we are shooting for"))), ".35"),
                                                
                                                textInput('or', 
                                                          div(h5(tags$span(style="color:blue", "Odds Ratio we are shooting for"))), "1.5"),
                                      
                                                textInput('allocation', 
                                                          div(h5(tags$span(style="color:blue", "Randomisation allocation"))), "0.5")
                                         
                                            ),
                                            
                                          actionButton("sim","simulate"),
                                            h4("Power via simulation"),    
                                            withSpinner(verbatimTextOutput("pow1")),
                                            h4("Power via Frank Harrell Hmisc function"),    
                                            withSpinner(verbatimTextOutput("pow2")),
                                            h4("simulate one data set, columns are treatment groups, rows observed response"),  
                                            tableOutput("obs"),
                                            
                                           # verbatimTextOutput("pow") %>% withSpinner(color="#0dc5c1"))
                                          
                                            h4(""),                           
                                   
                                   ),
                                   
                                   
                                   tabPanel("4 Notes & references", value=3, 
                                            
                                            h4("First, a power calculation function in R for proportions, using the intercept proportion, treatment effect, alpha and power is used to determine the sample size.") ,
                                            
                                            h4("Tab 1, presents the results of simulation where we investigate (i) adjusting for true prognostic covariates (i.a) ignoring them in the analysis. We investigate (ii)
                                       adjusting for covariates unrelated to the outcome (ii.a) ignoring them in the analysis. We investigate (iii) adjusting for covariates some unrelated and some related to the outcome (iii.a) 
                                       ignoring them in the analysis. We investigate (iv) adjusting for covariates related to the outcome which are correlated with each other (iv.a) ignoring them in the analysis (v) 
                                       adjusting for covariates of prognostic value that are imbalanced (v.a) ignoring them in the analysis and finally (vi) adjusting for imbalanced covariates of no prognostic value 
                                       (vi.a) ignoring them in the analysis. Plots of treatment effect estimates and standard error estimates are presented as well as a summary of simulations from which we can draw conclusions.") ,
                                            
                                            
                                            
                                            h4("
                                         The next input determines the intercept (probability of response in the baseline/placebo group). The treatment effect a log odds ratio is entered determined in the next input box. Power and the alpha level two sided can be selected using the next two boxes. A power calculation for a proportion is executed resulting in a sample size used in the app. The next input box is used to determine the number of simulations used in tab 1. The last input box determines the 
range over which the covariate beta coefficients are randomly selected, using +/- multiples of the true treatment effect. The bottom of the user input section allows control over which simulated scenario results are presented graphically. We can also simulate a new sample or check the code behind the app by hitting the orange buttons.
                                          We can also simulate a new sample or check the code behind the app by hitting the orange buttons.
                                          "),
                                            h4("
                                          The next tab present an opportunity to load in pre-run simulations that can be examined graphically.
                                          "),
                                            h4("
                                          The standard error of the treatment for displaying on the plot is calculated using the intercept probability, odds ratio and sample size assuming a 50:50 randomisation. 
                                          "),
                                            
                                            h4("
                                               The following statements are taken from from Frank Harrell and  James C Slaughter's 
                                               'Biostatistics for
Biomedical Research chapter 13': 'Use of binary logistic model for covariable adjustment will result in an increase in the S.E. of the treatment effect (log odds ratio)'. 
'But even with perfect balance, adjusted OR not equal to unadjusted OR'. 'Adjusted OR will be greater than unadjusted OR'. 'It is always more efficient to adjust for predictive covariates when logistic models are used, and thus in this regard the behavior of logistic regression is the same as that of classic linear regression.'
'...for logistic, Cox, and paired survival models unadjusted treatment effects are asymptotically biased low in absolute value.'
                                               
                                               "),
                                            
                                            
                                            
                                            column(width = 12, offset = 0, style='padding:1px;',
                                                   
                                                   tags$hr(),
                                                   div(h4("References:")),  
                                                   tags$a(href = "https://www.bmj.com/content/bmj/340/bmj.c869.full.pdf", tags$span(style="color:blue", "[1] CONSORT 2010 Explanation and Elaboration: updated guidelines for reporting parallel group randomised trials"),),   
                                                   div(p(" ")),
                                                   tags$a(href = "https://www.linkedin.com/pulse/stop-obsessing-balance-stephen-senn/", tags$span(style="color:blue", "[2] Stephen Senn, Stop obsessing about balance"),),   
                                                   div(p(" ")),
                                                   tags$a(href = "https://discourse.datamethods.org/t/should-we-ignore-covariate-imbalance-and-stop-presenting-a-stratified-table-one-for-randomized-trials/547/32", tags$span(style="color:blue", "[3] Stephen Senn, point 4, Should we ignore covariate imbalance and stop presenting a stratified table one for randomized trials"),),  
                                                   div(p(" ")),
                                                   tags$a(href = "https://twitter.com/f2harrell/status/1298640944405807105",  tags$span(style="color:blue", "[4]  Frank Harrell, twitter, 'unadjusted analysis makes the most severe assumptions of all (that risk factors do not exist)'."),),   
                                                   div(p(" ")),
                                                   tags$a(href = "https://statistics.fas.harvard.edu/files/statistics/files/21_stephen_senn.pdf", tags$span(style="color:blue", "[5] Randomisation isn’t perfect but doing better is harder than you think "),),   
                                                   div(p(" ")),
                                                   tags$a(href = "https://onlinelibrary.wiley.com/doi/epdf/10.1002/sim.8570", tags$span(style="color:blue", "[6] Graphical calibration curves and the integrated calibration index (ICI) for survival models, Statistics in Medicine. 2020;1–29 "),),  
                                                   div(p(" ")),
                                                   tags$a(href = "https://www.sciencedirect.com/science/article/abs/pii/S0002870300900012?via%3Dihub", tags$span(style="color:blue", "[7] Steyerberg, E. W., Bossuyt, P. M. M., & Lee, K. L. (2000). Clinical trials in acute myocardial infarction: Should we adjust for baseline characteristics? American Heart Journal, 139(5), 745–751. doi:10.1016/s0002-8703(00)90001-2"),),   
                                                   div(p(" ")),
                                                   tags$a(href = "http://clinicalpredictionmodels.org/", tags$span(style="color:blue", "[8] Steyerberg, E. W., Clinical Prediction Models, 2019 p459"),),   
                                                   div(p(" ")),
                                                   tags$a(href = "https://twitter.com/f2harrell/status/1299755896319475712", tags$span(style="color:blue", "[9] Frank Harrell, twitter, Adjusted analysis"),),   
                                                   div(p(" ")),
                                                   tags$a(href = "https://discourse.datamethods.org/t/guidelines-for-covariate-adjustment-in-rcts/2814/2", tags$span(style="color:blue", "[10 Frank Harrell, Guidelines for covariate adjustment in rcts"),),  
                                                   div(p(" ")),
                                                   tags$a(href = "https://www.fharrell.com/post/covadj/", tags$span(style="color:blue", "[11] E.Steyerberg explains some of the advantages of conditioning on covariates"),),  
                                                   div(p(" ")),
                                                   tags$a(href = "https://hbiostat.org/doc/bbr.pdf", tags$span(style="color:blue", "[12] Biostatistics for Biomedical Research Frank E Harrell Jr. James C Slaughter Updated August 5, 2020"),),  
                                                   div(p(" ")),
                                                   
                                                   tags$hr()
                                            ) 
                                            
                                   )
                                   
                                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   END NEW   
                               )
                               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                     )
                 ) 
                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end tab panels 
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

server <- shinyServer(function(input, output   ) {
    
    #############################
    
    shinyalert("Welcome! \nAdjusting for covariates in binary response RCT!",
               "Best to do it!", 
               type = "info")
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This is where a new sample is instigated and inputs converted to numeric
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    random.sample <- reactive({
        
        foo <- input$resample
        
        K <- as.numeric(unlist(strsplit(input$K,",")))
        
        Kp <- as.numeric(unlist(strsplit(input$Kp,",")))
        
        pow <- as.numeric(unlist(strsplit(input$pow,",")))
        
        alpha <- as.numeric(unlist(strsplit(input$alpha,",")))
        
        p1 <- as.numeric(unlist(strsplit(input$p1,",")))
        
        theta <- as.numeric(eval(parse(text= (input$theta))))  # note the code here
        
        simuls <- (as.numeric(unlist(strsplit(input$simuls,","))))    
        
        covar <- (as.numeric(unlist(strsplit(input$covar,","))))   
        
        Fact <- (as.numeric(unlist(strsplit(input$Fact,","))))
        
        return(list(  
            K=K,  
            Kp=Kp,  
            pow=pow/100,
            p1=p1,
            alpha=alpha/100, 
            theta=theta,
            simuls=simuls,
            covar=covar,
            Fact=Fact
        ))
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # tab 1 simulate data (covariates and response)  
    # create  response with prognostic covariate
    # create covariates that are not prognostic
    # create a mix of above 2
    # alos look at the difference of the covariates across arms
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mcmc <- reactive({
        
        sample <- random.sample()
        
        K=sample$K
        Kp=sample$Kp
        pow=sample$pow
        p1=sample$p1
        theta=sample$theta
        alpha=sample$alpha
        covar=sample$covar
        Fact=sample$Fact
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #power the examples
        post.odds <- p1/(1-p1) * exp(theta)
        post.prob <- post.odds/(1+post.odds)
        Po <- bsamsize(p1, post.prob, fraction=.5, alpha=alpha, power=pow)  # getting the sampe size
        
        MM <- N <-ceiling(Po[1][[1]])*2  # total will always be even
        bigN <- MM  
        N1 <- MM/2
        N2 <- N1
        
        # not used, as I want se for log odds ratio 
        se. <- sqrt(  p1*(1-p1) /(N1) +    post.prob*(1- post.prob)/(N2) ) 
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # For the se of diff, n1=n2 , a simulation may have different numbers allocated to treatments 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        return(list(  se.=se.,
                      Na=N1,
                      Nb=N2,
                      N=N, 
                      bigN=bigN
        ))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # SIMULATION CODE STARTS HERE
    # here is code to simulate scenarios prognostic covariates, covariates unrelated to y, mix of pro and unrelated to y covariates, 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    simul <- reactive({
        
        sample <- random.sample()
        # need to rename to avoid recursive issues
        K1=sample$K
        Kp=sample$Kp
        pow=sample$pow
        theta1=sample$theta        
        alpha=sample$alpha  
        covar=sample$covar
        Fact=sample$Fact
        p1=sample$p1
        simuls=sample$simuls
        
        N1 <- mcmc()$N # 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        # making up some beta coefficients, fixed for all simulations as it is outside loop
        tx <- abs(-theta1)
        b1 <- round(sort(runif(K1, -tx*Fact, tx*Fact)), digits=2) 
        
        simfun <- function(N=N1, K=K1, a=log(p1/(1-p1)), sigma=sigma1, theta=(theta1), b=b1) {
            
            # randomi <- runif(N)  
            # we can select this, does not seem to have a big impact
            if (covar==1) {  
                X <- array(runif(N*K , -1,1), c(N,K))     # initially covars were uniform dist
            } else {
                X <- array(rnorm(N*K, 0, sd1), c(N,K))   
            }
            
            z <- sample(c(0,1), N, replace=T)            # treatment indicators
            
            # NEW CODE
            lp<-NULL
            lp = a+ X %*% b + theta*z 
            y <- ifelse(runif(N)  < plogis(lp), 1, 0)    # one liner RANDOM!!!
            lp<-NULL
            lp = a + theta*z 
            y2 <- ifelse(runif(N)  < plogis(lp), 1, 0)   # one liner RANDOM!!!
            lp<-NULL
            lp <- a+ X[,1:Kp] %*% b[1:Kp] + theta*z      # just have 3 variables associated with response
            y3 <- ifelse(runif(N)  < plogis(lp), 1, 0)   # one liner RANDOM!!!
            
            data.frame(X=X, y=y, z=z, y2=y2, y3=y3)
            
        }
        
        #https://stackoverflow.com/questions/5251507/how-to-succinctly-write-a-formula-with-many-variables-from-a-data-frame
        
        # function to run the analyses
        statfun <- function(d) {
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            zz <- glm(y~.-y2-y3, data=d,   family = "binomial")     ## adjusting for prognostic X, y2 is not included by use of the '-'
            f <-  summary(zz)
            
            zz1 <- glm(y~z, data=d,   family = "binomial")          ## not adjusting for prognostic X, only trt. indicator included
            f1 <-  summary(zz1)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            zz2 <- glm(y2~.-y-y3,data=d,   family = "binomial")     ## adjusting for X which are not prognostic
            f2 <-  summary(zz2)
            
            zz3 <- glm(y2~z, data=d,   family = "binomial")         ## not adjusting for X which are not prognostic
            f3 <-  summary(zz3)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            zz4 <- glm(y3~.-y-y2,data=d,   family = "binomial")     ## adjusting some X  are prognostic
            f4 <-  summary(zz4)
            
            zz5 <- glm(y3~z, data=d,   family = "binomial")         ## not adjusting when some X  are prognostic
            f5 <-  summary(zz5)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
            
            cbind(
                
                coef(f)["z", "Estimate"],
                coef(f)["z", "Std. Error"],
                
                coef(f1)["z", "Estimate"],
                coef(f1)["z", "Std. Error"],
                
                coef(f2)["z", "Estimate"],
                coef(f2)["z", "Std. Error"],
                
                coef(f3)["z", "Estimate"],
                coef(f3)["z", "Std. Error"], #8
                
                coef(f4)["z", "Estimate"],
                coef(f4)["z", "Std. Error"],
                
                coef(f5)["z", "Estimate"],
                coef(f5)["z", "Std. Error"],
                
                # collect p values
                coef(f)["z", "Pr(>|z|)"]  < alpha,  #13    
                coef(f1)["z", "Pr(>|z|)"] < alpha,
                coef(f2)["z", "Pr(>|z|)"] < alpha,
                coef(f3)["z", "Pr(>|z|)"] < alpha,
                coef(f4)["z", "Pr(>|z|)"] < alpha,
                coef(f5)["z", "Pr(>|z|)"] < alpha  , #18
                ## mse
                #mse for logistic regression https://rdrr.io/cran/dvmisc/man/get_mse.html
                sum(zz$residuals^2) / zz$df.residual ,
                sum(zz1$residuals^2) / zz1$df.residual, 
                sum(zz2$residuals^2) / zz2$df.residual ,
                sum(zz3$residuals^2) / zz3$df.residual ,
                sum(zz4$residuals^2) / zz4$df.residual ,
                sum(zz5$residuals^2) / zz5$df.residual ,
                
                (quantile( sum(zz$residuals^2) / zz$df.residual , .025)), #25
                (quantile( sum(zz$residuals^2) / zz$df.residual , .975)), 
                (quantile( sum(zz1$residuals^2) / zz1$df.residual , .025)), #25
                (quantile( sum(zz1$residuals^2) / zz1$df.residual , .975)), 
                (quantile( sum(zz2$residuals^2) / zz2$df.residual , .025)), #25
                (quantile( sum(zz2$residuals^2) / zz2$df.residual , .975)), 
                (quantile( sum(zz3$residuals^2) / zz3$df.residual , .025)), #25
                (quantile( sum(zz3$residuals^2) / zz3$df.residual , .975)), 
                (quantile( sum(zz4$residuals^2) / zz4$df.residual , .025)), #25
                (quantile( sum(zz4$residuals^2) / zz4$df.residual , .975)), 
                (quantile( sum(zz5$residuals^2) / zz5$df.residual , .025)), #25
                (quantile( sum(zz5$residuals^2) / zz5$df.residual , .975)), 
                
                f$aic,  #37
                f1$aic,
                f2$aic,
                f3$aic,
                f4$aic,
                f5$aic ,   #42
                
                # https://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/
                1-(-f$deviance/2)/(-f$null.deviance/2),
                1-(-f1$deviance/2)/(-f1$null.deviance/2),
                1-(-f2$deviance/2)/(-f2$null.deviance/2),
                1-(-f3$deviance/2)/(-f3$null.deviance/2),
                1-(-f4$deviance/2)/(-f4$null.deviance/2),
                1-(-f5$deviance/2)/(-f5$null.deviance/2)
            )
        }
        
        library(plyr)
        res <- raply(simuls, statfun(simfun())) # run the model many times
        # summarize
        result <- apply(res,2,mean)
        q1.result <- apply(res,2, quantile, probs=c(0.025), na.rm=TRUE)
        q2.result <- apply(res,2, quantile, probs=c(0.975), na.rm=TRUE)
        # collect
        return(list(  
            
            res=res,
            result=result ,
            q1.result=q1.result,
            q2.result=q2.result,
            b1=b1
            
        )) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$sim1 <- renderPrint({
        
        return(simul()$result)
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # simulation code 2 do the same as the first simulation code, but this time correlated covariates are created
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    simul2 <- reactive({
        
        sample <- random.sample()
        # need to rename to avoid recursive issues
        K1=sample$K
        Kp=sample$Kp
        pow=sample$pow
        theta1=sample$theta        
        alpha=sample$alpha  
        covar=sample$covar
        Fact=sample$Fact
        p1=sample$p1
        simuls=sample$simuls
        
        N1 <- mcmc()$N # 
        
        b1=simul()$b1  #cal in same beta coefficient as first simulation
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # simulate models many times collect estimate and SE
        
        simfun2 <- function(N=N1, K=K1, a=log(p1/(1-p1)), theta=(theta1), b=b1) {
            # randomi <- runif(N)  
            x <- Matrix(runif(K*K,-RR,RR), K)   # create a correlation matrix randomly , wont allow very high correlations
            
            A <- forceSymmetric(x)
            
            diag(A) <- 1
            
            M <- A
            
            M <- nearPD(M, conv.tol = 1e-7)$mat # default
            # Cholesky decomposition
            L = chol(M)
            nvars = dim(L)[1]
            
            # Random variables that follow an M correlation matrix
            r = t(L) %*% matrix(rnorm(nvars*N, 0,sd1), nrow=nvars, ncol=N)  #2,2
            r = t(r)
            
            r <- as.matrix(r)#
            rdata <- as.data.frame(r)
            XX<- as.matrix(rdata)
            z <- sample(c(0,1), N, replace=T)                # treatment indicator
            lp<-NULL
            lp <- a+ XX %*% b + theta*z     # betas created earlier
            y <- ifelse(runif(N)  < plogis(lp), 1, 0)   # one line
            
            data.frame(X=XX, y=y, z=z)
            
        }
        
        #https://stackoverflow.com/questions/5251507/how-to-succinctly-write-a-formula-with-many-variables-from-a-data-frame
        
        statfun2 <- function(d) {
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            zz <- glm(y~. , data=d,   family = "binomial")      ## adjusting for prognostic X, y2 is not included by use of the '-'
            f <-  summary(zz)
            
            zz1 <- glm(y~z , data=d,   family = "binomial")      ## adjusting for prognostic X, y2 is not included by use of the '-'
            f1 <-  summary(zz1)
            
            cbind(
                coef(f)["z", "Estimate"],
                coef(f)["z", "Std. Error"],
                coef(f1)["z", "Estimate"],
                coef(f1)["z", "Std. Error"], 
                
                coef(f)["z", "Pr(>|z|)"]  < alpha,   
                coef(f1)["z", "Pr(>|z|)"]  < alpha,
                
                ## mse
                sum(zz$residuals^2) / zz$df.residual ,
                sum(zz1$residuals^2) / zz1$df.residual, 
                
                (quantile( sum(zz$residuals^2) / zz$df.residual , .025)), #25
                (quantile( sum(zz$residuals^2) / zz$df.residual , .975)), 
                
                (quantile( sum(zz1$residuals^2) / zz1$df.residual , .025)), #25
                (quantile( sum(zz1$residuals^2) / zz1$df.residual , .975)), 
                
                f$aic,  #37
                f1$aic,
                
                # https://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/
                1-(-f$deviance/2)/(-f$null.deviance/2),
                1-(-f1$deviance/2)/(-f1$null.deviance/2)
                
            )
            
        }
        
        library(plyr)
        res <- raply(simuls, statfun2(simfun2())) # run the model many times
        result <- apply(res,2,mean)
        q1.result <- apply(res,2, quantile, probs=c(0.025), na.rm=TRUE)
        q2.result <- apply(res,2, quantile, probs=c(0.975), na.rm=TRUE)
        
        return(list(  
            
            res=res,
            result=result,  # means
            q1.result=q1.result,
            q2.result=q2.result,
            betas=b1
            
        )) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # do the same as the first simulation code, but this time imbalanced covariates are created
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    simul3 <- reactive({
        
        sample <- random.sample()
        # need to rename to avoid recursive issues
        K1=sample$K
        Kp=sample$Kp
        pow=sample$pow
        theta1=sample$theta        
        alpha=sample$alpha  
        covar=sample$covar
        Fact=sample$Fact
        p1=sample$p1
        simuls=sample$simuls
        
        b1=simul()$b1  #cal in same beta coefficient as first simulation
        
        N1 <- mcmc()$N # 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # simulate models many times collect estimate and SE
        
        simfun3<- function(N=N1, K=K1, a=log(p1/(1-p1)),  theta=(theta1), b=b1) {
            
            MM = N
            N2=MM/2
            N1=N2 
            
            if (covar==1) {  
                X1 <- array(runif(N1*K , -1,1), c(N1,K))  
                X2 <- array(runif(N2*K , -.8,1.2), c(N2,K))   ##imbalance compared to above
                XY <- X <- rbind(X1,X2)
            } else {
                X1 <- array(rnorm(N1*K, 0,  1), c(N1,K))  
                X2 <- array(rnorm(N2*K, .3, 1), c(N2,K))   ##imbalance compared to above
                XY <- X <- rbind(X1,X2)
            }
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #randomi <- runif(N) 
            z <- rep(0:1, c(N1,N2))  #assign 1 so we maintain shift in arms
            lp<-NULL
            lp = a+ X %*% b + theta*z 
            y <- ifelse(runif(N)  < plogis(lp), 1, 0)   # one liner RANDOM!!!
            lp<-NULL
            lp = a + theta*z 
            y2 <- ifelse(runif(N)  < plogis(lp), 1, 0)   # one liner RANDOM!!!
            
            data.frame(X=XY, y=y, z=z, y2=y2)
            
        }
        
        #https://stackoverflow.com/questions/5251507/how-to-succinctly-write-a-formula-with-many-variables-from-a-data-frame
        
        statfun3 <- function(d) {
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            zz <- glm(y~.-y2,  data=d,   family = "binomial")      ## adjusting for prognostic X, y2 is not included by use of the '-'
            f <-  summary(zz)
            
            zz1 <- glm(y~z,  data=d,   family = "binomial")      ## adjusting for prognostic X, y2 is not included by use of the '-'
            f1 <-  summary(zz1)
            
            zz2 <- glm(y2~.-y,  data=d,   family = "binomial")      ## adjusting for prognostic X, y2 is not included by use of the '-'
            f2 <-  summary(zz2)
            
            zz3 <- glm(y2~z,  data=d,   family = "binomial")      ## adjusting for prognostic X, y2 is not included by use of the '-'
            f3 <-  summary(zz3)
            
            cbind(
                coef(f)["z", "Estimate"],
                coef(f)["z", "Std. Error"],
                
                coef(f1)["z", "Estimate"],
                coef(f1)["z", "Std. Error"],
                
                coef(f2)["z", "Estimate"],
                coef(f2)["z", "Std. Error"],
                
                coef(f3)["z", "Estimate"],
                coef(f3)["z", "Std. Error"], #8
                
                # collect p values for power
                coef(f)["z", "Pr(>|z|)"]  < alpha,  #9
                coef(f1)["z", "Pr(>|z|)"] < alpha,
                coef(f2)["z", "Pr(>|z|)"] < alpha,
                coef(f3)["z", "Pr(>|z|)"] < alpha,
                
                ## mse
                #mse for logistic regression https://rdrr.io/cran/dvmisc/man/get_mse.html
                sum(zz$residuals^2) / zz$df.residual ,
                sum(zz1$residuals^2) / zz1$df.residual, 
                sum(zz2$residuals^2) / zz2$df.residual ,
                sum(zz3$residuals^2) / zz3$df.residual ,
                
                
                (quantile( sum(zz$residuals^2) / zz$df.residual , .025)), #25
                (quantile( sum(zz$residuals^2) / zz$df.residual , .975)), 
                
                (quantile( sum(zz1$residuals^2) / zz1$df.residual , .025)), #25
                (quantile( sum(zz1$residuals^2) / zz1$df.residual , .975)), 
                
                (quantile( sum(zz2$residuals^2) / zz2$df.residual , .025)), #25
                (quantile( sum(zz2$residuals^2) / zz2$df.residual , .975)), 
                
                (quantile( sum(zz3$residuals^2) / zz3$df.residual , .025)), #25
                (quantile( sum(zz3$residuals^2) / zz3$df.residual , .975)), 
                
                f$aic,  #37
                f1$aic,
                f2$aic,
                f3$aic,
                
                # https://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/
                1-(-f$deviance/2)/(-f$null.deviance/2),
                1-(-f1$deviance/2)/(-f1$null.deviance/2),
                1-(-f2$deviance/2)/(-f2$null.deviance/2),
                1-(-f3$deviance/2)/(-f3$null.deviance/2) 
                
            )
            
        }
        
        library(plyr)
        res <- raply(simuls, statfun3(simfun3())) # run the model many times
        # summarize
        result <- apply(res,2,mean)
        q1.result <- apply(res,2, quantile, probs=c(0.025), na.rm=TRUE)
        q2.result <- apply(res,2, quantile, probs=c(0.975), na.rm=TRUE)
        # collect
        return(list(  
            
            res=res,
            result=result ,
            q1.result=q1.result,
            q2.result=q2.result,
            b1=b1
            
        )) 
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # simulation plots
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # collect simulation trt effect estimates from simulation and plot!
    
    output$reg.plotxx <- output$reg.plotx <- renderPlot({         #means
        
        # Get the  data from the three simulations
        res <- simul()$res
        res2 <- simul2()$res
        res3 <- simul3()$res
        
        sample <- random.sample()
        theta1=sample$theta        # get true trt effect
        
        d1 <-  density(res[,1] )
        d2 <-  density(res[,3] )
        d3 <-  density(res[,5] )
        d4 <-  density(res[,7] )
        d5 <-  density(res[,9] )
        d6 <-  density(res[,11] )
        d7 <-  density(res2[,1] )
        d8 <-  density(res2[,3] )
        d9 <-   density(res3[,1] )
        d10 <-  density(res3[,3] )
        d11 <-  density(res3[,5] )
        d12 <-  density(res3[,7] )
        
        dz <- max(c(d1$y, d2$y, d3$y, d4$y, d5$y, d6$y, d7$y, d8$y  , d9$y, d10$y, d11$y, d12$y  ))
        dx <- range(c(d1$x,d2$x,  d3$x, d4$x, d5$x, d6$x, d7$x, d8$x   , d9$x, d10$x, d11$x, d12$x  ))
        
        if (input$dist %in% "All") {
            
            plot((d1), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww,
                 xlab="Treatment effect log odds",  
                 ylab="Density")                           
            lines( (d2), col = "black", lty=w, lwd=ww)  
            lines( (d3), col = "red", lty=wz, lwd=ww)    
            lines( (d4), col = "red", lty=w, lwd=ww)          
            lines( (d5), col = "blue", lty=wz, lwd=ww)       
            lines( (d6), col = "blue", lty=w, lwd=ww)       
            lines( (d7), col = "purple", lty=wz, lwd=ww)       
            lines( (d8), col = "purple", lty=w, lwd=ww)       
            
            lines( (d9), col = "green", lty=wz, lwd=ww)       
            lines( (d10), col = "green", lty=w, lwd=ww)       
            lines( (d11), col = "grey", lty=wz, lwd=ww)       
            lines( (d12), col = "grey", lty=w, lwd=ww)  
            
        }
        
        else if (input$dist %in% "d1") {  #remove
            
            plot((d1), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww,
                 xlab="Treatment effect log odds",  
                 ylab="Density")  
            lines( (d2), col = "black", lty=w, lwd=ww)  
            
        }
        
        else if (input$dist %in% "d3") {  #remove
            
            plot((d3), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww,col="red",
                 xlab="Treatment effect log odds",  
                 ylab="Density")               
            lines( (d4), col = "red", lty=w, lwd=ww)          
            
        }
        
        else if (input$dist %in% "d5") {
            
            plot((d5), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww, col="blue",
                 xlab="Treatment effect log odds",  
                 ylab="Density")                    
            
            lines( (d6), col = "blue", lty=w, lwd=ww)       
            
        }
        
        else if (input$dist %in% "d7") {
            
            plot((d7), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww, col="purple",
                 xlab="Treatment effect log odds",  
                 ylab="Density") 
            
            lines( (d8), col = "purple", lty=w, lwd=ww)     
            
        }
        else if (input$dist %in% "d9") {
            
            plot((d9), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww, col="green",
                 xlab="Treatment effect log odds",  
                 ylab="Density")  
            
            lines( (d10), col = "green", lty=w, lwd=ww)     
            
        }
        
        else if (input$dist %in% "d11") {
            
            plot((d11), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww, col="grey",
                 xlab="Treatment effect log odds",  
                 ylab="Density")  
            
            lines( (d12), col = "grey", lty=w, lwd=ww)     
            
        }
        
        abline(v = theta1, col = "darkgrey")                
        legend("topright",       # Add legend to density
               legend = c(" adj. for true prognostic covariates", 
                          " not adj. for true prognostic covariates" ,
                          " adj. for covariates unrelated to outcome", 
                          " not adj. for covariates unrelated to outcome",
                          " adj. for mix of prognostic and unrelated to outcome", 
                          " not adj. mix of prognostic and unrelated to outcome", 
                          " adj. for correlated prognostic covariates", 
                          " not adj. for correlated prognostic covariates",
                          " adj. for imbalanced prognostic covariates", 
                          " not adj. for imbalanced prognostic covariates", 
                          " adj. for imbalanced covariates unrelated to outcome", 
                          " not adj. imbalanced covariates unrelated to outcome"
                          
               ),
               col = c("black", "black","red","red","blue", "blue", "purple", "purple", "green", "green", "grey", "grey"),
               lty = c(wz, w,wz,w,wz,w,wz,w,wz,w,wz,w)  ,lwd=ww
               , bty = "n", cex=1)
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    
    # Need to be able to calculate logistic regression se analytically base on p1 and OR?
    # so for now we do simulation, not ideal
    
    # simulx <- reactive({
    #   
    #   sample <- random.sample()
    #   # need to rename to avoid recursive issues
    #   K1=sample$K
    #   Kp=sample$Kp
    #   pow=sample$pow
    #   theta1=sample$theta        
    #   alpha=sample$alpha  
    #   covar=sample$covar
    #   Fact=sample$Fact
    #   p1=sample$p1
    #   simuls=sample$simuls
    # 
    #   N<- mcmc()$N
    # 
    #   ## function to simulate simple dataset #####################################################################
    #   
    #   fun.d<-function(nsample, drug.allocation, 
    #                   alpha,  beta.drug,
    #                   seed=NULL){ 
    #     
    #     if (!is.null(seed)) set.seed(seed)
    #     
    #     drug<- (rbinom(nsample, 1, prob =drug.allocation ))   
    #     
    #     Xmat <- model.matrix(~ drug )
    #     beta.vec <- c(alpha,  beta.drug )
    #     
    #     lin.pred <- Xmat[,] %*% beta.vec                 # Value of lin.predictor
    #     exp.p <- exp(lin.pred) / (1 + exp(lin.pred))     # Expected proportion
    #     y <- rbinom(n = nsample, size = 1, prob = exp.p) # Add binomial noise
    #     #y<- runif(nsample) <  exp.p                     # alternatively ads noise in this way
    #     
    #     d<-as.data.frame(cbind(y, drug))         # create a dataset
    #     
    #     return(d)
    #     
    #   }
    #   
    #  # function to pull out se######################################################################################
    #   simfunc <- function(d) {
    #     fit1 <- glm(y  ~ drug , d, family = binomial) 
    #     c( summary(fit1)$coef["drug","Std. Error"] )   # collect se
    #   }
    #   
    #  ################################################################################################################
    #  out <- replicate(1000, simfunc(fun.d( nsample=N, drug.allocation=0.5,    
    #                                       alpha=log(p1/(1-p1)),  # turn prob into odds and enter this here
    #                                       beta.drug=(theta1))))
    # 
    # ################################################################################################################
    #  se. <- mean(out)
    # 
    #    
    # return(list(   
    #   
    #   se.=se. 
    # ))
    
    ###############################################################################################################
    # True standard error 
    # not relying on simulation
    # logistic regression se by hand for log odds ratio of drug effect based on p1 and OR from which N is calculated
    ############################################################################################################### 
    simulx <- reactive({
        
        sample <- random.sample()
        
        pow=sample$pow
        theta1=sample$theta        
        alpha=sample$alpha  
        p1=sample$p1
        
        N<- mcmc()$N
        
        (post.odds <- p1/(1-p1) *  (theta1))          # calculate post odds
        (p2 <-post.prob <- post.odds/(post.odds+1))  # calculate post prob
        
        # calcualte N in each cell of 2by2 table
        cellA <-  (1-p1)*N/2
        cellC	<- p1*N/2
        cellB	<- (1-p2)*N/2
        cellD	<- p2*N/2
        
        # now se of log odds ratio
        se. <- log.odds.se. <- sqrt((1/cellA)+(1/cellB)+(1/cellC)+(1/cellD))
        
        return(list(   
            se.=se. 
        ))
        
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # SIMUALTION PLOT
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # collect simulatio standard error estimates from simulation and plot!
    
    output$reg.plotyy <- output$reg.ploty <- renderPlot({         #standard errors
        
        # Get the  data
        
        res <- simul()$res
        res2 <- simul2()$res
        res3 <- simul3()$res
        se. <- simulx()$se.
        
        
        
        sample <- random.sample()
        
        N1 <- mcmc()$N # 
        
        n1 <- mcmc()$Na
        n2 <- mcmc()$Nb
        
        theta<-sample$theta # not needed
        
        zz <- table.sim()$zz
        
        #######################
        # here we save simulation results
        # rename if important to keep
        save(list = c("wz","w","ww","se.","N1","n1","n2","res", "res2","res3","theta","zz"), file = "simulation_results.Rdata")  
        #######################
        
        d1 <-  density(res[,2] )
        d2 <-  density(res[,4] )
        d3 <-  density(res[,6] )
        d4 <-  density(res[,8] )
        d5 <-  density(res[,10] )
        d6 <-  density(res[,12] )
        d7 <-  density(res2[,2] )
        d8 <-  density(res2[,4] )
        
        d9 <-   density(res3[,2] )
        d10 <-  density(res3[,4] )
        d11 <-  density(res3[,6] )
        d12 <-  density(res3[,8] )
        
        
        dz <- max(c(d1$y, d2$y, d3$y, d4$y, d5$y, d6$y, d7$y, d8$y  , d9$y, d10$y, d11$y, d12$y    ))
        dx <- range(c(d1$x,d2$x,  d3$x, d4$x, d5$x, d6$x, d7$x, d8$x   , d9$x, d10$x, d11$x, d12$x    ))
        
        if (input$dist %in% "All") {
            
            plot( (d1), xlim = c(dx), main=paste0("Density of treatment standard error estimates, truth= ",p4(se.),""), ylim=c(0,dz),lty=wz, lwd=ww,
                  xlab="Standard error of log odds trt effect",  
                  ylab="Density")  
            lines( (d2), col = "black", lty=w, lwd=ww)  
            lines( (d3), col = "red", lty=wz, lwd=ww)    
            lines( (d4), col = "red", lty=w, lwd=ww)          
            lines( (d5), col = "blue", lty=wz, lwd=ww)       
            lines( (d6), col = "blue", lty=w, lwd=ww)       
            lines( (d7), col = "purple", lty=wz, lwd=ww)       
            lines( (d8), col = "purple", lty=w, lwd=ww)       
            
            lines( (d9), col = "green", lty=wz, lwd=ww)       
            lines( (d10), col = "green", lty=w, lwd=ww)       
            lines( (d11), col = "grey", lty=wz, lwd=ww)       
            lines( (d12), col = "grey", lty=w, lwd=ww)  
            
        }
        
        
        else if (input$dist %in% "d1") {  
            
            plot((d1), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww,
                 xlab="Standard error of log odds trt effect",  
                 ylab="Density") 
            lines( (d2), col = "black", lty=w, lwd=ww)  
            
        }
        
        else if (input$dist %in% "d3") {  
            
            
            plot((d3), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww,col="red",
                 xlab="Standard error of log odds trt effect",  
                 ylab="Density")  
            lines( (d4), col = "red", lty=w, lwd=ww)          
            
        }
        
        else if (input$dist %in% "d5") {
            
            plot((d5), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww, col="blue",
                 xlab="Standard error of log odds trt effect",  
                 ylab="Density")  
            
            lines( (d6), col = "blue", lty=w, lwd=ww)       
            
        }
        
        else if (input$dist %in% "d7") {
            
            plot((d7), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww, col="purple",
                 xlab="Standard error of log odds trt effect",  
                 ylab="Density") 
            
            lines( (d8), col = "purple", lty=w, lwd=ww)     
            
        }
        
        else if (input$dist %in% "d9") {
            
            plot((d9), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww, col="green",
                 xlab="Standard error of log odds trt effect",  
                 ylab="Density")  
            
            lines( (d10), col = "green", lty=w, lwd=ww)     
            
        }
        
        else if (input$dist %in% "d11") {
            
            plot((d11), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww, col="grey",
                 xlab="Standard error of log odds trt effect",  
                 ylab="Density")  
            
            lines( (d12), col = "grey", lty=w, lwd=ww)     
            
        }
        
        abline(v = se., col = "darkgrey")   
        legend("topright",           # Add legend to density
               legend = c(" adj. for true prognostic covariates", 
                          " not adj. for true prognostic covariates" ,
                          " adj. for covariates unrelated to outcome", 
                          " not adj. for covariates unrelated to outcome",
                          " adj. for mix of prognostic and unrelated to outcome", 
                          " not adj. mix of prognostic and unrelated to outcome", 
                          " adj. for correlated prognostic covariates", 
                          " not adj. for correlated prognostic covariates",
                          " adj. for imbalanced prognostic covariates", 
                          " not adj. for imbalanced prognostic covariates", 
                          " adj. for imbalanced covariates unrelated to outcome", 
                          " not adj. imbalanced covariates unrelated to outcome"
                          
               ),
               col = c("black", "black","red","red","blue", "blue", "purple", "purple", "green", "green", "grey", "grey"),
               lty = c(wz, w,wz,w,wz,w,wz,w,wz,w,wz,w) ,lwd=ww
               , bty = "n", cex=1)
    })
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # not used for binary logistic
    output$textWithNumber99 <- renderText({ 
        
        HTML(
            "Mean squared error (MSE: accuracy and precision) combines bias and
                    variance as (bias*bias+variance). It represents the total variation around the
                    true value, rather than the average estimated value. MSE gives an overall sense of the quality of the
                    estimator. As the MSE can be written as the sum of the variance of the estimator and the squared bias of the estimator, 
                    this implies that in the case of unbiased estimators, the MSE and variance are equivalent. So compare the calculated MSE to the 
                    true sigma squared 
                    on the left input and printed here:"
        )
        
    })  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # table for simulation summary, quite easy to make a mistake here! too much hard coding really, need to code better
    table.sim <- reactive({
        
        sample <- random.sample()
        theta1=sample$theta       
        
        res <- simul()$res  
        result <- simul()$result  
        result2 <- simul2()$result  
        result3 <- simul3()$result  
        
        q1.result <- simul()$q1.result  
        q2.result <- simul()$q2.result  
        
        q1.result2 <- simul2()$q1.result  
        q2.result2 <- simul2()$q2.result  
        
        q1.result3 <- simul3()$q1.result  
        q2.result3 <- simul3()$q2.result  
        
        zz <- rbind(
            (c( p4(result[1])   ,     p2(q1.result[1])  ,  p2(q2.result[1])   , p4(result[2] ) ,  p4(result[13] ) ,  p4(result[1] -theta1) ,      p4(result[37] )    ,  p4(result[43] )         )) ,
            (c( p4(result[3])   ,     p2(q1.result[3]) ,   p2(q2.result[3])   , p4(result[4] ) ,  p4(result[14] ) ,  p4(result[3] -theta1) ,      p4(result[38] )    ,  p4(result[44] )        )) ,
            (c( p4(result[5])   ,     p2(q1.result[5]) ,   p2(q2.result[5])   , p4(result[6] ) ,  p4(result[15] ) ,  p4(result[5] -theta1) ,      p4(result[39] )    ,  p4(result[45] )        )) ,
            (c( p4(result[7])   ,     p2(q1.result[7]) ,   p2(q2.result[7])   , p4(result[8] ) ,  p4(result[16] ) ,  p4(result[7] -theta1) ,      p4(result[40] )    ,  p4(result[46] )         )) ,
            (c( p4(result[9])   ,     p2(q1.result[9]) ,   p2(q2.result[9])   , p4(result[10] ) , p4(result[17] ) ,  p4(result[9] -theta1) ,      p4(result[41] )    ,  p4(result[47] )           )) ,
            (c( p4(result[11])  ,     p2(q1.result[11]) ,  p2(q2.result[11])  , p4(result[12] ) , p4(result[18] ) ,  p4(result[11] -theta1) ,      p4(result[42] )    ,  p4(result[48] )         )) ,
            (c( p4(result2[1])  ,     p2(q1.result2[1]),   p2(q2.result2[1])  , p4(result2[2] ) , p4(result2[5] ) ,  p4(result2[1] -theta1) ,      p4(result2[13] )   ,  p4(result2[15] )         )) ,
            (c( p4(result2[3])  ,     p2(q1.result2[3])  , p2(q2.result2[3])  , p4(result2[4] ) , p4(result2[6] ) ,  p4(result2[3] -theta1) ,      p4(result2[14] )   ,  p4(result2[16] )          )),
            (c( p4(result3[1])  ,     p2(q1.result3[1])  , p2(q2.result3[1])   , p4(result3[2] ) ,  p4(result3[9] ) ,  p4(result3[1] -theta1) , p4(result3[25] )   ,  p4(result3[29] )         )) ,
            (c( p4(result3[3])  ,     p2(q1.result3[3]) ,  p2(q2.result3[3])   , p4(result3[4] ) ,  p4(result3[10] ) ,  p4(result3[3] -theta1) ,p4(result3[26] )   ,  p4(result3[30] )         )) ,
            (c( p4(result3[5])  ,     p2(q1.result3[5]) ,  p2(q2.result3[5])   , p4(result3[6] ) ,  p4(result3[11] ) ,  p4(result3[5] -theta1) , p4(result3[27] )  ,  p4(result3[31] )         )) ,
            (c( p4(result3[7])  ,     p2(q1.result3[7]) ,  p2(q2.result3[7])   , p4(result3[8] ) ,  p4(result3[12] ) ,  p4(result3[7] -theta1) , p4(result3[28] )  ,  p4(result3[32]  )         )) 
        ) 
        
        zz <- as.data.frame(zz)
        
        colnames(zz) <- c("Mean  ", "Lower 95%CI", "Upper 95%CI", "Stand.error", "Power ","B" , "sigma","R2")
        
        zz <- data.frame(lapply(zz, function(x) as.numeric(as.character(x))))
        zz <- as.data.frame(zz)
        rownames(zz)<- c(
            " adj. for true prognostic covariates", 
            " not adj. for true prognostic covariates" ,
            " adj. for non prognostic covariates", 
            " not adj. for non prognostic covariates",
            " adj. for some non prognostic covariates", 
            " not adj. when some prognostic covariates", 
            " adj. for correlated prognostic covariates", 
            " not adj. for correlated prognostic covariates",
            " adj. for imbalanced prognostic covariates",
            " not adj. for imbalanced prognostic covariates",  
            " adj. for non prognostic imbalanced covariates",
            " not adj. for non prognostic imbalanced covariates"
        )
        zz <- zz[order(zz$B),]
        
        
        colnames(zz) <- c("Mean  ", "Lower 95%CI", "Upper 95%CI", "Std.error", "Power ","bias" , "AIC","McFadden's R2")
        zz <- zz[order(zz$AIC),]
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        return(list(  
            
            zz=zz
            
        )) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    })
    
    output$zzz <- output$zz <- renderPrint({
        
        d <- table.sim()$zz
        
        return(d)
    })
    ##~~~~~~~~~~~~~
    
    output$textWithNumber1a <- renderText({ 
        
        placebo <- mcmc()$placebo
        treated <- mcmc()$treated
        bigN <- mcmc()$bigN
        
        HTML(paste0(  #tags$hr(),
            "Figure 1 Simulation results. Randomised 1:1, we have  "  
            
            ,tags$span(style="color:red",  bigN  ),
            " total patients randomised 1:1 for each simulation. The true covariate coefficients are fixed at the same values for all simulations
                      and are selected randomly between +/- multiples of the treatment effect, as dictated by the input on left. The true covariate coefficients are printed at the bottom."
            
        ))    
        
    })
    
    
    output$betas <- renderPrint({
        
        d <- simul2()$betas
        return(print(d))
        
    })
    
    #####################################################################################################################
    # Here we code up how to click a button and automatically upload rdata file
    # tough to code again but duelling buttons link v helpful:  # https://shiny.rstudio.com/articles/action-buttons.html
    #####################################################################################################################
    
    # declare empty objects to populate
    content1 <- reactiveValues(tab1 = NULL)
    content2 <- reactiveValues(tab2 = NULL)
    content3 <- reactiveValues(tab3 = NULL)
    content4 <- reactiveValues(tab4 = NULL)
    content5 <- reactiveValues(tab5 = NULL)
    content6 <- reactiveValues(tab6 = NULL)
    
    # If upload button is pressed this will activate
    observeEvent(input$upload,    
                 
                 isolate({
                     isfar <-  load(url(pp))
                     
                     content1$tab1 <-  get((isfar)[12])  # fill object with specific Rdata dataset, this is the summary table  
                     content2$tab2 <-  get((isfar)[8])   # as above this is the res dataset
                     content3$tab3 <-  get((isfar)[9])
                     content4$tab4 <-  get((isfar)[10])
                     content5$tab5 <-  get((isfar)[11])
                     content6$tab6 <-  get((isfar)[4])
                 })
                 
    )
    
    # If upload2 button is pressed this will activate
    observeEvent(input$upload2, 
                 
                 isolate({
                     isfar <-  load(url(pp2))  # 2nd link
                     
                     content1$tab1 <-  get((isfar)[12])  # fill object with specific Rdata dataset, this is the summary table 
                     content2$tab2 <-  get((isfar)[8])   # as above this is the res dataset
                     content3$tab3 <-  get((isfar)[9])   
                     content4$tab4 <-  get((isfar)[10])
                     content5$tab5 <-  get((isfar)[11])
                     content6$tab6 <-  get((isfar)[4])
                 })
                 
    )
    
    # If upload3 button is pressed this will activate
    observeEvent(input$upload3, 
                 
                 isolate({
                     isfar <-  load(url(pp3))  # 3rdd link
                     
                     content1$tab1 <-  get((isfar)[12])  # fill object with specific Rdata dataset, this is the summary table 
                     content2$tab2 <-  get((isfar)[8])   # as above this is the res dataset
                     content3$tab3 <-  get((isfar)[9])
                     content4$tab4 <-  get((isfar)[10])
                     content5$tab5 <-  get((isfar)[11])
                     content6$tab6 <-  get((isfar)[4])
                 })
                 
    )
    
    # If upload3 button is pressed this will activate
    observeEvent(input$upload4, 
                 
                 isolate({
                     isfar <-  load(url(pp4))  # 4th link
                     
                     content1$tab1 <-  get((isfar)[12])  # fill object with specific Rdata dataset, this is the summary table 
                     content2$tab2 <-  get((isfar)[8])   # as above this is the res dataset
                     content3$tab3 <-  get((isfar)[9])
                     content4$tab4 <-  get((isfar)[10])
                     content5$tab5 <-  get((isfar)[11])
                     content6$tab6 <-  get((isfar)[4])
                     
                     
                     
                     
                 })
                 
    )
    
    # now we have put the data that we load into objects that can be used as inputs  
    
    # observeEvent(input$upload, {
    #     output$content1 <- renderPrint({
    #         if (is.null(content1$tab1)) return()
    #         content1$tab1
    #     })
    # })
    # 
    # output$content2 <- renderPrint({
    #     if (is.null(content2$tab2)) return()
    #     content2$tab2
    # })
    # output$content3 <- renderPrint({
    #     if (is.null(content3$tab3)) return()
    #     content3$tab3
    # })
    # output$content4 <- renderPrint({
    #     if (is.null(content4$tab4)) return()
    #     content4$tab4
    # })
    # output$content5 <- renderPrint({
    #     if (is.null(content5$tab5)) return()
    #     content5$tab5
    # })    
    # output$content6 <- renderPrint({
    #     if (is.null(content6$tab6)) return()
    #     content6$tab6
    # })
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  not working
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # values <- reactiveValues(one=0 )
    # 
    # observeEvent(input$upload, {
    #     values$one <- 1
    # })
    # observeEvent(input$upload2, {
    #     values$one <- 1
    # })
    # observeEvent(input$upload3, {
    #     values$one <- 1
    # })
    # observeEvent(input$upload4, {
    #     values$one <- 1
    # })
    # #################
    # end not working
    ##################
    #  https://groups.google.com/g/shiny-discuss/c/vd_nB-BH8sw
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # collect simulation trt effect estimates from upload dnd plot!
    # New function to plot basically copy of earlier function
    
    
      
    observeEvent(c(input$upload, input$upload2, input$upload3, input$upload4),{  
   
        output$content1 <- renderPrint({
            if (is.null(content1$tab1)) return()
            content1$tab1
        })
        
        output$content2 <- renderPrint({
            if (is.null(content2$tab2)) return()
            content2$tab2
        })
        output$content3 <- renderPrint({
            if (is.null(content3$tab3)) return()
            content3$tab3
        })
        output$content4 <- renderPrint({
            if (is.null(content4$tab4)) return()
            content4$tab4
        })
        output$content5 <- renderPrint({
            if (is.null(content5$tab5)) return()
            content5$tab5
        })    
        output$content6 <- renderPrint({
            if (is.null(content6$tab6)) return()
            content6$tab6
        })
        
        
       output$plot1   <- renderPlot({         #means    
        
        # pull in the loaded objects
        if (is.null(content2$tab2)) return()  # this stops red error messages before the first button is loaded
        if (is.null(content3$tab3)) return()
        if (is.null(content4$tab4)) return()
        if (is.null(content5$tab5)) return()
        
        res <- as.data.frame(content2$tab2)     # loaded objects assigned to objects
        res <- as.data.frame(lapply(res, as.numeric))
        
        res2 <- as.data.frame(content3$tab3)
        res2 <- as.data.frame(lapply(res2, as.numeric))
        
        res3 <- as.data.frame(content4$tab4)
        res3 <- as.data.frame(lapply(res3, as.numeric))
        
        theta1 <- (content5$tab5)
        
        ## below here code is the same 
        
        sample <- random.sample()
        
        d1 <-  density( res[,1])
        d2 <-  density(res[,3] )
        d3 <-  density(res[,5] )
        d4 <-  density(res[,7] )
        d5 <-  density(res[,9] )
        d6 <-  density(res[,11] )
        
        d7 <-  density(res2[,1] )
        d8 <-  density(res2[,3] )
        
        d9 <-   density(res3[,1] )
        d10 <-  density(res3[,3] )
        d11 <-  density(res3[,5] )
        d12 <-  density(res3[,7] )
        
        dz <- max(c(d1$y, d2$y, d3$y, d4$y, d5$y, d6$y, d7$y, d8$y  , d9$y, d10$y, d11$y, d12$y  ))
        dx <- range(c(d1$x,d2$x,  d3$x, d4$x, d5$x, d6$x, d7$x, d8$x   , d9$x, d10$x, d11$x, d12$x  ))
        
        if (input$dist %in% "All") {
            
            plot((d1), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww,
                 xlab="Treatment effect log odds",
                 ylab="Density")
            lines( (d2), col = "black", lty=w, lwd=ww)
            lines( (d3), col = "red", lty=wz, lwd=ww)
            lines( (d4), col = "red", lty=w, lwd=ww)
            lines( (d5), col = "blue", lty=wz, lwd=ww)
            lines( (d6), col = "blue", lty=w, lwd=ww)
            lines( (d7), col = "purple", lty=wz, lwd=ww)
            lines( (d8), col = "purple", lty=w, lwd=ww)
            
            lines( (d9), col = "green", lty=wz, lwd=ww)
            lines( (d10), col = "green", lty=w, lwd=ww)
            lines( (d11), col = "grey", lty=wz, lwd=ww)
            lines( (d12), col = "grey", lty=w, lwd=ww)
            
        }
        
        else if (input$dist %in% "d1") {  #remove
            
            
            plot((d1), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww,
                 xlab="Treatment effect log odds",
                 ylab="Density")
            lines( (d2), col = "black", lty=w, lwd=ww)
            
        }
        
        else if (input$dist %in% "d3") {  #remove
            
            
            plot((d3), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww,col="red",
                 xlab="Treatment effect log odds",
                 ylab="Density")
            lines( (d4), col = "red", lty=w, lwd=ww)
            
        }
        
        else if (input$dist %in% "d5") {
            
            plot((d5), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww, col="blue",
                 xlab="Treatment effect log odds",
                 ylab="Density")
            
            lines( (d6), col = "blue", lty=w, lwd=ww)
            
        }
        
        else if (input$dist %in% "d7") {
            
            plot((d7), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww, col="purple",
                 xlab="Treatment effect log odds",
                 ylab="Density")
            
            lines( (d8), col = "purple", lty=w, lwd=ww)
            
        }
        else if (input$dist %in% "d9") {
            
            plot((d9), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww, col="green",
                 xlab="Treatment effect log odds",
                 ylab="Density")
            
            lines( (d10), col = "green", lty=w, lwd=ww)
            
        }
        
        else if (input$dist %in% "d11") {
            
            plot((d11), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww, col="grey",
                 xlab="Treatment effect log odds",
                 ylab="Density")
            
            lines( (d12), col = "grey", lty=w, lwd=ww)
            
        }
        
        abline(v = theta1, col = "darkgrey")
        legend("topright",       # Add legend to density
               legend = c(" adj. for true prognostic covariates",
                          " not adj. for true prognostic covariates" ,
                          " adj. for covariates unrelated to outcome",
                          " not adj. for covariates unrelated to outcome",
                          " adj. for mix of prognostic and unrelated to outcome",
                          " not adj. mix of prognostic and unrelated to outcome",
                          " adj. for correlated prognostic covariates",
                          " not adj. for correlated prognostic covariates",
                          " adj. for imbalanced prognostic covariates",
                          " not adj. for imbalanced prognostic covariates",
                          " adj. for imbalanced covariates unrelated to outcome",
                          " not adj. imbalanced covariates unrelated to outcome"
                          
               ),
               col = c("black", "black","red","red","blue", "blue", "purple", "purple", "green", "green", "grey", "grey"),
               lty = c(wz, w,wz,w,wz,w,wz,w,wz,w,wz,w)  ,lwd=ww
               , bty = "n", cex=1)
    })
    
    
    
    
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # We do the same but for the se estimates
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
     output$plot2<- renderPlot({        
         
        
        #################################
        if (is.null(content2$tab2)) return()
        if (is.null(content3$tab3)) return()
        if (is.null(content4$tab4)) return()
        if (is.null(content6$tab6)) return()
        
        res <- as.data.frame(content2$tab2)
        res <- as.data.frame(lapply(res, as.numeric))
        
        res2 <- as.data.frame(content3$tab3)
        res2 <- as.data.frame(lapply(res2, as.numeric))
        
        res3 <- as.data.frame(content4$tab4)
        res3 <- as.data.frame(lapply(res3, as.numeric))
        
        se. <- (content6$tab6)                     ## this is the only difference to previous plot function, here we pull in se
        #################################
        
        sample <- random.sample()
        
        d1 <-  density(res[,2] )
        d2 <-  density(res[,4] )
        d3 <-  density(res[,6] )
        d4 <-  density(res[,8] )
        d5 <-  density(res[,10] )
        d6 <-  density(res[,12] )
        
        d7 <-  density(res2[,2] )
        d8 <-  density(res2[,4] )
        
        d9 <-   density(res3[,2] )
        d10 <-  density(res3[,4] )
        d11 <-  density(res3[,6] )
        d12 <-  density(res3[,8] )
        
        dz <- max(c(d1$y, d2$y, d3$y, d4$y, d5$y, d6$y, d7$y, d8$y  , d9$y, d10$y, d11$y, d12$y    ))
        dx <- range(c(d1$x,d2$x,  d3$x, d4$x, d5$x, d6$x, d7$x, d8$x   , d9$x, d10$x, d11$x, d12$x    ))
        
        if (input$dist %in% "All") {
            
            plot( (d1), xlim = c(dx), main=paste0("Density of treatment standard error estimates, truth= ",p4(se.),""), ylim=c(0,dz),lty=wz, lwd=ww,
                  xlab="Standard error of log odds trt effect",  
                  ylab="Density")  
            lines( (d2), col = "black", lty=w, lwd=ww)  
            lines( (d3), col = "red", lty=wz, lwd=ww)    
            lines( (d4), col = "red", lty=w, lwd=ww)          
            lines( (d5), col = "blue", lty=wz, lwd=ww)       
            lines( (d6), col = "blue", lty=w, lwd=ww)       
            lines( (d7), col = "purple", lty=wz, lwd=ww)       
            lines( (d8), col = "purple", lty=w, lwd=ww)       
            
            lines( (d9), col = "green", lty=wz, lwd=ww)       
            lines( (d10), col = "green", lty=w, lwd=ww)       
            lines( (d11), col = "grey", lty=wz, lwd=ww)       
            lines( (d12), col = "grey", lty=w, lwd=ww)  
            
        }
        
        
        else if (input$dist %in% "d1") {  
            
            plot((d1), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww,
                 xlab="Standard error of log odds trt effect",  
                 ylab="Density") 
            lines( (d2), col = "black", lty=w, lwd=ww)  
            
        }
        
        else if (input$dist %in% "d3") {  
            
            
            plot((d3), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww,col="red",
                 xlab="Standard error of log odds trt effect",  
                 ylab="Density")  
            lines( (d4), col = "red", lty=w, lwd=ww)          
            
        }
        
        else if (input$dist %in% "d5") {
            
            plot((d5), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww, col="blue",
                 xlab="Standard error of log odds trt effect",  
                 ylab="Density")  
            
            lines( (d6), col = "blue", lty=w, lwd=ww)       
            
        }
        
        else if (input$dist %in% "d7") {
            
            plot((d7), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww, col="purple",
                 xlab="Standard error of log odds trt effect",  
                 ylab="Density") 
            
            lines( (d8), col = "purple", lty=w, lwd=ww)     
            
        }
        
        else if (input$dist %in% "d9") {
            
            plot((d9), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww, col="green",
                 xlab="Standard error of log odds trt effect",  
                 ylab="Density")  
            
            lines( (d10), col = "green", lty=w, lwd=ww)     
            
        }
        
        else if (input$dist %in% "d11") {
            
            plot((d11), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww, col="grey",
                 xlab="Standard error of log odds trt effect",  
                 ylab="Density")  
            
            lines( (d12), col = "grey", lty=w, lwd=ww)     
            
        }
        
        abline(v = se., col = "darkgrey")   
        legend("topright",           # Add legend to density
               legend = c(" adj. for true prognostic covariates", 
                          " not adj. for true prognostic covariates" ,
                          " adj. for covariates unrelated to outcome", 
                          " not adj. for covariates unrelated to outcome",
                          " adj. for mix of prognostic and unrelated to outcome", 
                          " not adj. mix of prognostic and unrelated to outcome", 
                          " adj. for correlated prognostic covariates", 
                          " not adj. for correlated prognostic covariates",
                          " adj. for imbalanced prognostic covariates", 
                          " not adj. for imbalanced prognostic covariates", 
                          " adj. for imbalanced covariates unrelated to outcome", 
                          " not adj. imbalanced covariates unrelated to outcome"
                          
               ),
               col = c("black", "black","red","red","blue", "blue", "purple", "purple", "green", "green", "grey", "grey"),
               lty = c(wz, w,wz,w,wz,w,wz,w,wz,w,wz,w) ,lwd=ww
               , bty = "n", cex=1)
    })
    
    ##################################################################################################################################################   
    # Here we start coding to browse (coded out as we are not using). We start coding how to automatically load in Rdata from pre run simulations
    # took a while to get this right. 
    # https://github.com/rstudio/shiny-examples/tree/master/066-upload-file
    # v tough to get this working !!! this is calling in Rdata and loading particular data from the Rdata
    # load(url(pp<- "https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/A%205000%20default%20settings%20theta%20log1.5%20-1.00%20-0.67%20-0.43.Rdata"))
    # Initial idea to browse to Rdata, but it is better to click and auto load so I dropped this idea but leaving code here as a resource as it works
    ##################################################################################################################################################   
    #    output$contents3 <- renderTable({  # THIS IS THE SUMMARY TABLE
    #        
    #        inFile <- input$file1
    #        if (is.null(inFile))
    #            return(NULL)
    #       # return(validate(need(input$file1, ""))) #I'm sending an empty string as message
    #        isfar <- (load(inFile$datapath))
    #        (get((isfar)[12]))
    #        
    #    }, rownames = TRUE ,striped=TRUE, bordered=TRUE)
    #    
    # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #    trt.effect1  <- reactive({  #RES
    #        
    #        inFile <- input$file1
    #        if (is.null(inFile))
    #            return(NULL)
    #        isfar <- (load(inFile$datapath))
    #        get((isfar)[8])
    #    })
    #    
    #    trt.effect2  <- reactive({
    #        
    #        inFile <- input$file1
    #        if (is.null(inFile))
    #            return(NULL)
    #        isfar <- (load(inFile$datapath))
    #        get((isfar)[9])
    #    })
    #    trt.effect3  <- reactive({
    #        
    #        inFile <- input$file1
    #        if (is.null(inFile))
    #            return(NULL)
    #        isfar <- (load(inFile$datapath))
    #        get((isfar)[10])
    #    })
    #    
    #    trt.effect4  <- reactive({
    #        
    #        inFile <- input$file1
    #        if (is.null(inFile))
    #            return(NULL)
    #        isfar <- (load(inFile$datapath))
    #        get((isfar)[11])
    #    })
    #    
    #    trt.effect5  <- reactive({
    #        
    #        inFile <- input$file1
    #        if (is.null(inFile))
    #            return(NULL)
    #        isfar <- (load(inFile$datapath))
    #        get((isfar)[4])
    #    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # collect simulation trt effect estimates from simulation and plot!
    # output$reg.plotL   <- renderPlot({         #means
    #     
    #     
    #     res <- trt.effect1()
    #     
    #     res2 <- trt.effect2()
    #     
    #     res3 <- trt.effect3()
    #     
    #     theta1 <- trt.effect4()
    #     
    #     sample <- random.sample()
    #     
    #     
    #     d1 <-  density( res[,1]) 
    #     d2 <-  density(res[,3] )
    #     d3 <-  density(res[,5] )
    #     d4 <-  density(res[,7] )
    #     d5 <-  density(res[,9] )
    #     d6 <-  density(res[,11] )
    #     
    #     d7 <-  density(res2[,1] )
    #     d8 <-  density(res2[,3] )
    #     
    #     d9 <-   density(res3[,1] )
    #     d10 <-  density(res3[,3] )
    #     d11 <-  density(res3[,5] )
    #     d12 <-  density(res3[,7] )
    #     
    #     dz <- max(c(d1$y, d2$y, d3$y, d4$y, d5$y, d6$y, d7$y, d8$y  , d9$y, d10$y, d11$y, d12$y  ))
    #     dx <- range(c(d1$x,d2$x,  d3$x, d4$x, d5$x, d6$x, d7$x, d8$x   , d9$x, d10$x, d11$x, d12$x  ))
    #     
    #     if (input$dist %in% "All") {
    #         
    #         plot((d1), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww,
    #              xlab="Treatment effect log odds",  
    #              ylab="Density")                           
    #         lines( (d2), col = "black", lty=w, lwd=ww)  
    #         lines( (d3), col = "red", lty=wz, lwd=ww)    
    #         lines( (d4), col = "red", lty=w, lwd=ww)          
    #         lines( (d5), col = "blue", lty=wz, lwd=ww)       
    #         lines( (d6), col = "blue", lty=w, lwd=ww)       
    #         lines( (d7), col = "purple", lty=wz, lwd=ww)       
    #         lines( (d8), col = "purple", lty=w, lwd=ww)       
    #         lines( (d9), col = "green", lty=wz, lwd=ww)       
    #         lines( (d10), col = "green", lty=w, lwd=ww)       
    #         lines( (d11), col = "grey", lty=wz, lwd=ww)       
    #         lines( (d12), col = "grey", lty=w, lwd=ww)  
    #         
    #     }
    #     
    #     else if (input$dist %in% "d1") {  #remove
    #         
    #         
    #         plot((d1), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww,
    #              xlab="Treatment effect log odds",  
    #              ylab="Density")  
    #         lines( (d2), col = "black", lty=w, lwd=ww)  
    #         
    #     }
    #     
    #     else if (input$dist %in% "d3") {  #remove
    #         
    #         
    #         plot((d3), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww,col="red",
    #              xlab="Treatment effect log odds",  
    #              ylab="Density")               
    #         lines( (d4), col = "red", lty=w, lwd=ww)          
    #         
    #     }
    #     
    #     else if (input$dist %in% "d5") {
    #         
    #         plot((d5), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww, col="blue",
    #              xlab="Treatment effect log odds",  
    #              ylab="Density")                    
    #         
    #         lines( (d6), col = "blue", lty=w, lwd=ww)       
    #         
    #     }
    #     
    #     else if (input$dist %in% "d7") {
    #         
    #         plot((d7), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww, col="purple",
    #              xlab="Treatment effect log odds",  
    #              ylab="Density") 
    #         
    #         lines( (d8), col = "purple", lty=w, lwd=ww)     
    #         
    #     }
    #     else if (input$dist %in% "d9") {
    #         
    #         plot((d9), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww, col="green",
    #              xlab="Treatment effect log odds",  
    #              ylab="Density")  
    #         
    #         lines( (d10), col = "green", lty=w, lwd=ww)     
    #         
    #     }
    #     
    #     else if (input$dist %in% "d11") {
    #         
    #         plot((d11), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww, col="grey",
    #              xlab="Treatment effect log odds",  
    #              ylab="Density")  
    #         
    #         lines( (d12), col = "grey", lty=w, lwd=ww)     
    #         
    #     }
    #     
    #     abline(v = theta1, col = "darkgrey")                
    #     legend("topright",       # Add legend to density
    #            legend = c(" adj. for true prognostic covariates", 
    #                       " not adj. for true prognostic covariates" ,
    #                       " adj. for covariates unrelated to outcome", 
    #                       " not adj. for covariates unrelated to outcome",
    #                       " adj. for mix of prognostic and unrelated to outcome", 
    #                       " not adj. mix of prognostic and unrelated to outcome", 
    #                       " adj. for correlated prognostic covariates", 
    #                       " not adj. for correlated prognostic covariates",
    #                       " adj. for imbalanced prognostic covariates", 
    #                       " not adj. for imbalanced prognostic covariates", 
    #                       " adj. for imbalanced covariates unrelated to outcome", 
    #                       " not adj. imbalanced covariates unrelated to outcome"
    #                       
    #            ),
    #            col = c("black", "black","red","red","blue", "blue", "purple", "purple", "green", "green", "grey", "grey"),
    #            lty = c(wz, w,wz,w,wz,w,wz,w,wz,w,wz,w)  ,lwd=ww
    #            , bty = "n", cex=1)
    # })
    # 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # output$reg.plotM <- renderPlot({         #standard errors
    #     
    #     # Get the  data
    #     
    #     res <- trt.effect1()
    #     
    #     res2 <- trt.effect2()
    #     
    #     res3 <- trt.effect3()
    #     
    #     se. <- trt.effect5()
    #     
    #     sample <- random.sample()
    #     
    #     d1 <-  density(res[,2] )
    #     d2 <-  density(res[,4] )
    #     d3 <-  density(res[,6] )
    #     d4 <-  density(res[,8] )
    #     d5 <-  density(res[,10] )
    #     d6 <-  density(res[,12] )
    #     d7 <-  density(res2[,2] )
    #     d8 <-  density(res2[,4] )
    #     
    #     d9 <-   density(res3[,2] )
    #     d10 <-  density(res3[,4] )
    #     d11 <-  density(res3[,6] )
    #     d12 <-  density(res3[,8] )
    #     
    #     
    #     
    #     dz <- max(c(d1$y, d2$y, d3$y, d4$y, d5$y, d6$y, d7$y, d8$y  , d9$y, d10$y, d11$y, d12$y    ))
    #     dx <- range(c(d1$x,d2$x,  d3$x, d4$x, d5$x, d6$x, d7$x, d8$x   , d9$x, d10$x, d11$x, d12$x    ))
    #     
    #     if (input$dist %in% "All") {
    #         
    #         plot( (d1), xlim = c(dx), main=paste0("Density of treatment standard error estimates, truth= ",p4(se.),""), ylim=c(0,dz),lty=wz, lwd=ww,
    #               xlab="Standard error of log odds trt effect",  
    #               ylab="Density")  
    #         lines( (d2), col = "black", lty=w, lwd=ww)  
    #         lines( (d3), col = "red", lty=wz, lwd=ww)    
    #         lines( (d4), col = "red", lty=w, lwd=ww)          
    #         lines( (d5), col = "blue", lty=wz, lwd=ww)       
    #         lines( (d6), col = "blue", lty=w, lwd=ww)       
    #         lines( (d7), col = "purple", lty=wz, lwd=ww)       
    #         lines( (d8), col = "purple", lty=w, lwd=ww)       
    #         lines( (d9), col = "green", lty=wz, lwd=ww)       
    #         lines( (d10), col = "green", lty=w, lwd=ww)       
    #         lines( (d11), col = "grey", lty=wz, lwd=ww)       
    #         lines( (d12), col = "grey", lty=w, lwd=ww)  
    #         
    #     }
    #     
    #     
    #     else if (input$dist %in% "d1") {  
    #         
    #         plot((d1), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww,
    #              xlab="Standard error of log odds trt effect",  
    #              ylab="Density") 
    #         lines( (d2), col = "black", lty=w, lwd=ww)  
    #         
    #     }
    #     
    #     else if (input$dist %in% "d3") {  
    #         
    #         
    #         plot((d3), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww,col="red",
    #              xlab="Standard error of log odds trt effect",  
    #              ylab="Density")  
    #         lines( (d4), col = "red", lty=w, lwd=ww)          
    #         
    #     }
    #     
    #     else if (input$dist %in% "d5") {
    #         
    #         plot((d5), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww, col="blue",
    #              xlab="Standard error of log odds trt effect",  
    #              ylab="Density")  
    #         
    #         lines( (d6), col = "blue", lty=w, lwd=ww)       
    #         
    #     }
    #     
    #     else if (input$dist %in% "d7") {
    #         
    #         plot((d7), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww, col="purple",
    #              xlab="Standard error of log odds trt effect",  
    #              ylab="Density") 
    #         
    #         lines( (d8), col = "purple", lty=w, lwd=ww)     
    #         
    #     }
    #     
    #     else if (input$dist %in% "d9") {
    #         
    #         plot((d9), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww, col="green",
    #              xlab="Standard error of log odds trt effect",  
    #              ylab="Density")  
    #         
    #         lines( (d10), col = "green", lty=w, lwd=ww)     
    #         
    #     }
    #     
    #     else if (input$dist %in% "d11") {
    #         
    #         plot((d11), xlim = dx, main=paste0("Density of treatment standard error estimates, truth= ",p3(se.),""), ylim=c(0,dz),lty=wz, lwd=ww, col="grey",
    #              xlab="Standard error of log odds trt effect",  
    #              ylab="Density")  
    #         
    #         lines( (d12), col = "grey", lty=w, lwd=ww)     
    #         
    #     }
    #     
    #     abline(v = se., col = "darkgrey")   
    #     legend("topright",           # Add legend to density
    #            legend = c(" adj. for true prognostic covariates", 
    #                       " not adj. for true prognostic covariates" ,
    #                       " adj. for covariates unrelated to outcome", 
    #                       " not adj. for covariates unrelated to outcome",
    #                       " adj. for mix of prognostic and unrelated to outcome", 
    #                       " not adj. mix of prognostic and unrelated to outcome", 
    #                       " adj. for correlated prognostic covariates", 
    #                       " not adj. for correlated prognostic covariates",
    #                       " adj. for imbalanced prognostic covariates", 
    #                       " not adj. for imbalanced prognostic covariates", 
    #                       " adj. for imbalanced covariates unrelated to outcome", 
    #                       " not adj. imbalanced covariates unrelated to outcome"
    #                       
    #            ),
    #            col = c("black", "black","red","red","blue", "blue", "purple", "purple", "green", "green", "grey", "grey"),
    #            lty = c(wz, w,wz,w,wz,w,wz,w,wz,w,wz,w) ,lwd=ww
    #            , bty = "n", cex=1)
    # })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
        
    })
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This is where a new sample is instigated and inputs converted to numeric
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    random.sample2 <- reactive({
        
        foo <- input$resample2   # button is labelled this
        
        NN <- as.numeric(unlist(strsplit(input$NN,",")))
 
        pp1 <- as.numeric(unlist(strsplit(input$pp1,",")))
        
        pp2 <- as.numeric(unlist(strsplit(input$pp2,",")))
        
        or <- as.numeric(unlist(strsplit(input$or,",")))    
        
        allocation <- (as.numeric(unlist(strsplit(input$allocation,","))))   

        return(list(  
            NN=NN,  
            pp1=pp1,
            pp2=pp2,
            or= or, 
            allocation=allocation
            
        ))
        
    })
    
    # function to simulate logistic regression with one categorical covariate 'drug'
    
    twobytwo <- eventReactive(input$sim,{
        
        sample2 <- random.sample2()
       
        NN=sample2$NN
        pp1=sample2$pp1
        pp2=sample2$pp2
        or= sample2$or 
        allocation=sample2$allocation
        
        odds<- pp1 / (1-pp1)
        
            fun.d<-function(nsample, drug.allocation, 
                            alpha,  beta.drug,
                            seed=NULL){ 
                
                if (!is.null(seed)) set.seed(seed)
                
                drug<- (rbinom(nsample, 1, prob =drug.allocation ))   
                
                Xmat <- model.matrix(~ drug )
                beta.vec <- c(alpha,  beta.drug )
                
                lin.pred <- Xmat[,] %*% beta.vec                 # Value of lin.predictor
                exp.p <- exp(lin.pred) / (1 + exp(lin.pred))     # Expected proportion
                y <- rbinom(n = nsample, size = 1, prob = exp.p) # Add binomial noise
                #y<- runif(nsample) <  exp.p                     # alternatively ads noise in this way
                
                d<-as.data.frame(cbind(y, drug))         # create a dataset
                
                return(d)
            }
            
            # functions to fit model and pull out a stat
            # lrtest
            simfunc <- function(d) {
                fit1 <- glm(y  ~ drug , d, family = binomial) 
                fit2 <- glm(y  ~ 1, d, family = binomial) 
                c(anova(fit1,fit2, test='Chisq')[2,5] )
            }
            
            # wald test
            simfunc <- function(d) {
                fit1 <- glm(y  ~ drug , d, family = binomial) 
                c( summary(fit1)$coef["drug","Pr(>|z|)"] )
            }
    
 
    
    out <- replicate(100, simfunc(fun.d( nsample=NN,   #NN 
                                          drug.allocation=allocation,  
                                          alpha=log(odds),  
                                          beta.drug=log(or))))
    
    
    
    pow <-  mean( out < 0.05, na.rm=TRUE )   
    
    
    # frank harrell
    
    FH <- Hmisc::bpower(p1=pp1, odds.ratio=c(or ), n=NN, alpha=c(0.05))
    
   
   return(list(  
       
       pow=pow, FH=FH
       
       )) 
   
   
    })
    
    
  #  twobytwo1 <- eventReactive(input$sim,{
        
    
 
    output$obs <- renderTable( {
        
       sample2 <- random.sample2()
        
        NN=sample2$NN
        pp1=sample2$pp1
        pp2=sample2$pp2
        or= sample2$or 
        allocation=sample2$allocation
        # one dataset
        d <- fun.d( nsample=NN, drug.allocation=allocation,  
                    alpha=log(pp1/(1-pp1)),  
                    beta.drug=log(or))
        
        
        
        df <- (table(d))      
        
        N_metrics <- matrix(c(df[1,1], df[2,1], df[1,2], df[2,2]), ncol = 2)
        
        colnames(N_metrics) <- c("placebo", "drug")
        row.names(N_metrics) <- c ("no response", "response")
        
        N <- addmargins(N_metrics)
        
        N <- round(N,0)
        
    }, rownames = TRUE, digits=0)
    
  #  })
 
    
   output$pow1 <- renderPrint({
       
       return(twobytwo()$pow)
       
   }) 
   
   
   output$pow2 <- renderPrint({
       
       return(twobytwo()$FH)
       
   }) 
   
   
   
   
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# loading in user data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


# Run the application 
shinyApp(ui = ui, server = server)