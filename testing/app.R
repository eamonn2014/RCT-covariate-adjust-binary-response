#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# stripped down code to improve plot presentation on loading
# before any button is pressed we have ugly looking tab
# because plot code is run before button is pressed
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

options(max.print=1000000)    

fig.width8 <- 1380
fig.height7 <- 500
 

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

 
RR=.37  # used to limit correlations between variables
sd1= 3  # for X covariates sd

# links to Rdata objects uploaded to Git, these are pre run simulations see the save function
# see line 1224 for the save function
pp<- "https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/A%205000%20default%20settings%20theta%20log1.5%20-1.00%20-0.67%20-0.43.Rdata" # 5000 default log1.5 -1 -.67 -.43
pp2<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/B%205000%20default%20settings%20theta%20log0.5%20-1.68%20-1.39%20%200.71.Rdata"
pp3<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/C%205000%20default%20settings%20theta%20log2%20-3.46%20-1.05%20%201.15%20p1=.75.Rdata"
pp4<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/D_10Ksims_5covariates_p1_0.12_theta_log1.3_covariates_-1.02%20_0.42_0.43_0.61%20_1.01_3_prog.Rdata"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ui <-  fluidPage(theme = shinytheme("journal"), #https://www.rdocumentation.org/packages/shinythemes/versions/1.1.2 , paper another option to try
                 # paper
                 useShinyalert(),  # Set up shinyalert
                 setBackgroundColor(
                     color = c( "#2171B5", "#F7FBFF"), 
                     gradient = "linear",
                     direction = "bottom"
                 ),
                 
              
                 h3("  "), 
                 sidebarLayout(
                     
                     sidebarPanel( width=3 ,
                                   
                                   tags$style(type="text/css", ".span8 .well { background-color: #00FFFF; }"),
                                   
                                   actionButton(inputId='ab1', label="R Shiny ",   icon = icon("th"),   
                                                onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/RCT-covariate-adjust-binary-response/master/cov-adj-binary-response/app.R', '_blank')"), 
                                   actionButton(inputId='ab1', label="R code",   icon = icon("th"),   
                                                onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/RCT-covariate-adjust-binary-response/master/cov-adj-binary-response/R%20code.R', '_blank')"),  
                                #   actionButton("resample", "Simulate a new sample"),
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
                                       
                                       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
 
                                       
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
          
                                   
                                   # here have action buttons clicking will load a pre run simulation that can be examined
                                   tabPanel( "Load in pre-run simulations",
                                             
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
                                             
                                             ###HERE WE OUTPUT
                                             # this spinner indicating something is loading does not seem to work
                                             # shinycssloaders::withSpinner(
                                             #     div(plotOutput("reg.plotLL",  width=fig.width8, height=fig.height7)),  #trt est plot
                                             # ) ,
                                             # 
                                             withSpinner(plotOutput("reg.plotLL",  width=fig.width8, height=fig.height7),6),
                                             
                                             # this spinner indicating something is loading does not seem to work
                                             withSpinner(plotOutput("reg.plotMM",  width=fig.width8, height=fig.height7),6),
                                             # this spinner indicating something is loading does work
                                             shinycssloaders::withSpinner( 
                                                 verbatimTextOutput('content1'),  #summary table
                                                 
                                             ),
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
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    values <- reactiveValues(pressed=0)
    
    observeEvent(input$upload, {    
        values$pressed <- 1
        })
    
    observeEvent(input$upload2, {    
        values$pressed <- 1
    })
    
    observeEvent(input$upload3, {    
        values$pressed <- 1
    })
    
    observeEvent(input$upload4, {    
        values$pressed <- 1
    })
    
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # collect simulation trt effect estimates from upload dnd plot!
    # New function to plot basically copy of earlier function
    
    
  
    
    
    
   output$reg.plotLL   <- renderPlot({         #means
            
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
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
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
    
 
        
        output$reg.plotMM <- renderPlot({         #standard errors
            
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
        

 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    
})
 

# Run the application 
shinyApp(ui = ui, server = server)