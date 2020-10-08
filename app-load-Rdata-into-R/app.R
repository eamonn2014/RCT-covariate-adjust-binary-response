#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

rm(list=ls())
pp<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/A%205000%20default%20settings%20theta%20log1.5%20-1.00%20-0.67%20-0.43.Rdata"
#load(url(pp))


library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("Uploading Files"),
    sidebarLayout(
        # sidebarPanel(
        #     fileInput('file1', 'Choose file to upload',
        #               accept = c(
        #                   '.Rdata'
        #               )
        #     ),
        #     tags$hr(),
        #     
        #     tags$hr(),
        #     
        # ),
        # mainPanel(
        #     tableOutput('contents'),
        #     tableOutput('contents2'),
        #     tableOutput('contents3'),
        #     #div( verbatimTextOutput("x") )
        # )
        # 
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        sidebarPanel(
           # textInput("path", "File:"),
           # actionButton("browse", "Browse"),
            tags$br(),
            actionButton("upload", "Upload Data")
        ),
        
        mainPanel(
            verbatimTextOutput('content')
        )
        
        
        
        
        
        
        
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        
        
        
        
        
    ))



# Define server logic required to draw a histogram
options(shiny.maxRequestSize = 9*1024^2)

server <- function(input, output) {
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    contentInput <- reactive({ 
        
       if(input$upload == 0) return()
        
       isolate({
           isfar <-  load(url(pp))
                get((isfar)[12])
      })
    })
    
    output$content <- renderPrint({
        contentInput()
    })
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    
    
    
    
    
    # 
    # 
    # 
    # 
    # output$contents <- renderTable({
    #     # input$file1 will be NULL initially. After the user selects
    #     # and uploads a file, it will be a data frame with 'name',
    #     # 'size', 'type', and 'datapath' columns. The 'datapath'
    #     # column will contain the local filenames where the data can
    #     # be found.
    #     
    #     inFile <- input$file1
    #     
    #     
    # })
    # 
    # output$contents2 <- renderTable({
    #     
    #     inFile <- input$file1
    #     if (is.null(inFile))
    #         return(NULL)
    #     load(inFile$datapath)
    #     
    # })
    # 
    # 
    # output$contents3 <- renderTable({
    #     
    #     inFile <- input$file1
    #     if (is.null(inFile))
    #         return(NULL)
    #     isfar <- (load(inFile$datapath))
    #     get((isfar)[12])
    # })
    # 
    
    
    
    # 
    # base <- reactive ({
    #     
    #     inFile <- input$file1
    #     
    #     x <- get(inFile$w) 
    #     
    #     return(x)
    #     
    #     })
    
    
    
    
    
    # base <- reactive ({(get(load(input$file1$datapath)))})
    # 
    # zz <-renderPrint({
    # mydata <- base()$zz
    # })
    
    
    
    #  zz <-renderPrint({
    #     
    #      inFile <- input$file1
    #      
    #      if (is.null(inFile))
    #          
    #          return(NULL)
    #      
    #      isfar <- (load(inFile$datapath))
    #      d <- print((isfar)[12])
    #      print(contents()$d)
    #      
    # })
    
    
    
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)