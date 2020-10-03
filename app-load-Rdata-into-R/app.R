#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("Uploading Files"),
    sidebarLayout(
        sidebarPanel(
            fileInput('file1', 'Choose file to upload',
                      accept = c(
                          '.Rdata'
                      )
            ),
            tags$hr(),
            
            tags$hr(),
            
        ),
        mainPanel(
            tableOutput('contents'),
            tableOutput('contents2'),
            tableOutput('contents3'),
            #div( verbatimTextOutput("x") )
        )
        
        
        
        
    ))



# Define server logic required to draw a histogram
options(shiny.maxRequestSize = 9*1024^2)

server <- function(input, output) {
    
    
    output$contents <- renderTable({
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, it will be a data frame with 'name',
        # 'size', 'type', and 'datapath' columns. The 'datapath'
        # column will contain the local filenames where the data can
        # be found.
        
        inFile <- input$file1
        
        
    })
    
    output$contents2 <- renderTable({
        
        inFile <- input$file1
        if (is.null(inFile))
            return(NULL)
        load(inFile$datapath)
        
    })
    
    
    output$contents3 <- renderTable({
        
        inFile <- input$file1
        if (is.null(inFile))
            return(NULL)
        isfar <- (load(inFile$datapath))
        get((isfar)[12])
    })
    
    
    
    
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