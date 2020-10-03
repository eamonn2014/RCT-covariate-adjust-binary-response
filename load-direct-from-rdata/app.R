
#shiny app that loads rdata by one click just need to supply url
#https://stackoverflow.com/questions/21813773/r-shiny-uploading-a-file-on-button-click
rm(list=ls())
pp<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/A%205000%20default%20settings%20theta%20log1.5%20-1.00%20-0.67%20-0.43.Rdata"
#load(url(pp))

library(shiny)

# Define UI  
ui <- fluidPage(
    titlePanel("Uploading a file from rdata directly"),
    sidebarLayout(
        
        sidebarPanel(
         
            tags$br(),
            actionButton("upload", "Upload Data")
        ),
        
        mainPanel(
            verbatimTextOutput('content')
        )
        ))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define server logic  
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

}

# Run the application 
shinyApp(ui = ui, server = server)