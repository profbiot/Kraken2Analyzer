#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(seqinr)
library(stringr)

#Define UI functions for file input tools
TAB_UI <- function(id) {
    ns = NS(id)
    
    list(
        textOutput(ns("output_area")), 
        fileInput("file1", "Choose .tabular report from Galaxy version of Kraken2",
                  multiple = FALSE,
                  accept = ".tabular")
    )
}



# Define UI for application that draws a histogram
ui <- fluidPage(
    title="Shannon Diversity Calculator",
    # Include css style sheet
    includeCSS("styles.css"),
    
    # App title ----
    titlePanel(h1("Endicott College Bioinformatics: Microbiome Diversity Calculator")),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            # Input: Select a file ----
            TAB_UI("file1"),
            # Horizontal line, provides spacing
            tags$hr(),
            actionButton("pvButton2", "View Species Only Dataframe"),
            # Horizontal line, provides spacing
            tags$hr(),
            actionButton("calcButton", "Calculate Shannon Diversity"),
            # Horizontal line, provides spacing
            tags$hr(),
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            h4("Shannon Diversity Calculations"),
            div(style="width:500px;",verbatimTextOutput("value")),
            
            # Output: Preview Original File
            h4("Preview Data from file"),
            tableOutput("contents")
        )
    )
) 


# Define server logic required to draw a histogram
server <- function(input, output) {

    output$contents <- renderTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        
        req(input$file1)
        req(input$pvButton2)
        
        # when reading semicolon separated files,
        # having a comma separator causes `read.csv` to error
        tryCatch(
            {
                df <- read.table(input$file1$datapath,sep="\t", header=FALSE)
                names(df)<-c('Classification','Count')
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        #Create dataframe with just species rows using logical/grepl filter
        sRows <- df[grepl("\\|s", df[["Classification"]]),]
        
        #Use for loop with gsub to simplify Species name
        for (i in 1:nrow(sRows)){
            sRows[i,1]<-gsub("d__.*s__","",sRows[i,1],perl=TRUE)
        }
        
        #Format column names and write file
        names(sRows) <- c('Species','Count')
        
        return(sRows)
    })

    
    output$value <- renderText({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        
        req(input$file1)
        req(input$calcButton)
        
        # when reading semicolon separated files,
        # having a comma separator causes `read.csv` to error
        tryCatch(
            {
                df <- read.table(input$file1$datapath,sep="\t", header=FALSE)
                names(df)<-c('Classification','Count')
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        
        sRows <- df[grepl("\\|s", df[["Classification"]]),]
        
        myS <- nrow(df)
        myTotal <- sum(df$Count)
        myP <- df$Count/myTotal 
        
        myNewColumn = myP*log(myP)
        
        myNew2 = myP*log(myP)^2
        
        SDI <- -1*sum(myNewColumn)
        
        Sum2 <- sum(myNew2)
        
        Hmax <- -myTotal*(1/myTotal)*log(1/myTotal)
        
        SE <- SDI/log(myS)
        
        var_pt1 <- (Sum2-(SDI^2))/myTotal
        
        var_pt2 <- (myS-1)/(2*myTotal^2)
        
        SVar <- (var_pt1 + var_pt2)^0.5
        
        myString <- paste("Shannon Diversity Index:",format(SDI,digits = 2, nsmall = 3)," Evenness:",format(SE,digits = 2, nsmall = 3)," Variance:",format(SVar,digits = 2, nsmall = 3))
        
        return(myString)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
