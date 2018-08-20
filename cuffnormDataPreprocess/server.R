library(shiny)
library(shinyDND)

#upload file max size = 300 MB
options(shiny.maxRequestSize=300*1024^2) 



#
shinyServer(function(input, output) {

  observeEvent(input$fCount,{
    if(! "SampleList" %in% colnames(input))
    {
      req(input$fCount)
      tryCatch(
        { 
          df <- read.table(input$fCount$datapath)
          output$SampleList <- renderTable( head(df) )
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
      
      output$stage <- renderText( 'group samples' )
      outputOptions(output, "stage", suspendWhenHidden = FALSE)
      
      
    }
  })
  
  observeEvent(input$div2,{
    output$foo = renderText(
      paste("The dropUI element currently contains:", input$div2))
  })
})
