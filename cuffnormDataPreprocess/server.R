library(shiny)

#upload file max size = 300 MB
options(shiny.maxRequestSize=300*1024^2) 



#
shinyServer(function(input, output) {

  # observeEvent(input$fCount,{
  #   if(! "SampleList" %in% colnames(input))
  #   {
  #     req(input$fCount)
  #     tryCatch(
  #       { 
  #         df <- read.table(input$fCount$datapath)
  #         output$SampleList <- renderTable( head(df) )
  #       },
  #       error = function(e) {
  #         # return a safeError if a parsing error occurs
  #         stop(safeError(e))
  #       }
  #     )
  #     
  #     output$stage <- renderText( 'group samples' )
  #     outputOptions(output, "stage", suspendWhenHidden = FALSE)
  #     
  #     
  #   }
  # })
  # 
  # observeEvent(input$div2,{
  #   output$foo = renderText(
  #     paste("The dropUI element currently contains:", input$div2))
  # })
  
  STATE.UserHaveReadCount = FALSE
  observeEvent(input$userHaveReadCount, {
    STATE.UserHaveReadCount = TRUE
    output$MessageStep1 <- renderText('Great! Do you need rename your columns?')
    output$showPanel2<- reactive(STATE.UserHaveReadCount)
    outputOptions(output, "showPanel2", suspendWhenHidden = FALSE)
  })
  
  observeEvent(input$userDoNotHaveReadCount, {
    STATE.UserHaveReadCount = FALSE
    output$MessageStep1 <- renderText('OK. We need get the read count first.')
    output$showPanel2<- reactive(STATE.UserHaveReadCount)
    outputOptions(output, "showPanel2", suspendWhenHidden = FALSE)
  })
  
})
