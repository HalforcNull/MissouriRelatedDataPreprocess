library(shiny)


shinyUI(  
  fluidPage(
    # h1("Upload File"),
    # fileInput("fCount", label = "Count Table"), 
    # fileInput("fAttr", label = "Attr Table"), 
    # conditionalPanel(condition="true",
    #                  fluidRow(
    #                     column( 1, h2("Samples"),
    #                             dragSetUI("div1", textval = c("a","b"))
    #                     ),
    #                     column( 2,
    #                             textInput("textCtrl1", label=NULL, value = "Control 1"),
    #                             dropUI("ctrl1"),
    #                             textInput("textCtrl2", label=NULL, value = "Control 2"),
    #                             dropUI("ctrl2")
    #                     )
    #                  )
    # )
    # 
    
    
      column(3, h2("Step 1"),
             p("Do you have read count data ?"),
             actionButton("userHaveReadCount", "Yes"),
             actionButton("userDoNotHaveReadCount", "No"),
             textOutput("MessageStep1")
             ),
      
      # Panel 2
      conditionalPanel(condition="output.showPanel2 == true",
                       column(3, h2("Go To Idep") ))
      
      
    
    
  )
)
