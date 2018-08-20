library(shiny)
library(shinyDND)


shinyUI(  
  fluidPage(
    h1("Upload File"),
    fileInput("fCount", label = "Count Table"), 
    fileInput("fAttr", label = "Attr Table"), 
    conditionalPanel(condition="true",
                     fluidRow(
                        column( 1, h2("Samples"),
                                dragSetUI("div1", textval = c("a","b"))
                        ),
                        column( 2,
                                textInput("textCtrl1", label=NULL, value = "Control 1"),
                                dropUI("ctrl1"),
                                textInput("textCtrl2", label=NULL, value = "Control 2"),
                                dropUI("ctrl2")
                        )
                     )
    )
    
  )
)
