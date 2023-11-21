library(shiny)
library(luna)
library(shinybusy)
library(shinyWidgets)
library(shinyjs)
library(dplyr)
library(datamods)
library(lubridate)
library(shinydashboard)

ui <- navbarPage(
  useShinyjs(),
  id = "tabset",
  title = "Hypnoscope",

  # Component 1
  tabPanel(
    title = "N = 1",
    fluidRow(
      add_busy_spinner(spin = "fading-circle"),
      column( # Sidebar
        2,
        # style = "background-color:#f2f0eb; border-top: 100px; height: 800px; width: 250px;",
        # Input: Select a file ----
        fileInput("upload1",
          label = NULL,
          multiple = TRUE,
          accept = c(
            ".annot",
            ".eannot",
            ".xml",
            ".hypnos"
          )
        ),
        textOutput("text.header1a"),
        hr(style = "border-color: #d9d9d9"),
        actionButton("load.default", "Example hypnogram",width="100%"),
        hr(style = "border-color: #d9d9d9"),
        import_copypaste_ui(id = "myid", title = NULL, name_field = F)
      ),
      column( # Main
        10,
        plotOutput("hypno1",
          width = "100%", height = "180px",
        ),
        div(style = "margin-top: 50px"),
        tabsetPanel(
          id = "maintabs",
          tabPanel("Summaries", DT::dataTableOutput("table.hypno", width = "100%")),
          tabPanel("Stages", DT::dataTableOutput("table.hypno.stages")),
          tabPanel("Times", DT::dataTableOutput("table.hypno.times")),
          tabPanel("Cycles", DT::dataTableOutput("table.hypno.cycles")),
          tabPanel("Epochs", DT::dataTableOutput("table.hypno.epochs")),
        ),
      )
    )
  ),

  # Component 2
  tabPanel(
    title = "N > 1",
    fluidRow(
      add_busy_spinner(spin = "fading-circle"),
      column(
        2,
        # Input: Select a file ----
        fileInput("upload2",
          label = NULL,
          multiple = TRUE,
          accept = c(
            ".csv",
            ".tsv",
	    ".txt",
            ".hypnos",".gz"
          )
        ),
        textOutput("text.header2a"),
        hr(style = "border-color: #d9d9d9"),
        actionButton("load.default2", "Example hypnograms",width="100%"),
        hr(style = "border-color: #d9d9d9"),
        fileInput("upload3",
          label = "Optional covariates",
          multiple = TRUE,
          accept = c( ".csv",".tsv",".txt") ) ,
        hr(style = "border-color: #d9d9d9"),
        selectInput(inputId = "ultradian2", label = "Align by", choices = c("Clock-time", "Elapsed-time"), selected = "Clock-time", multiple = F),
        selectInput(inputId = "sort", label = "Sort by", choices = c("Unsorted","Sleep Onset", "Start of recording"), selected = "Unsorted", multiple = F),
        selectInput(inputId = "color", label = "Color scheme",
	 choices = c("All", "W", "WASO", "N1", "N2", "N3", "NR", "R", "S", "Cycle","L", "U"), selected = "All", multiple = F),
        hr(style = "border-color: #d9d9d9"),
        uiOutput(outputId = "n"),
        hr(style = "border-color: #d9d9d9"),
        selectInput(inputId = "procn1", label = "Process individual", choices = NULL , selected = NULL,multiple=F),
	actionButton( "runn1", "Generate hypnogram stats",width="100%"),
      ),

      column(
        10,
        align = "left",
	plotOutput("hyp1",  width = "1280px", height = "50px") , 
         box(
            style = "width:1280px; overflow-x: scroll; height:800px; overflow-y: scroll;",
            imageOutput("myImage", hover = hoverOpts(id="hover1" , delay=150 )) 
             )
	    
          )
    )
  ),

  # Component 3
  tabPanel(
    title = "Covariates",
    DT::dataTableOutput("covar.table") 
  ),	

  # Component 4
  tabPanel(
    title = "Help",
    tags$head(
      tags$style(HTML("
            code {
                display:block;
                padding:9.5px;
                margin:0 0 10px;
                margin-top:10px;
                font-size:13px;
                line-height:20px;
                word-break:break-all;
                word-wrap:break-word;
                white-space:pre-wrap;
                background-color:#F5F5F5;
                border:1px solid rgba(0,0,0,0.15);
                border-radius:4px;
                font-family:monospace;
            }"))
    ),
    h3("N=1 mode: view a single hypnogram and generate hypnogram statistics"),
    p("Upload a stage annotation file in any of the following formats (see Luna website for details):"),
    strong("   .annot"),
    p("Luna standard annotation format"),
    strong("   .eannot"),
    p("Epoch-based one stage code per row"),
    strong("   .xml"),
    p("NSRR-style XML annotation file format"),
    p(
      "Alternatively, in the Copy & Paste box, you can paste the staging annotations",
      span("(W,N1,N2,N3,R,L or ?)", style = "color:blue"),
      ". Format: one row per epoch, i.e each row contains a single label for that epoch, similar to an .eannot file"
    ),
    em("By default, Luna assumes epochs are 30 seconds in duration."),
    hr(),
    h3("N>1 mode: viewing and aligning multiple hypnograms"),
    p("Upload a '.hypnos' file in the following format: 4 tab-delimited columns with a header:"),
    p("  - first column: subject identifier (must be labelled ID)"),
    p("  - second column: epoch number (can be labelled anything in the header)"),
    p("  - third column: clock-time in hh:mm:ss 24-hour format (can be labelled anything in the header)"),
    p("  - fourth column: stage code as N1,N2,N3,R,W,L or ? (can be labelled anything in the header)"),
    br(),
    p("To use Luna to create the required inputs (given a sample list and records with stage annotations):"),
    code("luna s.lst -o out.db -s HYPNO epoch verbose=F"),
    code("destrat out.db +HYPNO -r E -v OSTAGE CLOCK_TIME > multiple.hypnos"),
    p("Note that using OSTAGE instead of STAGE takes the original stage, i.e. before any edits imposed by Luna" ),
    p("You can optionally include the sleep cycle code (1,2,3 or NA is not in a cycle) in the 5th column:"),
    code("destrat out.db +HYPNO -r E -v OSTAGE CLOCK_TIME CYCLE > multiple.hypnos"),
    p("You can also work with compressed files directly:"),
    code("destrat out.db +HYPNO -r E -v OSTAGE CLOCK_TIME | gzip > multiple.hypnos.gz"),
    p("The file extension should be either .hypno or .hypnos.gz"),
    hr(),
    h5("Attaching other covariate data"),
    p("You can also attached artitrary covariate data by uploading a tab-delimited text file in Luna vars format (first column must be ID)")
  )
)
