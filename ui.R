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
        actionButton("load.default", "Example"),
        hr(style = "border-color: #d9d9d9"),
        import_copypaste_ui(id = "myid", title = NULL, name_field = F)
      ),
      column( # Main
        10,
        plotOutput("hypno1",
          width = "100%", height = "200px",
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
            ".hypnos"
          )
        ),
        textOutput("text.header2a"),
        hr(style = "border-color: #d9d9d9"),
        actionButton("load.default2", "Example"),
        hr(style = "border-color: #d9d9d9"),
        selectInput(inputId = "ultradian2", label = "Align by", choices = c("CLOCK_TIME", "ONSET"), selected = "CLOCK_TIME", multiple = F),
        selectInput(inputId = "sort", label = "Sort by", choices = c("Default"), selected = "Default", multiple = F),
        selectInput(inputId = "color", label = "Color scheme", choices = c("All", "W", "N1", "N2", "N3", "R","L","U"), selected = "All", multiple = F),
        hr(style = "border-color: #d9d9d9"),
        uiOutput(outputId = "n")
      ),
      column(
        10,
        align = "center",
        box(
          style = "width:1200px; height:800px; overflow-y: scroll;",
          imageOutput("myImage")
        )
      )
    )
  ),

  # Component 3
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
    p("1) For viewing a single hypnogram, Please upload a file in any of the following formats."),
    strong(" .annot"),
    br(),
    strong(" .eannot"),
    br(),
    strong(" .xml"),
    br(),
    br(),
    p(
      "In copy & paste box, you can paste the staging annotations",
      span("(W,N1,N2,N3,R)", style = "color:blue"),
      ". The format is one row per epoch, i.e each row contains a single label, that is attached to that epoch."
    ),
    em("By default, Luna assumes epochs are 30-seconds in duration."),
    # div(" In copy & paste box, you can paste staging annotations (W,N1,N2,N3,R).
    #  The format is one row per epoch, i.e each row contains a single label.", style = "color:blue"),
    # p(".eannot", style = "font-family: 'times'; font-si16pt"),
    # p(".xml", style = "font-family: 'times'; font-si16pt"),
    hr(),
    p("2) For viewing mutiple hypnograms, Please upload a file in the following format."),
    strong(" .hypnos"),
    br(),
    br(),
    p("To create this format, Please use the following Luna command."),
    code("luna s.lst -o out.db -s HYPNO epoch verbose=F"),
    code("destrat out.db +HYPNO -r E -v OSTAGE CLOCK_TIME > mutiple.hypnos"),
    br(),
  )
)
