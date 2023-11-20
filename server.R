# server.R
library(shiny)
library(luna)
library(shinybusy)
library(shinyWidgets)
library(shinyjs)
library(dplyr)
library(datamods)
library(lubridate)
library(shinydashboard)
source("./helpers.R")
options(shiny.maxRequestSize = 2000 * 1024^2)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  # Values, an object that is created outside a reactive expression is updated depending
  # on inputs but is not re-executed whenever an input changes.
  values <- reactiveValues(opt = list())
  values2 <- reactiveValues(opt = list())

  #-------------------------------- Component 1 -------------------------------#
  observeEvent(input$upload1, {
    # Set to null
    values$file.details <- NULL
    annot.names <- annot.paths <- NULL

    idx <- c(
      grep("\\<eannot\\>", ignore.case = T, input$upload1$name),
      grep("\\<annot\\>", ignore.case = T, input$upload1$name),
      grep("\\<xml\\>", ignore.case = T, input$upload1$name)
    )

    annot.names <- input$upload1$name[idx]
    annot.paths <- input$upload1$datapath[idx]

    # Input file details
    values$file.details <- list(
      annot.names = annot.names,
      annot.paths = annot.paths
    )
  })

  # --------- Copy and Paste module----------------

  imported <- import_copypaste_server("myid",
    btn_show_data = FALSE, fread_args = list(col.names = "Annots", header = F, blank.lines.skip = T)
  )

  toListen <- reactive({
    list(imported$status(), imported$data())
  })

  observeEvent(toListen(), {
    if (!is.null(imported$status()) && !is.null(imported$data())) {
      reset(id = "upload1")
      txtPath <- tempfile(fileext = ".eannot")
      for (stage in imported$data()) cat(stage, "\n", file = txtPath, sep = "\n", append = TRUE)

      annot.names <- basename(txtPath)
      annot.paths <- txtPath
      lset("epoch-check", "100000000")

      # Input file details
      values$file.details <- list(
        annot.names = annot.names,
        annot.paths = annot.paths
      )
    }
  })

  #--------------------------------------------

  # initiate view w/ an example data set
  observeEvent(input$load.default, {
    reset(id = "upload1")

    values$file.details <-
      list(
        annot.names = "learn-nsrr02.xml",
        annot.paths = "data/learn-nsrr02.xml"
      )
    # Reset Shiny UI components
  })

  # Load data
  observeEvent(values$file.details, {
    # pull file names
    annot.names <- values$file.details[["annot.names"]]
    annot.paths <- values$file.details[["annot.paths"]]

    # clear old data out
    try(lrefresh())
    try(ldrop())

    # clear all
    values$opt <- NULL
    values$elen <- 30

    updateTabsetPanel(inputId = "tabset", selected = "N = 1")

    # register that we have new data attached
    values$hasannots <- !is.null(annot.names)


    if (values$hasannots) {
      values$opt[["annotnames"]] <- annot.names
      values$opt[["annotpaths"]] <- annot.paths
    }

    # some brief console output
    cat(" has annotations?", values$hasannots, "\n")

    # If in .annot format, Get start time from file
    file_ext <- endsWith(annot.paths, ".annot")
    if (isTRUE(file_ext)) {
      startTime <- extract_hms(annot.paths)
      lempty.edf(start = startTime)
    } else {
      lempty.edf()
    }

    # read all EDF+ annotations as class-level
    # so that they show in the display
    lset("edf-annot-class-all", "T")


    # Add any annotations
    eannot_file <- endsWith(values$opt[["annotpaths"]], ".eannot")
    if (eannot_file) {
      lset("epoch-check", "100000000")
    }

    for (a in values$opt[["annotpaths"]]) {
      cat("attaching", a, "\n")
      ladd.annot.file(a)
    }

    # Restructure EDF to match the length of annotations
    ret <- leval("ANNOTS")
    end_sec <- sum(subset(ret$ANNOTS$ANNOT, ANNOT %in% c("N1", "N2", "N3", "W", "R", "L", "?"))$DUR)
    cmd <- paste("MASK sec=0-", end_sec, " & RE", sep = "")
    leval(cmd)

    # kick off initial analyses
    init()
  })

  init <- function() {
    values$elen <- 30
    ne <- lepoch()
    if (ne == 0) values$elen <- 1

    ret <- leval(paste("EPOCH dur=", values$elen, " verbose & SEGMENTS", sep = ""))
    values$opt[["init.epochs"]] <- ret$EPOCH$E
    values$opt[["ne"]] <- dim(ret$EPOCH$E)[1]
    values$opt[["init.secs"]] <- max(ret$EPOCH$E$STOP)
    values$opt[["init.segidx"]] <- ret$SEGMENTS$SEG[, c("START", "STOP")]
    session$resetBrush("hypno_brush")

    # Get stage-aligned epochs and hypnogram
    ret <- leval(paste("EPOCH align verbose dur=", values$elen, sep = ""))
    values$opt[["ne.aligned"]] <- dim(ret$EPOCH$E)[1]
    values$opt[["init.epochs.aligned"]] <- ret$EPOCH$E
    ret <- leval("HYPNO epoch")

    stgs <- leval("STAGE force")$STAGE$E$STAGE
    values$hasstaging <- !is.null(stgs)
    values$variable.staging <- F
    if (values$hasstaging) {
      stgs <- stgs[stgs == "N1" | stgs == "N2" | stgs == "N3" | stgs == "R" | stgs == "W"]
      nstgs10 <- sum(table(stgs) >= 10)
      values$variable.staging <- nstgs10 >= 2
    }
    if (values$hasstaging) {
      values$opt[["hypno.stats"]] <- ret$HYPNO$BL
      values$opt[["hypno.epochs"]] <- ret$HYPNO$E
      values$opt[["all.hypno.epochs"]] <- ret$HYPNO$E
      values$opt[["ss"]] <- ret$HYPNO$E[, c("E", "STAGE", "START_SEC")]
      values$opt[["ss"]]$STOP_SEC <- values$opt[["ss"]]$START_SEC + values$elen
      values$opt[["hypno.cycles"]] <- ret$HYPNO$C
      values$opt[["hypno.stages"]] <- ret$HYPNO$SS
    }

    # report to console
    cat(" # epochs (raw)", values$opt[["ne"]], "\n")
    cat(" # epochs (stage-aligned)", values$opt[["ne.aligned"]], "\n")
    cat(" has-staging?", values$hasstaging, "\n")

    output$hypno1 <- renderPlot({
      req(values$hasstaging)
      par(mar = c(5, 2.5, 0, 0))
      lhypno2(values$opt[["hypno.epochs"]],
        cycles = values$opt[["hypno.epochs"]]$CYCLE,
        times = values$opt[["init.epochs.aligned"]]$START
      )
    })

    output$text.header1a <- renderText({
      req(values$hasstaging)
      values$opt[["annotnames"]]
    })

    # Base-level EDF headers output
    output$table.header3 <- DT::renderDataTable({
      df <- values$opt[["header1"]]
      df$ID <- df$EDF_ID <- NULL
      df <- df[, c("EDF_TYPE", "NS", "START_DATE", "START_TIME", "STOP_TIME", "REC_DUR_HMS", "REC_DUR_SEC", "EPOCH", "TOT_DUR_HMS", "TOT_DUR_SEC", "NR", "REC_DUR")]
      df <- data.frame(t(df))
      df$VAR <- rownames(df)
      names(df) <- c("Value", "Variable")
      df <- df[, c("Variable", "Value")]
      df$Variable <- c(
        "EDF type", "Number of signals", "Start date", "Start time", "Stop time", "Duration (h:m:s)", "Duration (sec)", "Duration (epochs)",
        "Total duration (h:m:s)", "Total duration (sec)", "Number of records", "Record duration (sec)"
      )

      DT::datatable(df,
        extensions = c("Buttons"),
        options = list(
          scrollY = "375px",
          paging = F,
          info = FALSE,
          searching = FALSE,
          dom = "tB", buttons = list(list(extend = "copy", text = "Copy")),
          columnDefs = list(list(className = "dt-left", targets = "_all"))
        ),
        rownames = FALSE
      )
    })

    # Channel-wise EDF headers output
    output$table.header2 <- DT::renderDataTable({
      DT::datatable(values$opt[["header2"]],
        extensions = c("Buttons"),
        options = list(
          scrollY = "375px",
          paging = F,
          info = FALSE,
          searching = FALSE,
          dom = "tB", buttons = list(list(extend = "copy", text = "Copy")),
          columnDefs = list(list(className = "dt-center", targets = "_all"))
        ),
        rownames = FALSE
      )
    })
  }

  # -------------------Hypnogram statistics------------------------------------#

  output$table.hypno <- DT::renderDataTable({
    req(values$hasstaging, values$variable.staging)

    m <- as.data.frame(matrix(
      c(
        "TRT", "Total Recording Time, based on scored epochs (T0 – T6) (mins)",
        "TIB", "Time In Bed: Lights Off to Lights On (mins) (T1 – T5) (mins)",
        "SPT", "Sleep period time: Sleep Onset to Final Wake (T2 – T4) (mins)",
        "SPT_PER", "Persistent Sleep Period time: Persistent Sleep Onset to Final Wake (mins)",
        "LOT", "Lights On Time (mins)",
        "TWT", "Total Wake time during Lights Off = SLP_LAT + WASO + FWT (mins)",
        "TST", "Total Sleep Time (mins)",
        "TST_PER", "Total Persistent Sleep Time (mins)",
        " ", "",
        "SME", "Sleep maintenance efficiency, TST / SPT (denom. = T2 – T4)",
        "SE", "Sleep efficiency, TST / TIB (denom. = T1 – T5)",
        "WASO", "Wake time between sleep onset and final wake onset (T2 – T4)",
        "FWT", "Duration of wake from final wake onset to Lights on (T4 – T5)",
        "SOL", "Sleep latency (T1 – T2)",
        "SOL_PER", "Persistent sleep latency (T1 to onset of persistent sleep)",
        "REM_LAT", "REM latency (sleep onset T2 – first REM epoch)",
        "REM_LAT2", "REM latency excluding W, i.e. elapsed NR at REM onset",
        "NREMC", "Number of NREM cycles",
        "NREMC_MINS", "Mean NREM cycle duration (mins)",
        "  ", "",
        "CONF", "Number of epochs w/ conflicting stages (should be 0)",
        "SINS", "Recording starts in sleep (0=N, 1=Y)",
        "EINS", "Recording ends in sleep (0=N, 1=Y)",
        "OTHR", "Duration of non-sleep/non-wake annotations (unknown/movement)",
        "FIXED_WAKE", "Excessive leading/trailing Wake epochs set to L",
        "FIXED_LIGHTS", "Number of epochs set or L before/after lights out/on",
        "FIXED_SLEEP", "Sleep epochs set to L due to extreme WASO intervals",
        "LOT", "Lights On duration (mins)",
        "LOST", "Lights On Sleep duration (mins) : sleep set to L (should be 0)",
        "   ", "",
        "SFI", "Sleep Fragmentation Index: Sleep/W transition count / TST",
        "TI_S", "Stage Transition Index (excludes W): N1/N2/N3/R transition count / TST",
        "TI_S3", "3-class Stage Transition Index: NR/R/W transition count / SPT",
        "TI_RNR", "REM/NREM Transition Index, NR/R transition count / TST",
        "LZW", "LZW complexity index",
	"LZW3", "3-class LZW complexity index"
      ),
      ncol = 2, byrow = T
    ))

    # spacers
    values$opt[["hypno.stats"]][, " "] <- NA
    values$opt[["hypno.stats"]][, "  "] <- NA
    values$opt[["hypno.stats"]][, "   "] <- NA

    # add in values
    m$VALUE <- round(as.numeric(values$opt[["hypno.stats"]][, m[, 1]]), 3)

    DT::datatable(m,
      options = list(
        scrollY = "380px",
        dom = "tB", buttons = list(list(extend = "copy", text = "Copy")),
        ordering = F, paging = F,
        info = FALSE,
        searching = FALSE
      ),
      rownames = FALSE, colnames = rep("", 3)
    )
  })

  output$table.hypno.times <- DT::renderDataTable({
    req(values$hasstaging, values$variable.staging)

    m <- as.data.frame(matrix(
      c(
        "X0_START", "Study Start",
        "X1_LIGHTS_OFF", "Lights Off time (or start of recording)",
        "X2_SLEEP_ONSET", "Sleep Onset time",
        "X3_SLEEP_MIDPOINT", "Mid-point of T2 & T4",
        "X4_FINAL_WAKE", "Final Wake Onset time",
        "X5_LIGHTS_ON", "Lights On time (or end of recording)",
        "X6_STOP", "Study Stop"
      ),
      ncol = 2, byrow = T
    ))
    m$HMS <- as.character(values$opt[["hypno.stats"]][, gsub("X", "HMS", m[, 1])])
    m$E <- round(as.numeric(values$opt[["hypno.stats"]][, gsub("X", "E", m[, 1])]), 3)
    m$EPOCH <- floor(m$E * 2) + 1

    DT::datatable(m,
      options = list(
        ordering = F, pageLength = dim(m)[1],
        buttons = list(list(extend = "copy", text = "Copy")),
        lengthChange = FALSE,
        dom = "tB", info = FALSE,
        searching = FALSE,
        columnDefs = list(list(className = "dt-left", targets = "_all"))
      ),
      rownames = FALSE,
      colnames = c("Time-point", "Description", "Clock-time", "Elapsed (mins)", "Epoch")
    )
  })

  # ------------------------------------------------------------
  # Hypnogram stage statistics
  #

  output$table.hypno.stages <- DT::renderDataTable({
    req(values$hasstaging, values$variable.staging)

    dt <- values$opt[["hypno.stages"]]
    dt <- dt[dt$SS %in% c("?", "N1", "N2", "N3", "R", "S", "W", "WASO"), ]
    dt$EDUR <- as.integer(dt$MINS * 2)
    dt <- dt[, c("SS", "MINS", "EDUR", "PCT", "BOUT_MD", "BOUT_N")]
    dt$PCT <- round(100 * dt$PCT, 2)
    names(dt) <- c("Stage", "Duration (m)", "Duration (e)", "Duration (%)", "Median-bout(m)", "N-bouts")
    DT::datatable(dt,
      options = list(
        scrollY = "300px",
        dom = "tB", buttons = list(list(extend = "copy", text = "Copy")),
        info = FALSE,
        searching = FALSE,
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      ),
      rownames = FALSE
    )
  })

  # ------------------------------------------------------------
  # Hypnogram cycle statistics
  #

  output$table.hypno.cycles <- DT::renderDataTable({
    req(values$hasstaging, values$variable.staging)
    dt <- values$opt[["hypno.cycles"]]
    dt$NUM <- 1:(dim(dt)[1])
    dt <- dt[, c("NUM", "NREMC_START", "NREMC_N", "NREMC_MINS", "NREMC_NREM_MINS", "NREMC_REM_MINS")]
    names(dt) <- c("Cycle", "Start(E)", "Duration(E)", "Duration(m)", "NREM-duration(m)", "REM-duration(m)")
    DT::datatable(dt,
      options = list(
        scrollY = "300px",
        dom = "tB", buttons = list(list(extend = "copy", text = "Copy")),
        info = FALSE,
        searching = FALSE,
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      ),
      rownames = FALSE
    )
  })

  # ------------------------------------------------------------
  # Hypnogram epoch-statistics

  output$table.hypno.epochs <- DT::renderDataTable({
    req(values$hasstaging, values$variable.staging)

    dt <- values$opt[["hypno.epochs"]]
    dt <- dt[, c("E", "CLOCK_TIME", "MINS", "STAGE", "CYCLE", "PERSISTENT_SLEEP", "WASO", "E_N1", "E_N2", "E_N3", "E_REM", "E_SLEEP", "E_WASO")]
    names(dt) <- c("E", "Clock", "Mins", "Stage", "Cycle", "Per-Sleep", "WASO", "E(N1)", "E(N2)", "E(N3)", "E(REM)", "E(S)", "E(WASO)")
    #    dt$PCT <- round( dt$PCT , 2 )

    DT::datatable(dt,
      options = list(
        scrollY = "350px",
        dom = "tB", buttons = list(list(extend = "copy", text = "Copy")),
        info = FALSE, paging = F,
        searching = FALSE,
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      ),
      rownames = FALSE
    )
  })
  # -------------------- End of Hypnogram statistics ---------------------------#


  #--------- Initiate N=1 run from a N>1 panel selection ----------------------#

  observeEvent( input$runn1 , {
      req( input$procn1 )
      updateTabsetPanel(session, "tabset", selected = "N = 1" )
      reset(id = "upload1")

      # make a temp annot file (w/ times)
      txtPath <- tempfile(fileext = ".annot")

      # class(stage) . . start stop .
      stages <- values2$opt[["data"]][ values2$opt[["data"]]$ID %in% input$procn1 , c("SS","CLOCK_TIME") ] 
      stages <- data.frame( stages[,1] , INST = "." , CH = "." , START = stages[,2] , STOP = "+30" , META = "." ) 
      write.table( stages , file = txtPath , sep="\t" , row.names=F, col.names=F, quote=F ) 

      annot.names <- input$procn1
      annot.paths <- txtPath
      lset("epoch-check", "100000000")

      # Input file details
      values$file.details <- list(
        annot.names = annot.names,
        annot.paths = annot.paths
      )    
  })


  #-------------------------------- End of component1--------------------------#


  #-------------------------------- Component2---------------------------------#

  observeEvent(input$load.default2, {
    reset(id = "upload2")

    values2$file.details <-
      list(
        hypnos.names = "example.hypnos",
        hypnos.file = "data/example.hypnos.gz"
      )

   values2$covar.details <-
      list(
        covar.names = "example.tsv",
        covar.file = "data/example.tsv"
      )
    # Reset Shiny UI components
  })

  observeEvent(input$upload2, {
    values2$file.details <-
      list(
        hypnos.names = input$upload2$name,
        hypnos.file = input$upload2$datapath
      )
  })




  #---- Covariate data ------------------------------------------------------#

  observeEvent(input$upload3, {
   values2$covar.details <-
      list(
        covar.names = input$upload3$name,
        covar.file = input$upload3$datapath
      )
  })


  # Load covariates
  observeEvent(values2$covar.details, {
    values2$opt[["covarnames"]] <- values2$covar.details[["covar.names"]]

    # col 1 must be ID
    covars <- read.table( values2$covar.details[[ "covar.file" ]], header = T, stringsAsFactors = F , sep="\t" )
    if ( names( covars )[1] != "ID" ) stop( "need ID as first col" )

    covars <- covars[ covars$ID %in% values2$ids , ]
    covars <- covars[ order( covars$ID ) , ] 
    values2$covars <- covars
    
    init.covars()
  })


  init.covars <- function() {
    srts <- c( "Unsorted" , "Sleep Onset", "Start of recording" )
    updateSelectizeInput(session, inputId = "sort", choices = c( srts, names(values2$covars)[-1] ) , selected = NULL)
  }


  output$covar.table <- DT::renderDataTable({
    req(values2$covars )

    df <- values2$covars 
    df[is.na(df)] <- "."

    DT::datatable( df , escape = F, rownames = F ,
       options = list(
        scrollY = "480px",
        scrollX = "100%",
        buttons = list(list(extend = "copy", text = "Copy")),
        paging = F, ordering = F,
        info = T,
        searching = T) ) 
#        columnDefs = list(list(className = "dt-center", targets = 0:2) ) )
})



  #---- Process .hypnos files ------------------------------------------------#

  # Load data
  observeEvent(values2$file.details, {
    
    values2$opt[["hypnonames"]] <- values2$file.details[["hypnos.names"]]
    # assume 4 cols, ID E CLOCK_TIME STAGE;
    #  first col must be "ID" , otherwise names can change

    # optional: 5th column C == cycle
    #  use +10, +20 etc if epoch in 1st, 2nd, etc cycle
    #  first col must be "ID" , otherwise names can change

    d <- read.table(values2$file.details[["hypnos.file"]], header = T, stringsAsFactors = F)
    if ( names(d)[1] != "ID" ) stop( "first col must be ID" )
    if ( names(d)[2] != "E" ) stop( "second col must be E (epoch)" )

    # ensure original order is ID-sorted
    d <- d[ order(d$ID,d$E) , ] 

    if ( dim(d)[2] == 4 ) {
     names(d)[3] <- "CLOCK_TIME"
     names(d)[4] <- "SS"
     values2$cycles = 0
     d$CYCLE <- 0 
    } else {
     names(d)[3] <- "CLOCK_TIME"
     names(d)[4] <- "CYCLE"
     names(d)[5] <- "SS"
     values2$cycles = length( unique( d$CYCLE ) )
     d$CYCLE[ is.na( d$CYCLE ) ] <- 0 
    } 

    # save IDs (in order)
    values2$ids <- unique( d$ID ) 
    updateSelectizeInput(session, inputId = "procn1", choices = values2$ids , selected = NULL)

    ni <- length( values2$ids )
    nf <- length(d$ID[d$E == 1])
    if (ni != nf) {
      cat("not all indivs have first epoch (E==1)...\n")
      d <- d[d$ID %in% d$ID[d$E == 1], ]
    }

    # wipe any covariates
    values2$covars <- NULL
    init.covars()

    # Get first epoch for each individual
    d1 <- d[d$E == 1, ]

    # Convert HMS to seconds
    secs <- lubridate::period_to_seconds(lubridate::hms(d1$CLOCK_TIME))

    # for clocktime frames
    all.secs <- lubridate::period_to_seconds(lubridate::hms(d$CLOCK_TIME))
    sixes <- lubridate::period_to_seconds(lubridate::hms( c("18:00:00","06:00:00") ) )

    # Align to 30 sec epoch
    secs <- 30 * floor(secs / 30)
    all.secs <- 30 * floor( all.secs / 30)
    sixes <- 30 * floor( sixes / 30 )
    
    # Need clarification ( secs less than 12 hrs )
    secs[secs < 43200] <- secs[secs < 43200] + 86400
    all.secs[all.secs < 43200] <- all.secs[all.secs < 43200] + 86400
    sixes[ sixes < 43200 ] <- sixes[ sixes < 43200 ] + 86400
    
    # Get earliest and latest time points
    td <- seconds_to_period(min(secs))
    values$first.epoch <- sprintf("%02d:%02d:%02d", td@hour, minute(td), second(td))

    # Get epoch-wise starts for all people
    d1$E1 <- (secs - min(secs)) / 30
    values2$opt[["start_recording"]] <- d1$E1

    # Adjust epoch counts : EA = aligned epochs (by clock)
    d <- merge(d, d1[, c("ID", "E1")], by = "ID")
    d$EA <- d$E + d$E1

    # get last epoch (clock-time)
    first_clock <- min(secs) 
    last_clock <- max( all.secs )
#    cat( "fistlast clock " , first_clock , last_clock , "\n" )

    # scale 6pm/6am given secs
#    cat( "min( secs ) " , min(secs) , "\n" )
#    cat( "max( secs ) " , max(secs) , "\n" )
#    cat( "sx " , sixes , "\n" )
    values$sixpm <- ( sixes[1] - first_clock ) / ( last_clock - first_clock )
    values$sixam <- ( sixes[2] - first_clock ) / ( last_clock - first_clock ) 
#    cat( "sixes", values$sixpm , values$sixam , "\n" ) 

    # Get key anchors for each individual: sleep onset / offset, lights, sleep midpoint
    d$SLEEP <- as.integer(d$SS %in% c("N1", "N2", "N3", "R"))
    d1 <- as.data.frame(tapply(d$EA[d$SLEEP == 1], d$ID[d$SLEEP == 1], min))
    names(d1) <- "T2"
    d1$ID <- rownames(d1) # Timed epoch when subject starts to sleep

    # nb. do not assume all individuals will have sleep... thus all.x = T
    d <- merge(d, d1[, c("ID", "T2")], by = "ID", all.x = T)

    # Align to T2 == 0 (sleep onset) 
    d$E2 <- d$EA - d$T2

    # T2 alignment ticks - 
    values$onset_anchor <- ( 0 - min( d$E2 ) ) / ( max( d$E2 ) - min( d$E2 ) ) 
    values$onset_onehr  <-  120 / ( max( d$E2 ) - min( d$E2 ) ) 

    values2$opt[["tmp_data"]] <- d1
    values2$opt[["data"]] <- d
    values2$opt[["ni"]] <- ni
    values2$opt[["stgpal"]] <- c("#3D6D88", "#579DAF", "#0B2F38", "#CA7647", "#EBE3DD", lstgcols("?"))
  }) # End of observe Event (input$upload)

  # Create reactive expressions with reactive({ })
  # Reaction expression for Elapsed-time (ONSET) computation
  # Results are cached the first time the expression is called
  onset <- reactive({
    dmin <- tapply(values2$opt[["data"]]$E2, values2$opt[["data"]]$ID, min)
    dmax <- tapply(values2$opt[["data"]]$E2, values2$opt[["data"]]$ID, max)
    ids <- unique(values2$opt[["data"]]$ID)
    mindmin <- min(dmin)
    dmin <- dmin - mindmin + 1
    dmax <- dmax - mindmin + 1
    ne <- max(values2$opt[["data"]]$E2) - min(values2$opt[["data"]]$E2) + 1
    m <- matrix(NA, nrow = ne, ncol = values2$opt[["ni"]])
    for (i in 1:values2$opt[["ni"]])
     m[(dmin[i]):(dmax[i]), i] <- lstgn.hypnoscope( values2$opt[["data"]][ values2$opt[["data"]]$ID == ids[i] ,c( "SS", "CYCLE" ) ] )
    values$clocktime_mode = F
    m
  })

  # Reactive expression for Clock-time computation
  # Results are cached the first time the expression is called
  clockTime <- reactive({
    dmin <- tapply(values2$opt[["data"]]$EA, values2$opt[["data"]]$ID, min)
    dmax <- tapply(values2$opt[["data"]]$EA, values2$opt[["data"]]$ID, max)
    ids <- unique(values2$opt[["data"]]$ID)
    ne <- max(values2$opt[["data"]]$EA) - min(values2$opt[["data"]]$EA) + 1
    m <- matrix(NA, nrow = ne, ncol = values2$opt[["ni"]])
    for (i in 1:values2$opt[["ni"]])
     m[(dmin[i]):(dmax[i]), i] <- lstgn.hypnoscope(values2$opt[["data"]][values2$opt[["data"]]$ID == ids[i], c("SS","CYCLE")  ])
    values$clocktime_mode = T
    m
  })

  # note : cycle encoding : +10, +20 etc if in 1st, 2nd, etc cycle
  #  thus add dummy NAs to end... should never have codes 8/9/10 anyway
  stgpalW <- reactive({
    stcol <- c("#EBE3DD", "#EBE3DD", "#EBE3DD", "#EBE3DD", lstgcols("W"), "#EBE3DD", lstgcols("?"),NA,NA,NA)
    rep( stcol , values2$cycles + 1 ) 
  })

  stgpalN1 <- reactive({
    stcol <- c(lstgcols("N1"), "#EBE3DD", "#EBE3DD", "#EBE3DD", "#EBE3DD", "#EBE3DD", lstgcols("?"),NA,NA,NA)
    rep( stcol , values2$cycles + 1 )    
  })

  stgpalN2 <- reactive({
    stcol <- c("#EBE3DD", lstgcols("N2"), "#EBE3DD", "#EBE3DD", "#EBE3DD", "#EBE3DD", lstgcols("?"),NA,NA,NA)
    rep( stcol , values2$cycles + 1 )    
  })

  stgpalN3 <- reactive({
    stcol <- c("#EBE3DD", "#EBE3DD", lstgcols("N3"), "#EBE3DD", "#EBE3DD", "#EBE3DD", lstgcols("?"),NA,NA,NA)
    rep( stcol , values2$cycles + 1 )
   })

  stgpalNR <- reactive({
    stcol <- c( "black", "black", "black", "#EBE3DD", "#EBE3DD", "#EBE3DD", lstgcols("?"),NA,NA,NA)
    rep( stcol , values2$cycles + 1 )
   })

 stgpalR <- reactive({
    stcol <- c("#EBE3DD", "#EBE3DD", "#EBE3DD", lstgcols("R"), "#EBE3DD", "#EBE3DD", lstgcols("?"),NA,NA,NA)
    rep( stcol , values2$cycles + 1 )
  })

  stgpalL <- reactive({
    stcol <- c("#EBE3DD", "#EBE3DD", "#EBE3DD", "#EBE3DD", "#EBE3DD", lstgcols("L"), "#EBE3DD",NA,NA,NA)
    rep( stcol , values2$cycles + 1 )
  })

  stgpalU <- reactive({
    stcol <- c("#EBE3DD", "#EBE3DD", "#EBE3DD", "#EBE3DD", "#EBE3DD", "#EBE3DD", "#FF0000" ,NA,NA,NA)  # lstgcols("?") )
    rep( stcol , values2$cycles + 1 )
  })

  stgpalS <- reactive({
    stcol <- c("#000000", "#000000", "#000000", "#000000", "#EBE3DD", "#EBE3DD", "#FF0000" ,NA,NA,NA)  # lstgcols("?") )
    rep( stcol , values2$cycles + 1 )
  })

  stgpalC <- reactive({
    stcol <- rep( "#FFFFFF",10 ) 
    for (c in 1:values2$cycles ) stcol <- c( stcol , rep( c , 10 ) )
    stcol
  })

  default <- reactive({
    stcol <- c("#3D6D88", "#579DAF", "#0B2F38", "#CA7647", "#EBE3DD", lstgcols("L"), lstgcols("?"),NA,NA,NA)
    rep( stcol , values2$cycles + 1 )
  })

  output$text.header2a <- renderText({
    req(values2$file.details)
    values2$opt[["hypnonames"]]
  })


 #  
 # Row sort options
 #

  observeEvent(c(input$sort), {
    req(values2$file.details )    

    if ( input$sort == "Sleep Onset" )
     values2$sort_by <- order(values2$opt[["tmp_data"]]$T2, decreasing = T )
    else if ( input$sort == "Start of recording" )
     values2$sort_by <- order(values2$opt[["start_recording"]], decreasing = T )
    else if ( input$sort == "Unsorted" )
     values2$sort_by <- NULL
    else
     values2$sort_by <- order( values2$covars[ , input$sort ] , decreasing = T )
  })


  #
  # Observe Inputs and plot
  #

  observeEvent(c(values2$file.details, input$ultradian2, input$color, input$sort), {
    req(values2$file.details)

    values$clocktime_mode = input$ultradian2 == "Clock-time"

    data <- switch(input$ultradian2,
      "Clock-time" = clockTime(), # Call reactive expression
      "Elapsed-time" = onset()
    )

   # sort rows?
   if ( ! is.null( values2$sort_by ) )
    data <- data[, values2$sort_by ] 

   # palette
    stcol <- switch(input$color,
      W = stgpalW(),
      N1 = stgpalN1(),
      N2 = stgpalN2(),
      N3 = stgpalN3(),
      NR = stgpalNR(),
      R = stgpalR(),
      S = stgpalS(),
      L = stgpalL(),
      U = stgpalU(),
      Cycle = stgpalC(),
      All = default()
    )

    # Set image properties
    values2$opt["width"] <- 1150
    values2$opt["res"] <- 72

    if (ncol(data) > 500) {
      values2$opt["height"] <- ncol(data) * 1
    }
    if (between(ncol(data), 401, 500)) {
      values2$opt["height"] <- ncol(data) * 1.5
    }
    if (between(ncol(data), 301, 400)) {
      values2$opt["height"] <- ncol(data) * 1.5
    }
    if (between(ncol(data), 201, 300)) {
      values2$opt["height"] <- ncol(data) * 2
    }
    if (between(ncol(data), 151, 200)) {
      values2$opt["height"] <- ncol(data) * 3
    }
    if (between(ncol(data), 101, 150)) {
      values2$opt["height"] <- ncol(data) * 4
    }
    if (between(ncol(data), 51, 100)) {
      values2$opt["height"] <- ncol(data) * 5
    }
    if (between(ncol(data), 21, 50)) {
      values2$opt["height"] <- ncol(data) * 12
      values2$opt["width"] <- 1000
      values2$opt["res"] <- 60
    }
    if (between(ncol(data), 2, 20)) {
      values2$opt["height"] <- ncol(data) * 14
      values2$opt["width"] <- 800
      values2$opt["res"] <- 24
    }

    output$myImage <- renderImage(
      {

        # Generate the PNG
        outfile <- tempfile(fileext = ".jpeg")

        jpeg(outfile,
          width = as.numeric(values2$opt["width"]),
          height = as.numeric(values2$opt["height"]),
          res = as.numeric(values2$opt["res"]),
          quality = 100
        )

	# use full space
	par(mar=c(0,0,0,0))

        plot.new()

        # plot.window(xlim = c(-10, nrow(data)), ylim = c(-10, ncol(data)))
        # axis(side = 1, pos = 0, at = seq(from = 0, to = nrow(data), by = 1), col = "gray20",
        #     lwd.ticks = 0.25, cex.axis = 1, col.axis = "gray20", lwd = 1)
        # axis(side = 2, pos = 25, at = seq(from = 0, to = ncol(data), by = 1),
        #     col = "gray20", las = 2, lwd.ticks = 0.5, cex.axis = 1,
        #     col.axis = "gray20", lwd = 1.5)

	#
	# determine palette (based on number,sequence of observed stages
	#

	ustg <- as.integer(unique(names(table(data))))
	
	brks <- c( ustg[1] - 0.5 , ustg + 0.5 )

        image(data, xlim=c(-0.05,1) , ylim=c(-0.1,1) , 
          useRaster = T, col = stcol[ustg],
          xaxt = "n", yaxt = "n", axes = F, frame.plot=F ,
	  breaks = brks 
        )


	#
	# side plot
	#

	if ( any( names( values2$covars ) == input$sort ) )
	{
	 z <- values2$covars[ , input$sort ]
	 z <- z[ order(z, decreasing=T) ] 
	 z <- ( z - min(z,na.rm=T) ) / ( max(z,na.rm=T) - min(z,na.rm=T) )
         ni <- length(z)
	 xx <- -( z / 22 ) # for approx -0.05 .. 0.00 scaling 
	 yy <- seq(0,1,length.out = ni )
         polygon( c(0,0,-0.05,xx) , c(1,0,0,yy) ,col="purple" , border=NA )
        }

        #
	# x-time axis (clock-time or elapsed time)
	#
	
	if ( values$clocktime_mode )
	{
 	 lines( c(0,1) , c(-0.02 , -0.02) , col="black", lwd=1 ) 	
	 onehr <- ( values$sixam - values$sixpm ) / 12 
         hrs <- values$sixpm + (0:17) * onehr
	 tms <- c( paste( 6:11,"pm",sep="") , "M" , paste( 1:11,"am",sep="") )
	 tms <- tms[ hrs >=0 & hrs <= 1 ]
	 hrs <- hrs[ hrs >=0 & hrs <= 1 ]
	 for (t in 1:length(tms)) {
	    lines( c(hrs[t],hrs[t]),c(-0.015,-0.03) ) 
            text( hrs[t], -0.03 , tms[t],pos=1,adj=0.5)
            } 
	 lines( c(hrs[ which( tms == "8pm" ) ] , hrs[ which( tms == "8pm" ) ] ) , c(0,1), lwd=0.5, lty=2 )
	 lines( c(hrs[ which( tms == "M" ) ] , hrs[ which( tms == "M" ) ] ) , c(0,1), lwd=1, lty=1 )
	 lines( c(hrs[ which( tms == "6am" ) ] , hrs[ which( tms == "6am" ) ] ) , c(0,1), lwd=0.5, lty=2 )
	}

	#
	# elapsed time
	#

	if ( ! values$clocktime_mode ) 
        {
         lines( c(0,1) , c(-0.02 , -0.02) , col="black", lwd=1 )
         hrs <- values$onset_anchor + seq(-6,12) * values$onset_onehr
         tms <- paste( (-6):12,"hr",sep="") 
         tms <- tms[ hrs >=0 & hrs <= 1 ]
         hrs <- hrs[ hrs >=0 & hrs <= 1 ]
         for (t in 1:length(tms)) {
            lines( c(hrs[t],hrs[t]),c(-0.015,-0.03) )
            text( hrs[t], -0.03 , tms[t],pos=1,adj=0.5)
            }
#         lines( c(hrs[ which( tms == "6hr" ) ] , hrs[ which( tms == "0hr" ) ] ) , c(0,1), lwd=1, lty=1 )
        }


        # save
        dev.off()


        # Return a list containing information about the image
        list(
          src = outfile,
          contentType = "image/jpeg",
          width = as.numeric(values2$opt["width"]),
          height = as.numeric(values2$opt["height"]),
          alt = "This is alternate text"
        )
      },
      deleteFile = TRUE
    )
  })

  # server - Print number of observations plotted
  output$n <- renderUI({
    req(values2$file.details)

    data <- switch(input$ultradian2,
      "Clock-time" = clockTime(), # Call reactive expression
      "Elapsed-time" = onset()
    )
    HTML(paste0(
      "Number of observations: ", " <b>", ncol(data), "</b>"
    ))
  })



 observeEvent(input$hover1, {
  req(input$hover1)
#  cat( input$hover1$x , input$hover1$y , values2$opt[["height"]], "\n" )

  idx <- input$hover1$y / values2$opt[["height"]]

  # adjust height by 10% due to margin
  idx <- idx * 1.1
  if ( idx <= 1 )
   {
   ids <- values2$ids
   if ( ! is.null( values2$sort_by ) )
    ids <- ids[ values2$sort_by ] 
   # y-axis is upside in image()...
   ids <- rev( ids )  
   ni <- length( ids )
   idx <- min( max(1, as.integer(idx*ni) ) , ni )
#   cat( idx , ids[idx], "\n" )
   values2$h1.stgs <- values2$opt[["data"]]$SS[ values2$opt[["data"]]$ID == ids[idx] ]
   values2$h1.id <- ids[idx]
   } else values2$h1.stgs <- values2$h1.id <- NULL
})

 output$hyp1 <- renderPlot({
   req( values2$h1.stgs )
   req( length( values2$h1.stgs )> 2 ) 
   par(mar = c(0, 0, 0, 0))
   lhypno.mini( values2$h1.stgs , values2$h1.id  )
 })




}





