
# include cycle infor
#  no cycle = 0
#  cycle 1  +10
#  cycle 2  +20
#  etc

lstgn.hypnoscope <- function(ssc) {
    ss <- ssc[,1]

    ss[ss == "N1" | ss == "NREM1"] <- 1 
    ss[ss == "N2" | ss == "NREM2"] <- 2
    ss[ss == "N3" | ss == "NREM3"] <- 3
    ss[ss == "R" | ss == "REM"] <- 4
    ss[ss == "W" | ss == "wake"] <- 5
    ss[ss == "WASO" ] <- 6
    ss[ss == "L" ] <- 7
    ss[ss == "?" | ss == "U" ] <- 8
    ss[is.na(ss)] <- 8
    ss <- as.integer(ss)
    # add in cycle info
    ss <- ss + ssc[,2] * 10

    ss
}

lstgn.waso <- function (ss) 
{
    ss[ss == "N1" | ss == "NREM1"] <- -1
    ss[ss == "N2" | ss == "NREM2"] <- -2
    ss[ss == "N3" | ss == "NREM3"] <- -3
    ss[ss == "R" | ss == "REM"] <- 0
    ss[ss == "W" | ss == "wake" | ss == "WASO" ] <- 1
    ss[ss == "?" | ss == "L"] <- 2
    ss[is.na(ss)] <- 2
    as.numeric(ss)
}


lstgcols.waso <- function (s)  {
    as.vector(sapply(s, function(x) {
        ifelse(x == "NREM1" | x == "N1", rgb(0, 190, 250, 255, 
            maxColorValue = 255), ifelse(x == "NREM2" | x == 
            "N2", rgb(0, 80, 200, 255, maxColorValue = 255), 
            ifelse(x == "NREM3" | x == "N3", rgb(0, 0, 80, 255, 
                maxColorValue = 255), ifelse(x == "NREM4" | x == 
                "N3", rgb(0, 0, 50, 255, maxColorValue = 255), 
                ifelse(x == "REM" | x == "R", rgb(250, 20, 50, 
                  255, maxColorValue = 255), ifelse(x == "L", 
                  rgb(246, 243, 42, 255, maxColorValue = 255), 
                  ifelse(x == "wake" | x == "W" | x == "WASO" , rgb(49, 173, 
                    82, 255, maxColorValue = 255), rgb(100, 100, 
                    100, 100, maxColorValue = 255))))))))
    }))
}


lhypno.mini <- function( ss , ids )
{
  ss[is.na(ss)] <- "?"
  uids <- unique(ids)  
  mx <- tapply( ss , ids , length )
  plot( c(0,max(mx)+100) , c(0,length(mx)+1) , type="n" , ylab = "", yaxt = "n" , axes = F ) 
  for (i in 1:length(uids) ) {
    ee <- ss[ ids == uids[i] ]
    points( 100 + 1:length(ee) , rep( i , length(ee) ) , col = lstgcols.waso(ee) , pch = "|", cex=1 )
    text( 1 , i , paste( uids[i] , " (" , round(length(ee)/120,1)," hrs)", sep=""), pos=4,adj=0)
  }
}

lhypno2 <- function(hypno, cycles = NULL, times = seq(0, by = 30, length.out = length(ss)), start = 0, stop = max(times)) {
  ss <- hypno$STAGE
  ss[is.na(ss)] <- "?"
  e <- times / 3600
  sn <- lstgn.waso(ss)
  
  plot(e, sn, type = "n", lwd = 2, col = "gray", axes = F, ylim = c(-3, 3.5), ylab = "", yaxt = "n", xaxs = "i", xlim = c(start, stop) / 3600, xlab = "Time (hrs)")
  # change points
  chgs <- which(ss[1:(length(ss) - 1)] != ss[2:length(ss)])
  for (chg in chgs) {
    # do not plot connector if change spans a gap; gap define assuming 30-second epochs
    if (!(times[chg + 1] - times[chg] > 40)) {
      lines(rep(((times[chg] + times[chg + 1]) / 2) / 3600, 2), c(sn[chg], sn[chg + 1]), lwd = 2, col = "gray")
    }
  }
  points(e, sn, col = lstgcols.waso(ss), type = "p", cex = 1, pch = 20)

  start_spt <- head(which(hypno$SPT == 1), n = 1) * 30
  stop_spt <- tail(which(hypno$SPT == 1), n = 1) * 30

  axis(1, c(start, stop) / 3600, lab = c("   Recording Start", ""), col = "#A8A8A8")
  axis(1, c(start_spt, stop) / 3600, xlab = "Time (hrs)", lab = c("Sleep Onset", ""), line = 2.5, col = "#A8A8A8")


  axis(2, 2, "?", col.axis = "black", las = 2)
  axis(2, 1, "W", col.axis = lstgcols("W"), las = 2)
  axis(2, 0, "R", col.axis = lstgcols("R"), las = 2)
  axis(2, -1, "N1", col.axis = lstgcols("N1"), las = 2)
  axis(2, -2, "N2", col.axis = lstgcols("N2"), las = 2)
  axis(2, -3, "N3", col.axis = lstgcols("N3"), las = 2)
  if (!is.null(cycles)) {
    if (length(cycles) != length(ss)) stop("ss and cycles must be same length")
    cc <- unique(cycles)
    cc <- cc[!is.na(cc)]
    odd <- T
    for (i in cc) {
      xc <- range(e[cycles == i & !is.na(cycles)])
      if (odd) {
        rect(xc[1], 3, xc[2], 3.3, col = "orange")
      } else {
        rect(xc[1], 2.7, xc[2], 3, col = "purple")
      }
      odd <- !odd
    }
  }
}


extract_hms <- function(fileName) {
  sleep_stage <- c("W", "N1", "N2", "N3", "R", "L")
  n <- 3
  hms_regex <- "^([01]?[0-9]|2[0-3]):[0-5][0-9]:[0-5][0-9]$"
  pat <- paste0("^([^:]+(?::[^:]+){", n - 1, "}).*")

  conn <- file(fileName, open = "r")
  linn <- readLines(conn)
  for (i in 1:length(linn)) {
    line <- gsub("\t", " ", linn[i], fixed = TRUE)
    line2 <- strsplit(line, " ")[[1]]
    given_stage <- line2[1]
    isSleepStage <- given_stage %in% sleep_stage

    if (isSleepStage) {
      startTime <- line2[4]
      startTime <- sub(pat, "\\1", startTime)
      startTime <- gsub("\\..*", "", startTime)

      if (grepl(hms_regex, startTime)) {
        print(startTime)
        break
      }
    }
  }
  close(conn)
  return(startTime)
}
