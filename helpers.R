lhypno2 <- function(hypno, cycles = NULL, times = seq(0, by = 30, length.out = length(ss)), start = 0, stop = max(times)) {
  ss <- hypno$STAGE
  ss[is.na(ss)] <- "?"
  e <- times / 3600
  sn <- lstgn(ss)
  
  plot(e, sn, type = "n", lwd = 2, col = "gray", axes = F, ylim = c(-3, 3.5), ylab = "", yaxt = "n", xaxs = "i", xlim = c(start, stop) / 3600, xlab = "Time (hrs)")
  # change points
  chgs <- which(ss[1:(length(ss) - 1)] != ss[2:length(ss)])
  for (chg in chgs) {
    # do not plot connector if change spans a gap; gap define assuming 30-second epochs
    if (!(times[chg + 1] - times[chg] > 40)) {
      lines(rep(((times[chg] + times[chg + 1]) / 2) / 3600, 2), c(sn[chg], sn[chg + 1]), lwd = 2, col = "gray")
    }
  }
  points(e, sn, col = lstgcols(ss), type = "p", cex = 1, pch = 20)

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
