slog <- function (x, m = 10){
  absX = abs(x)
  s = sign(x)
  if (m == 0) {
    return(s * log1p(absX))
  }
  if (m == 2) {
    return(s * log2(absX + 1))
  }
  if (m == 10) {
    return(s * log10(absX + 1))
  }
  else {
    factor = log(m)
    return(s * log1p(absX) * log(m))
  }
}


nanmin <- function (Data){
  if (length(dim(Data)) == 2) {
    SpaltenMinima <- apply(Data, 2, function(x) min(x, na.rm = TRUE))
    SpaltenInd <- NaN
  }
  else {
    SpaltenMinima <- min(Data, na.rm = TRUE)
    SpaltenInd <- which(Data == SpaltenMinima)
  }
  return(SpaltenMinima)
}

nanmax <- function (Data){
  if (length(dim(Data)) == 2) {
    SpaltenMinima <- apply(Data, 2, function(x) max(x, na.rm = TRUE))
    SpaltenInd <- NaN
  }
  else {
    SpaltenMinima <- max(Data, na.rm = TRUE)
    SpaltenInd <- which(Data == SpaltenMinima)
  }
  return(SpaltenMinima)
}

findAttrCol <- function (AttrName, Header, Data){
  i <- which(Header == AttrName)
  if (length(i) == 1) {
    AttrCol = Data[, i]
  }
  else {
    AttrCol = NaN
    warning("findAttrCol: No or multiple hits for \"", AttrName,
            "\" !")
  }
  return(AttrCol)
}

findAttrInd <- function (AttrName, VarNames){
  if (is.vector(VarNames)) {
    VarNames <- as.matrix(VarNames)
  }
  if (is.matrix(VarNames)) {
    x <- dim(VarNames)
    AnzVar <- x[1]
    VarNameLength <- x[2]
  }
  else {
    stop("findAttrInd(): VarNames must be a vector or a matrix!")
  }
  AttrNameLength <- length(AttrName)
  AttrInd <- 0
  if (AttrNameLength > VarNameLength) {
    return(AttrInd)
  }
  for (i in c(1:AnzVar)) {
    z <- grep(VarName[i, ], AttrName)
    if (z == 1) {
      AttrInd <- i
      return(AttrInd)
    }
  }
}


qquniform <- function (x, Name = "Data", pstyle = 20, filename = "") {
  x <- x[c(which(!is.nan(x)))]
  y <- qunif((1:length(x))/length(x), min = min(x), max = max(x))
  if (filename != "") 
    postscript(paste(filename, ".eps", sep = ""))
  qqplot(x, y, xlab = "uniform", ylab = Name, main = "Q/Q Plot compared to uniform distribution", 
         pch = pstyle, col = "blue")
  if (filename != "") 
    dev.off()
  return(y)
}

qquniformfit <- function (x, FitSymbol = "p", pstyle = 1, xug = 0, xog = 0, Name = "Data", 
                          Plot = 1, filename = "") {
  if (filename != "") 
    postscript(paste(filename, ".eps", sep = ""))
  x <- sort(na.last = T, x)
  y <- qquniform(x, Name, pstyle)
  z <- lm(y ~ x)
  abline(z, col = "red", lwd = 2)
  if (filename != "") 
    dev.off()
}


qqlognorm <- function (x, VarName = "Data", pstyle = 1, MU, SIGMA) {
  if (missing(MU)) {
    if (nanmedian(x) > 0) {
      MU = log(nanmedian(x))
    }
    else {
      stop("No MU given and median of input data is 0 or less, so no estimation possible")
    }
  }
  if (missing(SIGMA)) {
    if (mean(x, na.rm = T) > 0) {
      SIGMA = sqrt(2 * (log(mean(x, na.rm = T)) - MU))
    }
    else {
      stop("No SIGMA given and mean of input data is 0 or less, so no estimation possible")
    }
  }
  n = length(x)
  RandLogNorm = rlnorm(n, MU, SIGMA)
  qqplot(x = x, y = RandLogNorm, xlab = cat("LogNorm(", MU, 
                                            " , ", SIGMA, " )", sep = ""), ylab = VarName, main = "QQplot for LogNormal distribution", 
         pch = pstyle)
  return(invisible(list(qqx = x, qqy = RandLogNorm, MU = MU, 
                        SIGMA = SIGMA)))
}

