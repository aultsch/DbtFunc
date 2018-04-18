slog <- function (x, m = 10)
{
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


nanmin <- function (Data)
{
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

nanmax <- function (Data)
{
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

findAttrCol <- function (AttrName, Header, Data)
{
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

findAttrInd <- function (AttrName, VarNames)
{
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
