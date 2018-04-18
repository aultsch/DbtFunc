\name{slog}
\alias{slog}

\title{ signed Logarithms}
\description{
 \code{slog} computes signed logarithms, by default natural logarithms,  common (i.e., base 10) logarithms, and binary (i.e., base 2) logarithms
}
\usage{
slog(x, m = 10)
}

\arguments{
  \item{x}{a numeric or complex vector}
  \item{m}{mantissa of log \cr if m=0: natural logarithms \cr if m=2: binary logarithms \cr if m=10: common logarithms \cr otherwise: other base}
}

\value{
  A vector of the same length as \code{x} containing the transformed values.
}
\details{
SignedLog = slog(x);  default log == log10

 SignedLog = slog(x,m);
 
 Logarithmus  mit Vorzeichen:
 
 SignedLog = sign(x) .* log(abs(x)+1) fuer LogBase =m, default m=10
 
 INPUT
 
 x           Data 
 
 OPTIONAL

 m           Basis des Log, wenn =0 oder nicht angegebenn => log10

 OUTPUT
 
SignedLog = sign(x) .* log(abs(x)+1) fuer LogBase =m, default m=10


}

\author{fm, MT}

