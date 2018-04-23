\name{qquniformfit}
\alias{qquniformfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ QQUniformFit }
\description{
  Quantile/Quantile-Plot compared to uniform / rectangular distribution with best-fit-line
}
\usage{
qquniformfit(x, FitSymbol, pstyle, xug, xog, Name, Plot, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Variable to be drawn }
  \item{FitSymbol}{ Optional: plot symbol for the interpolation line }
  \item{pstyle}{ Line type; '.' if not declared }
  \item{xug}{ Lower limit of the interval [xug, xog], in which a line will be interpolated; default=-3 }
  \item{xog}{ Upper limit of the interval [xug, xog], in which a line will be interpolated; default=3 }
  \item{Name}{ Description of the y-axis; 'data' if not declared }
  \item{Plot}{ Plot will be drawn only when this argument's value isn't NULL. Otherwise will be drawn in an .eps file. Plot=1 if not declared }
  \item{filename}{ Filename for the .eps file in which the plot will be saved }
}
\details{  }
\value{   }
\references{  }
\author{databionics UnivMarburg}
\note{   }
\seealso{  }
\examples{ }

