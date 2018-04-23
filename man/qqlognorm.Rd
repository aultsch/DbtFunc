\name{qqlognorm}
\alias{qqlognorm}
\title{
QQ-Plot im Vergleich. zur LogNormalverteilung
}
\description{
Quantile/Quantile = QQ-Plot im Vergleich. zur LogNormalverteilung
}
\usage{
qqlognorm(x, VarName = "Data", pstyle = 1, MU, SIGMA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of data which may contain Na/NaN
}
  \item{VarName}{
Name of the data shown in the plot
}

  \item{pstyle}{
Integer defining the display style of the plot
}

  \item{MU}{
The µ of the log normal distribution
}

  \item{SIGMA}{
  the σ of the log normal distribution
}
}
\value{
  \item{qqx}{
    Vector of x coordinates of the qq plot
  }
  \item{qqy}{
    Vector of y coordinates of the qq plot  
  }
  \item{MU}{
    The µ of the log normal distribution
  }
  \item{SIGMA}{
    the σ of the log normal distribution
  }
}
\author{
Matlab: ALU 
R: FP
}