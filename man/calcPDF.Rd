\name{calcPDF}
\Rdversion{1.1}
\alias{calcPDF}
\alias{calcQDepPDF}
\title{Functions to calculate the PDF}
\description{Functions to calculate the pair distribution function (PDF)
  and Q-dependent PDF given a matrix of atomic positions.  The latter is
  currently relatively time-consuming. 
}
\usage{
calcPDF(nanop, calpha = 1, dr = 0.01, minR = 1, maxR = 20, p = 1,
foranalcs=FALSE)

calcQDepPDF(nanop, dr=.1, minR=1, maxR=20, minQ=1, maxQ=20,
           a1 = 16.8819, b1=.4611, a2=18.5913, b2=8.6216,
           a3=25.5582, b3=1.48260, a4=5.86, b4=36.3956, c=12.0658) 
}
\arguments{
  \item{nanop}{numeric matrix in which each row gives the coordinates of
    an atomic position in the nanoparticle
}
  \item{calpha}{numeric indicating the concentration of atoms of the
type under consideration

}
  \item{dr}{numeric indicating the step size in r desired 

}
  \item{minR}{numeric indicating the minimum value of r for which the
  function should be evaluated

}
  \item{maxR}{numeric indicating the maximum value of r for which the
  function should be evaluated
}
  \item{p}{numeric indicating the atomic number density or an
  approximation thereof
}
\item{minQ}{numeric indicating the minimum value of Q for which the
  function should be evaluated
}
\item{maxQ}{numeric indicating the maximum value of Q for which the
  function should be evaluated
}
\item{a1}{numeric X-ray scattering factor parameter; see details
}
  \item{b1}{numeric X-ray scattering factor parameter; see details

}
  \item{a2}{numeric X-ray scattering factor parameter; see details

}
  \item{b2}{numeric X-ray scattering factor parameter; see details

}
  \item{a3}{numeric X-ray scattering factor parameter; see details

}
  \item{b3}{numeric X-ray scattering factor parameter; see details

}
  \item{a4}{numeric X-ray scattering factor parameter; see details

}
  \item{b4}{numeric X-ray scattering factor parameter; see details

}
\item{c}{numeric X-ray scattering factor parameter; see details

}
\item{foranalcs}{logical indicating whether the result will be used in
  combination with calculating analytical broadening for a core-shell
  particle using the function \code{\link{broadPDF}} }

}
\details{
The X-ray scattering factor is approximated by the function
\eqn{f(s)=a_1 * exp(-b_1 s) + a_2 * exp(-b_2 s) + a_3 * exp(-b_3 s) +
  a_4 * exp(-b_4 s) + c}
with \eqn{s = \frac{Q}{4\pi}}.  The constants in the function are
possible to specify as arguments.  Their default values correspond to
the values for Au.  
}
\value{List with elements
\item{Q}{numeric vector of values at which the function was evaluated}
\item{gQ}{numeric vector of function values}

}
\seealso{
\code{\link{simPart}}, \code{\link{displacePart}}
}
\note{
  \code{calcQDepPDF} and \code{calcPDF} do not use the same
  normalization; \code{calcQDepPDF} is currently very time-consuming. 
  }
\examples{
## simulate a particle 
aa <- simPart()
## uniform displacement of positions
ab <- displacePart(aa, sigma=.005)
## different displacement in core than shell 
ab1 <-  displacePart(aa, rcore=8, sigmacore=.005, sigmashell=.02)

## calculate and plot PDF associated with both particles 
ac <- calcPDF(ab)
ac1 <- calcPDF(ab1)
plot(ac$r, ac$gr, type="l")
lines(ac1$r, ac1$gr, col=2)

}

\keyword{manip}

