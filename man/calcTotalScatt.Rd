\name{calcTotalScatt}
\Rdversion{1.1}
\alias{calcTotalScatt}
\alias{calcRedTotalScatt}
\title{Functions to calculate the total scattering structure function
  and reduced total scattering structure function}
\description{
Functions to calculate the total scattering structure function 
and reduced total scattering structure function given a matrix in which
each row represents the position of an atom.
}
\usage{
calcTotalScatt(nanop, dQ=.01, minQ=1,maxQ=20, a1 = 16.8819, b1=.4611,
                              a2=18.5913, b2=8.6216, a3=25.5582,
                              b3=1.48260, a4=5.86, b4=36.3956,
                              c=12.0658)

calcRedTotalScatt(nanop, dQ=.01, minQ=1, maxQ=20, a1 = 16.8819,
           b1=.4611, a2=18.5913, b2=8.6216, a3=25.5582, b3=1.48260,
           a4=5.86, b4=36.3956, c=12.0658)
}

\arguments{
\item{nanop}{numeric matrix in which each row gives the coordinates of
    an atomic position in the nanoparticle
}
\item{dQ}{numeric indicating the step size in Q desired 
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
}
\details{The X-ray scattering factor is approximated by the function
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
\examples{
aa <- simPart()
ab <- displacePart(aa, rcore=8, sigmacore=.005, sigmashell=.01)
ac <- calcTotalScatt(ab)
plot(ac$Q, ac$gQ, type="l")

}

\keyword{manip}

