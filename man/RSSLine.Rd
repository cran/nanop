\name{RSSLine}
\Rdversion{1.1}
\alias{RSSLine}
\title{Function to generate 1-D RSS Lines}
\description{Function to generate 1-D RSS lines that represent the
  residual sum of squares as a single parameter (\code{r}, \code{rcore}
or \code{latticep}) is varied.  
}
\usage{
RSSLine(ref,  rcore, sigc, sigs, r, latticep=4.08, 
                    calpha=1, dr=.01, minR=1, maxR=20, p = 1,
                    dQ=.01, minQ=1,maxQ=20, a1 = 16.8819,
                    b1=.4611, a2=18.5913, b2=8.6216,
                    a3=25.5582, b3=1.48260, a4=5.86, b4=36.3956,
                    c=12.0658, calcPDF=TRUE, calcScatt=TRUE, norm=TRUE,
                    avRes = 1, rcenter=TRUE) 

}
\arguments{
  Only one of \code{r}, \code{rcore}, \code{latticep} can have length
  greater than 1; the vector with length greater than 1 represents the values
  over which to evaluate the RSS, with the other parameters fixed.  The
  argument \code{ref} indicates the value of the varied parameter to use
  in simulation of the data.
   \item{ref}{numeric value indicating the value of the varied parameter
   to be used in simulating the data.}
  \item{rcore}{either a single numeric value or a vector of numeric values for the
    parameter rcore}
  \item{sigc}{a single numeric value for the
    parameter sigc}
  \item{sigs}{a single numeric value  for the
    parameter sigs}
  \item{r}{either a single numeric value or a vector of numeric values for the
    parameter r}
  \item{latticep}{either a single numeric value or a vector of numeric values for the
    parameter latticep}
  \item{dr}{numeric indicating the step size in r desired }
  \item{minR}{numeric indicating the minimum value of r for which the
  function should be evaluated

}
  \item{maxR}{numeric indicating the maximum value of r for which the
  function should be evaluated
}
  \item{p}{numeric indicating the atomic number density or an
  approximation thereof
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

} \item{avRes}{numeric (integer) giving the ensemble size used to
  calculate both the data and the model}
\item{rcenter}{logical indication
  whether center should be chosen at random within the
   unit cell; note that use of a random center results in more
   computational time, since particle positions are not recycled}
 \item{calpha}{numeric scaling factor for PDF determination}
  \item{calcPDF}{logical indicating whether the RSS line should be
    determined for a PDF model}
  \item{calcScatt}{logical indicating whether the RSS line should be
    determined for a total scattering structure function model}
  \item{norm}{logical indicating whether the results should be normalized}
}
\details{
The X-ray scattering factor is approximated by the function
\eqn{f(s)=a_1 * exp(-b_1 s) + a_2 * exp(-b_2 s) + a_3 * exp(-b_3 s) +
  a_4 * exp(-b_4 s) + c}
with \eqn{s = \frac{Q}{4\pi}}.  The constants in the function are
possible to specify as arguments.  Their default values correspond to
the values for Au.  
}
\value{List in which the first element are the values of the parameter
  that is varied, and additional elements 
\item{ansR}{numeric vector of values at which the RSS was evaluated
  for the PDF RSS function}
\item{ansQ}{numeric vector of values at which the RSS was evaluated
  for the total scattering structure function RSS function}

}
\seealso{
\code{\link{simPart}}, \code{\link{displacePart}}
}

\keyword{manip}

