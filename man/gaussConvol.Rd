\name{gaussConvol}
\Rdversion{1.1}
\alias{gaussConvol}
\title{Gauss convolution}
\description{Function to calculate the convolution of a given function with the Gaussian function} 
\usage{		   
gaussConvol(SQ, Q, Qdamp=0.0457, err=1e-6)
}
\arguments{
  \item{SQ}{numeric vector containing function to be convoluted.}
  \item{Q}{numeric vector containing grid points corresponding to \code{SQ}.}
  \item{Qdamp}{numeric indicating standard deviation parameter for the Gaussian function.}
  \item{err}{numeric, relative accuracy requested.}
}

\value{
  numeric vector of function values.
}
\seealso{
  \code{\link{calcTotalScatt}}
}

\examples{
## simulate a particle 
Cu <- createAtom("Cu")
part <- simPart(atoms=list(Cu), atomsShell=list(Cu), r=20, 
    rcore=14, latticep=4.08, latticepShell=3.89)

## calculate total scattering function
gQ <- calcTotalScatt(part, type="neutron",  minQ=0.771, maxQ=18,
    dQ=0.01, scatterLength=c(4.87, 7.97), sigma=c(.01, .01))
plot(gQ$Q, gQ$gQ, type="l")
						  
## simulate instrumental resolution effect
Q <- gQ$Q
gQ <- gaussConvol(SQ=gQ$gQ, Q=gQ$Q, Qdamp=0.061, err=1e-5)
lines(Q, gQ, col=2)		   
}

\keyword{simulation}

