\name{broadPDF}
\Rdversion{1.1}
\title{Analytically broaden a PDF}
\alias{broadPDF}
\description{Analytically broaden a PDF using Gaussians} 
\usage{
broadPDF(pdfob, sigma=NA, sigmacore=NA, sigmashell=NA,
                     rcore=NA, delta=0)     
}
\arguments{
  \item{pdfob}{A list with elements \code{r} and \code{gr} that
    represent distances and values of the PDF, respectively.  A list of
    this form is returned from the function \code{\link{calcPDF}}.
  }
  \item{sigma}{numeric which, if not \code{NA}, determines the variance
    of the Gaussian displacements from the mean atomic position
    throughout the nanoparticle.  If this argument is not \code{NA} then
    the following three arguments should be \code{NA}.
  }
  \item{sigmacore}{numeric which, if not \code{NA}, determines the variance
    of the Gaussian displacements from the mean atomic position
    in the core of the nanoparticle; if this argument is not
    \code{NA} then both
    \code{rcore} and \code{sigmashell} must also be not \code{NA}.
  }
  \item{sigmashell}{numeric which, if not \code{NA}, determines the variance
    of the Gaussian displacements from the mean atomic position
    in the shell of the nanoparticle; if this argument is not
    \code{NA} then both
    \code{rcore} and \code{sigmacore} must also be not \code{NA}.
  }
  \item{rcore}{numeric which, if not \code{NA}, determines the radius of
    the core; if this argument is not \code{NA} then both
    \code{sigmacore} and \code{sigmashell} must also be not \code{NA}.
  }
  \item{delta}{Correlation parameter}
}
\value{A list with elements \code{r} and \code{gr} that
  represent distances and values of the PDF, respectively.}
\note{If you wish to use this function to broaden a core-shell particle,
  the particle should be simulated using \code{\link{simPart}} and
  \code{\link{calcPDF}} with \code{foranalcs=TRUE}, as shown
by example below. } 

\examples{


## Particle 
aa <- simPart()

## Using a stochastic model for displacements 
ax <- displacePart(aa, sigma=.02) 
axx <- calcPDF(ax) 

## Using analytical broadening
axx1 <- calcPDF(aa)
aB <- broadPDF(axx1, sigma=.02) 

# plot PDFs calculated using both methods
matplot(aB$r, cbind(axx$gr, aB$gr), type="l", lty=1, lwd=1:2)

#########
#  this example shows analytical broadening of the PDF for a core-shell
## particle; note the use of the foranalcs=TRUE argument 

## Particle 
aa <- simPart(r=9, rcore=5, foranalcs=TRUE)

## broaden the PDF
ax <- calcPDF(aa, foranalcs=TRUE) 
aB <- broadPDF(axx1, sigmashell=.02, sigmacore=.003) 

}

\keyword{manip}

