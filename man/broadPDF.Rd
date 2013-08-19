\name{broadPDF}
\Rdversion{1.1}
\title{Analytically broaden the PDF}
\alias{broadPDF}
\description{Analytically broaden the PDF using Gaussians} 
\usage{
broadPDF(pdfob, sigma=0,  delta=NA, n=NA, nAtomTypes=1)
}
\arguments{
  \item{pdfob}{A list with elements \code{r} and \code{gr} that represent distances and values of the PDF, respectively. For core/shell particles elements \code{gr_CCSS} and \code{gr_CS} should also be specified. A list of this form is returned by function \code{\link{calcPDF}}.
  }
  \item{sigma}{numeric vector which, if not \code{NA}, determines the variance of the Gaussian displacements from the mean atomic position throughout the nanoparticle. If the particle core and shell contain \code{Nc} and \code{Ns} different atom types, respectively, then the first \code{Nc} elements in vector \code{sigma} correspond to atoms within the core and the next \code{Ns} elements describe Gaussian displacements for the shell atoms. See examples section in \code{\link{simPart}} for details.
  }
  \item{delta, n}{numerics describing the correlation parameters \eqn{n} and \eqn{\delta} for thermal atomic displacements; see details.}
  \item{nAtomTypes}{number of different types of atoms in the particle.}
}
\value{A list with elements \code{r} and \code{gr} that
  represent distances and values of the PDF, respectively.}
\details{
The correlated atomic displacement parameter for the atoms \eqn{\mu} and \eqn{\nu} is calculated as 
\deqn{ \sigma^2_{\mu, \nu} = (\sigma^2_{\mu}  + \sigma^2_{\nu} ) [ 1 - \frac{\delta}{r^n}].  }
}
\note{
  This routine can be a faster way to account for thermal displacements than \code{\link{displacePart}}.
}
\examples{
## simulate particle
Cu1 <- createAtom("Cu", sigma=0.012) 
Cu2 <- createAtom("Cu", sigma=0.008) 
part <- simPart(atoms=list(Cu1), atomsShell=list(Cu2), r=20, 
    rcore=16, latticep=4.08, latticepShell=3.89)

## use a stochastic model for displacements 
partx <- displacePart(part, sigma=attributes(part)$sigma) 
gr1 <- calcPDF(partx, maxR=40) 

## use analytical broadening
gr2 <- calcPDF(part, maxR=40)
gr2 <- broadPDF(gr2, sigma=attributes(part)$sigma, nAtomTypes=2)

# plot PDFs calculated using both methods
matplot(gr1$r, cbind(gr1$gr, gr2$gr), type="l", lty=1, lwd=1:2)
}

\keyword{simulation}

