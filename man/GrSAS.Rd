\name{GrSAS}
\Rdversion{1.1}
\alias{GrSAS}
\alias{GrSASCS}
\title{Functions to calculate the gamma baseline term for PDF}
\description{Functions to calculate analytically the gamma baseline term given a particle lattice and size parameters.} 
\usage{		   
GrSAS(r, Rcore=NA, Rpart, latticep, latticepShell=NA, 
    N1, N2=NA, sym, symShell=NA)   
GrSASCS(r, Rcore=NA, Rpart, latticep, latticepShell=NA, 
    N1, N2=NA, sym, symShell=NA)
}
\arguments{
  \item{r}{numeric vector that contains grid points at which baseline term should be evaluated.}
  \item{Rcore}{numeric which, if not \code{NA}, determines the radius (radii) of the core.}
  \item{Rpart}{numeric indicating the radius (radii) of the particle.}
  \item{latticep, latticepShell}{numeric vectors indicating the lattice parameter(s) for the core( shell); see \code{\link{simPart}} for details.}
  \item{N1, N2}{numeric indicating number of atoms within the unit cell in the particle core (shell).}
  \item{sym, symShell}{characters describing the structure to be used in the particle core (shell) simulations; see \code{\link{simPart}}.}
}
\details{Function \code{GrSAS} can be used for both uniform and core/shell particles.  In the second case the uniform model is applied with scattering length density averaged through the nanoparticle.  Function \code{GrSASCS} calculates baseline term for core/shell particles using model described in Glatter, 1979.
}
\value{
  numeric vector of function values.
}
\seealso{
  \code{\link{calcPDF}}, \code{\link{calcQDepPDF}}
}
\references{
Glatter O. (1979): The interpretation of real-space information from small-angle scattering experiments. \emph{J. Appl. Cryst.}, \bold{12}, 166--175.
}
\examples{
## calculate baseline term for uniform particle			
gammaR <- GrSAS(r=seq(0.01, 30, 0.01), Rpart=15, 
    latticep=c(4.3, 7.02), sym="hcp", N1=4)
plot(seq(0.01, 30, 0.01), gammaR, type="l")

## compare with baseline computed as Fourier transform 
## of the total scattering function:
Zn <- createAtom("Zn")
S <- createAtom("S")
part <- simPart(atoms=list(Zn,S), r=15, latticep=c(4.3, 7.02), 
    sym="hcp")
gQSAS <- calcTotalScatt(part, type="neutron",  minQ=0.001, 
    maxQ=0.9, dQ=0.005)
gammaR2 <- calcQDepPDF(part, minR=0.01, maxR=30, dr=0.01, 
    maxQ=.85, minQ=0.001, verbose=20,
    preTotalScat=list(Q=gQSAS$Q, gQ=gQSAS$gQ))
lines(gammaR2$r, gammaR2$gr, col=2)
}

\keyword{simulation}

