\name{calcPDF}
\Rdversion{1.1}
\alias{calcPDF}
\alias{calcQDepPDF}
\title{Functions to calculate the PDF}
\description{Functions to calculate the pair distribution function (PDF) and Q-dependent PDF given a matrix of atomic positions.  The latter is relatively time-consuming. 
}
\usage{
calcPDF(nanop, dr=.01, minR=.01, maxR=20, 
        scatterLength=NA, scatterFactor=NA,  
        type="neutron", Qmin=1e-16)                        
calcQDepPDF(nanop=NA, dr=.1, minR=1, maxR=20, dQ=.01, minQ=1, maxQ=20,
            verbose=0, subdivisions = 100, order=1000,
            rel.tol=.Machine$double.eps^.7,
            addNoise=FALSE, noiseFun=NA, 
            totalScattParams=list(), preTotalScat=NA, ...)
}
\arguments{
  \item{nanop}{numeric matrix in which each row gives the coordinates of an atomic position in the nanoparticle. For \code{calcPDF} calculations if nanop is not an object returned by \code{\link{simPart}} or \code{\link{displacePart}} attributes \code{nAtomTypes}, \code{atomType}, \code{layer_start}, \code{layer_end}, \code{layerS_start}, \code{layerS_end}, \code{rowcore}, \code{scatterLength}, and \code{scatterFactor} should be set manually; see \code{\link{simPart}}. 
}
  \item{dr}{numeric indicating the desired step size in r.}
  \item{minR}{numeric indicating the minimum value of r for which function should be evaluated.}
  \item{maxR}{numeric indicating the maximum value of r for which function should be evaluated.}
  \item{scatterLength}{numeric vector describing neutron scattering lengths for all atom types in the particles. If \code{NA} the value is sought in \code{nanop} attributes.} 
  \item{scatterFactor}{list containing X-ray scattering factor parameters, see \code{\link{calcTotalScatt}}. If \code{NA} the value is sought in \code{nanop} attributes.}
  \item{type}{character; type of scattering. Either "X-ray" or "neutron".}
  \item{Qmin}{numeric used to approximate Q-dependent X-ray scattering factor via its value at point \code{Q=Qmin}.}
  \item{dQ}{numeric indicating the desired step size in Q in total scattering function calculations.}
  \item{minQ}{numeric indicating the minimum value of Q for which the total scattering function should be evaluated.}
  \item{maxQ}{numeric indicating the maximum value of Q for which the total scattering function should be evaluated.}
  \item{verbose}{numeric; if greater than zero a status report is given after computing the function every \code{verbose} steps in r.}
  \item{subdivisions}{numeric, the maximum number of subintervals for the Fourier integral calculation.}
  \item{order}{numeric, order of Gauss-Legendre quadrature.}
  \item{rel.tol}{numeric, relative accuracy requested for the Fourier integral calculation.}
  \item{addNoise}{logical indicating whether noise should be added to the total scattering structure function.  If \code{addNoise=TRUE} then function \code{noiseFun} is applied to the approximation of the total scattering structure function at points \code{Q} determined by the numerical integration routine.}
  \item{noiseFun}{function used to add noise to the total scattering structure function, applied only if \code{addNoise=TRUE}. The first argument should be a numeric vector representing the total scattering structure function. Additional arguments may be passed via \code{...}.}
  \item{totalScattParams}{list contacting objects \code{sigma}, \code{n}, \code{delta}, \code{kind}, \code{dr}, \code{del}, \code{eps}, \code{type}, \code{scatterLength}, and \code{scatterFactor} that are passed to \code{calcTotalScatt} function for total scattering function calculation. See \code{\link{calcTotalScatt}} for details. Should be specified if \code{preTotalScat=NA}.}  
  \item{preTotalScat}{list with elements \code{Q} (grid points) and \code{gQ} (values of the total scattering structure function at the corresponding gridpoints). If \code{NA} \code{calcTotalScatt} parameter should be specified.}  
  \item{...}{additional arguments for \code{noiseFun}.}
}
\details{
If \code{preTotalScat} is not \code{NA} \code{calcQDepPDF} function calculates Fourier transform of the vector \code{gQ$preTotalScat}.  Otherwise it uses parameters specified in \code{totalScattParams} to generate total scattering function first.  
}
\value{
list with elements
\item{r}{numeric vector of values at which the function was evaluated,}
\item{gr}{numeric vector of function values.}
}
\seealso{
\code{\link{simPart}}, \code{\link{displacePart}}
}
\note{
  To normalize \code{calcPDF} in the same way as \code{calcQDepPDF} it should be multiplied by \eqn{4\pi r}.

  \code{calcQDepPDF} is currently rather time-consuming. 
  }
\examples{
## simulate a particle 
Cu <- createAtom("Cu")
part <- simPart(atoms=list(Cu), atomsShell=list(Cu), r=8, 
    rcore=6, latticep=4.08, latticepShell=3.89)
## uniform displacement of positions
partx1 <- displacePart(part, sigma=c(.005, 0.005))
## different displacement in core than shell 
partx2 <-  displacePart(part, sigma=c(.005, .02))

## calculate and plot PDF associated with both particles 
gr1 <- calcPDF(partx1, maxR=24, scatterLength=c(4.87, 7.97), dr=.02)
gr2 <- calcPDF(partx2, maxR=24, scatterLength=c(4.87, 7.97), dr=.02)

plot(gr1$r, gr1$gr, type="l")
lines(gr2$r, gr2$gr, col=2)

## calculate scattering fucntion
gQ <- calcTotalScatt(part, type="neutron", dr=.02,
    scatterLength=c(4.87, 7.97), sigma=c(.005, .02))
t1 <- which(gQ$Q > 30)[1]
t2 <- which(gQ$Q > 34.9)[1]
cut <- gQ$Q[t1:t2][which(abs(gQ$gQ[t1:t2])
    ==min(abs(gQ$gQ[t1:t2]) ))[1]]
## calculate Q-dependent PDF
gr3 <- calcQDepPDF(part, minR=0, maxR=24, dr=0.02, minQ=.771, maxQ=cut,
    verbose=100, preTotalScat=list(Q=gQ$Q, gQ=gQ$gQ) )
## compare results with that obtained by calcPDF:
## ...normalization:
gr2 <- 4*pi*gr2$r*gr2$gr
## calculate and subtract gamma term:
gQSAS <- calcTotalScatt(part, type="neutron",  minQ=0.001, 
    maxQ=0.771, dQ=0.005,  dr=.02, 
    scatterLength=c(4.87, 7.97), sigma=c(.005, .02))
gammaR <- calcQDepPDF(part, minR=0.01, maxR=24,
    maxQ=.771, minQ=0.001, dr=.02,
    verbose=100, preTotalScat=list(Q=gQSAS$Q, gQ=gQSAS$gQ))
gr2 <- gr2 - gammaR$gr

## plot pair distribution functions associated with both technique:
plot(gr1$r, gr2, type="l")
lines(gr3$r, gr3$gr, col=2, lwd=2)

}

\keyword{simulation}
