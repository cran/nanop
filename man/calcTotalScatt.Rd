\name{calcTotalScatt}
\Rdversion{1.1}
\alias{calcTotalScatt}
\alias{IqSAS}
\alias{IqSASP}

\title{Functions to calculate the total scattering structure function}
\description{
Functions to calculate the total scattering structure function given a matrix in which each row represents the position of an atom or to simulate it in the SAS region using particle lattice and size parameters. 
}
\usage{
calcTotalScatt(nanop, dQ=.01, minQ=0.771, maxQ=35, type="neutron", 
    scatterFactor=NA, scatterLength=NA, sigma=NA, n=0, 
    delta=0, kind="fastHist", dr = 0.001,  del = 0.01, eps=1e-3)               
IqSAS(Q, Rcore=NA, Rpart, latticep, latticepShell=NA, 
    scatterLength, N1, N2=NA, pDimer=0, sym, symShell=NA)	
IqSASP(Q, shell=NA, Rpart, latticep, latticepShell=NA, 
    scatterLength, N1, N2=NA, pDimer=0, sym, symShell=NA, rsigma)
}
\arguments{
  \item{nanop}{numeric matrix in which each row gives the coordinates of an atomic position in the nanoparticle. If nanop is not an object returned by \code{\link{simPart}} or \code{\link{displacePart}} then attributes \code{dimer}, \code{nAtomTypes}, \code{atomType}, \code{layer_start}, \code{layer_end}, \code{layerS_start}, \code{layerS_end}, \code{r} and \code{rcore} should be set manually; see \code{\link{simPart}}. 
}
  \item{dQ}{numeric indicating the desired step size in Q.}
  \item{minQ}{numeric indicating the minimum value of Q for which the function should be evaluated.}
  \item{maxQ}{numeric indicating the maximum value of Q for which the function should be evaluated.}
  \item{type}{character indicating type of scattering. Either "X-ray" or "neutron".}
  \item{scatterFactor}{list containing X-ray scattering factor parameters; see details. If \code{NA} the value is sought in \code{nanop} attributes.}
  \item{scatterLength}{in \code{calcTotalScatt} function call, numeric vector describing neutron scattering lengths for all atom types in the particles. If \code{NA} the value is sought in \code{nanop} attributes. In \code{IqSAS} and {IqSASP} calls, numeric vector describing average neutron scattering length for the particle core and shell; see examples.}
  \item{sigma}{numeric vector which, if not \code{NA}, determines the variances of the Gaussian displacements from the mean atomic positions throughout the nanoparticle. If \code{NA} the value is sought in \code{nanop} attributes. See \code{\link{simPart}} for details.}
  \item{n, delta}{numerics describing the correlation parameters \eqn{n} and \eqn{\delta} for thermal atomic displacements; see details.}
  \item{kind}{character, can be set to "exact", "fast", "fast_av" or "fastHist". If "exact" the Debye sum is calculated as is, that can be relatively time-consuming. If "fast", "fast_av" or "fastHist" then the pseudogrid approach (Cervellino et al., 2006) is used. If "fastHist" the histogram bin approach is used to calculate interatomic distances and their multiplicities. Option "fast_av" should be used if \code{nanop} attribute \code{"r"} length greater than \code{1}, i.e. for polydisperse particles . Nonzero \code{dr} value in that case switches computational scheme to histogram bin approach. In most case "fastHist" is recommended for better performance.}
  \item{dr}{numeric. If \code{kind=}"fastHist" or "fast_av" describes histogram bin width.}
  \item{del}{numeric. If \code{kind}="fast" or "fast_av" describes the pseudolattice parameter \eqn{\Delta}.}
  \item{eps}{numeric. If kind="fast" or "fast_av" describes the step size \eqn{\epsilon} for interatomic distances.}
  \item{Q}{numeric vector, contains grid points on which the total scattering structure function should be evaluated.}
  \item{N1, N2}{numerics indicating number of atoms within the unit cell in the particle core and shell.}
  \item{Rcore}{numeric which, if not \code{NA}, determines the radius of the core.}
  \item{Rpart}{numeric indicating the radius of the nanoparticle}
  \item{latticep, latticepShell}{numeric vectors indicating the lattice parameter(s) for the core (shell); see \code{\link{simPart}}.}
  \item{pDimer}{numeric indicating probability of atom to form a cluster with its neighbour.}
  \item{sym, symShell}{characters describing the structure to be used; see \code{\link{simPart}}.}
  \item{shell}{numeric indicating shell thickness for the core/shell particle. For polydisperse particles shell is supposed to be of constant thickness.}
  \item{rsigma}{numeric indicating standard deviation in the log normal distribution of the particle (core) radius.}
}
\details{The X-ray scattering factor is approximated by the function
\deqn{f(s)=a_1 * exp(-b_1 s) + a_2 * exp(-b_2 s) + a_3 * exp(-b_3 s) + a_4 * exp(-b_4 s) + c}
where \eqn{s = (\frac{Q}{4\pi})^2}. The constants in the function are possible to specify as arguments.  In that case \code{scatterFactor} parameter should contain vectors \code{a1}, \code{b1}, \code{a2}, \code{b2}, \code{a3}, \code{b3}, \code{a4}, \code{b4}, \code{a5}, \code{b5} and \code{c}. Their default values correspond to the values for Au and Pd atoms.  

The correlated atomic displacement parameter for the atoms \eqn{\mu} and \eqn{\nu} is calculated as 
\deqn{ \sigma^2_{\mu, \nu} = (\sigma^2_{\mu}  + \sigma^2_{\nu} ) [ 1 - \frac{\delta}{r^n}].}

}

\note{IqSASP function calculates total scattering function for the polydisperse uniform particles and polydisperse core/shell particles with constant shell thickness.
}

\references{
Cervellino A, Giannini C, Guagliardi A. (2006): On the efficient evaluation of Fourier patterns for nanoparticles and clusters. \emph{J. Comput. Chem.} \bold{27}, 995--1008.
}

\value{\code{calcTotalScatt}: list with elements
\item{Q}{numeric vector of values at which the function was evaluated,}
\item{gQ}{numeric vector of function values.}

\code{IqSAS} and \code{IqSASP}: numeric vector of function values.
}

\seealso{
\code{\link{simPart}}, \code{\link{displacePart}}
}
\examples{

## simulate particle
Au <- createAtom("Cu")
Au$name <- "Au"
Pd <- createAtom("Cu")
Pd$name <- "Pd"

part <- simPart(list(Au), atomsShell=list(Pd), r=10, rcore=8)
gQ <- calcTotalScatt(part, type="neutron", sigma=c(.08, .012), 
    kind="fast", del=5e-4)
plot(gQ$Q, gQ$gQ, type="l")

## "fast_av" option to calculate scattering function
## for the array of polydisperse particles:
Zn <- createAtom("Zn")
S <- createAtom("S")
part <- simPart(atoms=list(Zn, S), r=c(8, 10, 11.2, 13.4, 14), 
    sym="hcp", latticep=c(4.3, 7.02))
gQ <- calcTotalScatt(part, type="neutron", sigma=c(.08, .012), 
    kind="fast_av")
plot(gQ$Q, gQ$gQ, type="l")

## total scattering function in small-angle region using 
##  parametric model:
Q=seq(0.001, 0.771, 0.005)
gQSAS <- IqSAS(Q=Q, Rpart=26, latticep=c(3.21, 5.21), 
    scatterLength=7.1, N1=2, sym="hcp")
plot(Q, log(abs(gQSAS)), type="l")

## total scattering function for polydisperse sample with 
## lognormal distribution
gQSAS2 <- IqSASP(Q=Q, shell=2.8, Rpart=11.2, latticep=c(4.08), 
    latticepShell=3.89, N1=4, N2=4, scatterLength=c(7.1, 8.3),
    sym="fcc", symShell="fcc", rsigma=1.1)
plot(Q, log(abs(gQSAS2)), type="l")

}

\keyword{simulation}

