\name{createAtom}
\Rdversion{1.1}
\alias{createAtom}
\title{Describe fractional coordinates and scattering parameters of the basis atoms}
\description{Function to create object of type \code{atoms} to be passed to \code{\link{simPart}}. } 
\usage{		   
createAtom(name, base=NA, sigma=0.01, scatterLength=NA, 
    scatterFactor=NA)
}
\arguments{
  \item{name}{character describing the name of chemical element; see details.}
  \item{base}{numeric matrix in which each row represents fractional coordinates of the corresponding basis atom within the primitive cell; see examples. If \code{NA} \code{\link{getBase}} is called with \code{name} argument.}
  \item{sigma}{numeric vector which, if not \code{NA}, determines the variance of the Gaussian displacements from the mean atomic positions for given atom.}
  \item{scatterLength}{numeric describing neutron scattering lengths for given atom. If \code{NA} \code{\link{getScatterLength}} is called with \code{name} argument.}
  \item{scatterFactor}{list containig X-ray scattering factor parameters. If \code{NA} \code{\link{getScatterFactor}} is called with \code{name} argument.}
}

\value{
  list with elements \code{name}, \code{base}, \code{sigma}, \code{scatterLength}, \code{scatterFactor}. 
}
\seealso{
  \code{\link{getBase}}, \code{\link{getScatterLength}}, \code{\link{getScatterFactor}}
}
\details{
Parameter \code{name} can be set up to any character value; if other arguments are not specified it is used as an argument in functions \code{\link{getBase}}, \code{\link{getScatterLength}} and \code{\link{getScatterFactor}} calls.  Functions \code{getScatterLength} and \code{getScatterFactor} contatain information for most chemical elements and some of their isotopes and ions.  Fractional coordinates, however, does depend not on the type of chemical element but on the structure under considirations.  Therefore argument \code{name} in function \code{getBase} indicates crystal structure (see \code{\link{getBase}}).  For example passing character "Cl" as an argument results in fractional coordinates of Cl atoms in the rock salt crystal structure.  
}
\examples{
## create CdSe nanoparticle:
## CdSe particle has wurtzite structure
## Zn and S atoms could be used to describe
## basis atoms fractional coordinates
## (see getBase() )

## get fractional coordinates for Zn
Cd_base <- getBase("Zn")  
## get scattering parameters for Cd
Cd_scL <- getScatterLength("Cd") 
Cd_scF <- getScatterFactor("Cd")
Cd <- createAtom("Cd", base=Cd_base, scatterLength=Cd_scL, 
    scatterFactor=Cd_scF, sigma=0.01)

Se_base <- getBase("S")
Se_scL <- getScatterLength("Se")
Se_scF <- getScatterFactor("Se")

Se <- createAtom("Se", base=Se_base, scatterLength=Se_scL, 
    scatterFactor=Se_scF, sigma=0.008)

## atoms Cd and Se now can be used as arguments in simPart
part <- simPart(atoms=list(Cd, Se), latticep=c(4.3, 7.02), 
    sym="hcp")
## uncomment to plot particle:
#plotPart(part)

## create rock salt structure
Na <- createAtom("Na")
Cl <- createAtom("Cl")
## name, scatterLength and scatterFactor parameters
## in Na and Cl are now set up to those of Na and Cl.
## if we are interested only in particle structure
## or actually simulating rock salt there is no need
## to change them
part <- simPart(atoms=list(Na, Cl))
#plotPart(part)

## set up fractional coordinates that cannot be
## simulated by getBase()
## (for example we have ZnS wurtzite structure 
## with atom S z-coordinate 0.378
## different from ideal value of 0.375 (3/8) ):

S_base <- matrix(c(
      c(0.5,  sqrt(3)/6, 0.378),  #S
	  c(0.5,  -sqrt(3)/6,  (0.378-0.5))),
      ncol=3, byrow=TRUE)
## each row represents fractional atomic coordinates in the primitive cell
S2 <- createAtom("S", base=S_base)
## use standard parameters for Zn atoms
Zn <- createAtom("Zn")
part <- simPart(atoms=list(Zn, S2), latticep=c(4.3,7.02), sym="hcp")
#plotPart(part)
   
}

\keyword{simulation}

