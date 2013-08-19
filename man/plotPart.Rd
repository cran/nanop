\name{plotPart}
\Rdversion{1.1}
\alias{plotPart}
\title{Draws a three-dimensional scatterplot}
\description{Function to visualize a nanoparticle using \href{http://cran.r-project.org/web/packages/rgl/index.html}{rgl} package.}
\usage{
plotPart(nanop, radius=0.4, legend=TRUE, col=NA, box=FALSE, 
    play=FALSE, atoms=NA, miller=NA, lattice=c(4.08))
}
\arguments{
  \item{nanop}{numeric matrix in which each row gives the coordinates of an atomic position in the nanoparticle. If nanop is not an object returned by \code{\link{simPart}} or \code{\link{displacePart}} attributes \code{nAtomTypes}, \code{atomType}, \code{r}, \code{sym}, and \code{symShell} must be set manually; see \code{\link{simPart}}. 
  }
  \item{radius}{numeric vector or single value. Each atom on a plot is represented by a sphere. \code{radius} defines the sphere radius (radii).}
  \item{legend}{logical indicating whether to print plot legend.}
  \item{col}{numeric vector defining colours to be used for plotted items. If vector \code{col} length does not correspond to number of atom types within the particle then standard colouring scheme is used. 
  }
  \item{box}{logical indicating whether to draw box and axes.}
  \item{play}{logical. If \code{TRUE} animation with constantly rotating particle is played.}
  \item{atoms}{character. If not \code{NA} specifies atoms to be displayed, see details.}
  \item{miller}{numeric vector, specifies Miller indices. If not \code{NA} only the plane in a particle core described by given indices is displayed. Should be given in a form \code{c(h, k, l)} for the non-hexagonal symmetry and \code{c(h, k, i, l)} for the hexagonal symmetry. Should be specified together with \code{lattice} parameter.}
  \item{lattice}{numeric vector indicating particle core lattice parameters. Should be given in the same form as in \code{\link{simPart}}.}
}

\value{a vector of object IDs.}
\details{
If only core (shell) atoms of a specific type to be plotted \code{atoms} should be set up to \code{"core X"} or \code{"shell X"}, respectively.  Character describing atom type \code{"X"} can be taken from attributes(part)$atomsCore or attributes(part)$atomsShell. 
}
\examples{
## rgl library demands graphical terminal to be available
## uncoment all plotPart() calls for 3D visualization
## simulate particle
Au <- createAtom("Cu")
Au$name <- "Au"
Pd <- createAtom("Cu")
Pd$name <- "Pd"
part <- simPart(atoms=list(Au), atomsShell=list(Pd), rcore=8)
## 3d scatter plot
#plotPart(part, col=c(2,4))

## increase number of atom types within the particle:
Zn <- createAtom("Zn")
S <- createAtom("S")
part <- simPart(atoms=list(Zn ,S), atomsShell=list(Au), r=14, 
    rcore=12, sym="hcp", symShell="fcc", latticep=c(4.3, 7.04), 
    latticepShell=4.08)
## 3d scatter plot
#plotPart(part, col=c(2,4,3))
## play animation:
#plotPart(part, col=c(2,4,3), play=TRUE)

## plot only shell particles
#plotPart(part, col=c(2,4,3), atoms="shell Au", play=TRUE)

part <- simPart(atoms=list(Zn ,S),r=20, sym="hcp", 
    latticep=c(4.3, 7.04))
## display plane normal to z-axis:
#plotPart(part, miller=c(0, 0, 0 ,1), lattice=c(4.3, 7.04))
##S atoms:
#plotPart(part, miller=c(0, 0, 0 ,1), lattice=c(4.3, 7.04), 
#    atoms = "core S")
## save picture in a file using rgl function:
#rgl.snapshot( filename = "plane0001 S atoms.png")

Na <- createAtom("Na")
Cl <- createAtom("Cl")
part <- simPart(atoms=list(Na,Cl), sym="fcc")
#plotPart(part, miller=c(1,0,1), box=TRUE, lattice=c(4.08))
## plot only Na atoms:
#plotPart(part, miller=c(1,0,1), box=TRUE, lattice=c(4.08), 
#    atoms = "core Na")
}
\keyword{visualization}
