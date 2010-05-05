\name{simPart}
\Rdversion{1.1}
\alias{simPart}
\alias{displacePart}
\title{Functions to simulate the atomic positions in a nanoparticle}
\description{Functions to simulate the deterministic
  atomic positions in a nanoparticle and displace those positions
  stochastically to model thermal effects}
\usage{
simPart(sym="fcc", latticep=4.08, r=10, latticepshell=NA,
rcore=NA, rcenter=FALSE, center=c(0,0,0), foranalcs=FALSE)

displacePart(nanop, sigma=NA, sigmacore=NA, sigmashell=NA, rcore=NA,
rcenter=FALSE, center=c(0,0,0), latticep = 4.08)

}
\arguments{
  \item{sym}{character vector describing the structure to be used;
    currently only "fcc" for face-centered cubic is available.
  }
  \item{latticep}{numeric indicating the lattice parameter(s).}
  \item{r}{numeric indicating the radius of the nanoparticle. (Currently,
    only spherical nanoparticle models are implemented).
  }
  \item{nanop}{numeric matrix in which each row represents the position
    of an atom, e.g., as returned by \code{simPart}.
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
  \item{latticepshell}{numeric indicating the lattice parameter(s) associated
    with the shell.}
  \item{rcenter}{logical value indicating whether to choose the center
    of the particle at random
    within the unit cell}
  \item{center}{If \code{rcenter=FALSE}, the position at the particle center.}
  \item{foranalcs}{logical indicating whether the result will be used in
  combination with calculating analytical broadening for a core-shell
  particle using the function \code{\link{broadPDF}}. } 

}

\value{numeric matrix with three columns
  in which each row represents an atomic position. 

}
\examples{

## Deterministic particle 
aa <- simPart()

## if the package 'rgl' is installed, the lines calling
## 'spheres3d' can be uncommented; they will provide a 3D
## visualization of a particle.
# spheres3d(aa, radius=.4)

## Particle with uniform displacements 
aa <- simPart(r=11, latticep=5)
ax <- displacePart(aa, sigma=.02) 
# spheres3d(ax, radius=.4)

## Particle with displacements in the core different from displacements
## in the shell
aa <- simPart(r=10, latticep=5)
ax <- displacePart(aa, sigmacore=.001, sigmashell=.1, rcore=10) 
# spheres3d(ax, radius=.4)

## Particle with radius drawn from a log-Normal size-distribution
ax <- simPart(r=exp(rnorm(1, log(10), log(1.01))), rcenter=TRUE)
}

\keyword{manip}

