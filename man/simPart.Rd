\name{simPart}
\Rdversion{1.1}
\alias{simPart}
\alias{displacePart}
\title{Functions to simulate the atomic positions in a nanoparticle}
\description{Functions to simulate the deterministic atomic positions in a nanoparticle and displace those positions stochastically to model thermal effects}
\usage{
simPart(atoms, sym = "fcc", latticep = 4.08, r=10,  
    atomsShell=NA, symShell = NA, latticepShell = NA, rcore = NA, 
    shell=NA, box=NA, ellipse=NA, pDimer=0, pStack=0, 
    rcenter=FALSE, center=c(0,0,0), 
    move=TRUE, rotShell=FALSE, rcenterShell=FALSE)
						
displacePart(nanop, sigma=NA, rcenter=FALSE, latticep=4.08)
}

\arguments{
  \item{atoms}{list with elements describing basis atoms in the primitive cell. Elements should be the same type as the return value of \code{\link{createAtom}}; see examples.}
  \item{sym}{character describing the structure to be used:
  
        "fcc" for face-centered cubic structure,
		
		"bcc" for body-centered cubic structure,
		
        "sc" for simple cubic structure,
		
        "hcp" for hexagonal close-packed structure.
  }
  \item{latticep}{numeric vector indicating the lattice parameter(s). Should be given in the form:
  
        \code{c(a)} for fcc, bcc and sc structure,
		
        \code{c(a, c)} for hcp structure.      
  }
  \item{r}{numeric vector indicating the radius (radii) of the nanoparticle. If length is greater than \code{1} particles for all given radii are simulated. See \code{Value} for details.}
  \item{atomsShell}{list with elements describing basis atoms in the primitive cell for the particle shell.}
  \item{symShell}{character describing the structure associated with the shell. See \code{sym} for details.}
  \item{latticepShell}{numeric vector indicating the shell lattice parameter(s). See \code{latticep} for details.}
  \item{rcore}{numeric vector which, if not \code{NA}, determines the radius (radii) of the core. Should be the same length as vector \code{r}.}
  \item{shell}{numeric which, if \code{rcore=NA}, determines the particle shell thickness.}
  \item{box}{numeric vector. If not \code{NA} particle is simulated in the form of parallelepiped with dimensions \code{box[1]}, \code{box[2]} and \code{box[3]}.}
  \item{ellipse}{numeric vector. If not \code{NA} particle is simulated in the form of ellipse; \code{ellipse[1]}, \code{ellipse[2]} and \code{ellipse[3]} give the three elliptic radii.}
  \item{pDimer}{numeric indicating probability of atom to form a cluster with its neighbour.}
  \item{pStack}{numeric indicating probability of stacking fault to appear in each atomic plane. Currently is not available for perovskite symmetry. Produces ABCACABC stacking for cubic lattice and ABABABCBC stacking for hexagonal lattice.}
  \item{rcenter}{logical value indicating whether to choose the center of the particle at random within the unit cell.}
  \item{center}{If \code{rcenter=FALSE}, the position at the particle center.}
  \item{move}{logical value indicating whether to shift atoms near the core-shell boundary to avoid unphysically small interatomic distances. Might not work correctly if \code{r} length is greater than \code{1}.}
  \item{rotShell}{logical value indicating whether to rotate the particle shell by a random angle with respect to the core.}
  \item{rcenterShell}{logical value indicating whether to choose the center of the particle shell at random. If \code{FALSE} then the shell center coincide with that of the core.}
  \item{sigma}{numeric vector which determines the variances of the Gaussian displacements from the mean atomic positions throughout the nanoparticle. If \code{NA} value is taken from \code{nanop} attributes. If the particle core and shell contain \code{Nc} and \code{Ns} different atom types, respectively, then the first \code{Nc} elements in vector \code{sigma} correspond to atoms within the core and the next \code{Ns} elements describe Gaussian displacements for the shell atoms. See \code{exapmles} for more details.
  }
  \item{nanop}{numeric matrix in which each row represents the position of an atom, e.g., as returned by \code{simPart}.}
}
\value{numeric matrix with three columns in which each row represents an atomic position and list of the following attributes that describe particle properties:
  \item{"nStacks"}{number of stacking faults simulated in the particle. }
  \item{"rowcore"}{number of atoms within the core. }
  \item{"rowshell"}{number of atoms within the shell. }
  \item{"center"}{position of the particle center. }
  \item{"nAtomTypes"}{number of different atom types in the particle, e.g. for a core/shell particle with "Cu"/"ZnS" symmetry \code{nAtomTypes=3}.  }
  \item{"atomType"}{vector of length \code{rowcore+rowshell} that flags atoms in the particle with their atom type numbers. Atoms in the core are numbered with notation \eqn{1,2, \ldots}, while atoms in the shell have numbers \eqn{-1,-2, \ldots}, e.g. for a core/shell particle with "Cu"/"ZnS" symmetry vector \code{atomType} consists of elements \code{{1, -1, -2}}   }
  \item{"scatterLength"}{numeric vector describing neutron scattering lengths for all atom types in the particle}. 
  \item{"scatterFactor"}{list describing X-ray scattering factor parameters. See \code{\link{calcTotalScatt}} for details.}
  \item{"sigma"}{ see \code{sigma} argument. }
  \item{"layer_end", "layer_start"}{numeric vectors that are used if vector \code{r} length is greater than \code{1}. In that case for optimization purposes only one particle with the biggest value of radius \code{r} is generated, together with the information on how to extract atomic positions for the subparticle with given radius. E.g. if \code{r=c(8, 11.1, 16)}, \code{nanop} is an object that contains resulting matrix, \code{"layer_start=c(1,11,72)"} and \code{"layer_end=c(1024,1065,2048)"} then \code{nanop[11:1065,]} gives atomic coordinates for the subparticle with radius \code{r=11.1}. } 
  \item{"layerS_end", "layerS_start"}{ numeric vectors. Indicates layers for the particle shell. } 
  \item{"dimer"}{logical indicating whether the particle is a cluster of two spherical particles.}
  \item{"r"}{If \code{shape="sphere"} particle radius (radii). If \code{shape="ellipse"} elliptic radii. If \code{shape="box"} parallelepiped half-dimensions. }
  \item{"rcore"}{If \code{shape="sphere"} core radius (radii). If \code{shape="ellipse"} core elliptic radii. If \code{shape="box"} core half-dimensions.}
  \item{"shape"}{particle shape. }
  \item{"sym"}{particle (core) symmetry. }
  \item{"symShell"}{particle shell symmetry. }
  \item{"atomsCore"}{names of the chemical elements in the particle core.}
  \item{"atomsShell"}{names of the chemical elements in the particle shell.}
  
}

\examples{
## Uncomment all plotPart() calls for 3D visualization

## create CdSe nanoparticle:
## CdSe particle has wurtzite structure
## Zn and S atoms could be used to create base matrix
## (see getBase() )
Cd_base <- getBase("Zn")  #get fractional coordinates for Zn
Cd_scL <- getScatterLength("Cd") #get scattering parameters for Cd
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
#plotPart(part)

## Deterministic particle 
## Particle with uniform displacements 
Cu <- createAtom("Cu")
part <- simPart(atoms=list(Cu), atomsShell=list(Cu), rcore=8, 
    latticep=5)
partx <- displacePart(part, sigma=.02) 
#plotPart(partx, radius=.4, box=TRUE)

## Particle with displacements in the core different 
## from displacements in the shell
## create rock salt structure
Zn <- createAtom("Zn")
S <- createAtom("S")
part <- simPart(atoms=list(Zn,S), atomsShell=list(Cd, Se), r=14, 
    rcore=10, sym="hcp", latticep=c(3.1, 4.1), symShell="hcp", 
    latticepShell=c(4.3, 7.02))
partx <- displacePart(part, sigma=c(.01, .005, .012, .012)) 
## first elements in sigma and correspond to Zn atoms,
## second - to S atoms, third - to Cd atoms, last - to Se atoms.
attributes(part)$atomType
## elements '1' indicate Zn atoms in the total matrix, '2' indicate  
## S atoms, '-1' indicate Cd atoms, and '-2' indicate Se atoms.
#plotPart(partx, radius=.4)

## Particles with radii drawn from a log-Normal size-distribution
## and constant thickness 0.8
r <- exp(rnorm(10, log(10), log(1.1)))
part <- simPart(r=r, shell=0.8, atoms=list(Cu), 
    atomsShell=list(Cu))
## particle attributes
attributes(part)
## Extract second subparticle
t1 <- attributes(part)$layer_start[2]
t2 <- attributes(part)$layer_end[2]
part2 <- part[t1:t2, ]
## In order to use part2 as an argument for calcTotalScatt() or  
## other functions certain attributes should be set up 
## (see calcTotalScatt).
## To avoid extraction of every subparticle use option 
##  kind="fast_av" in calcTotalScatt() function.

}
\keyword{simulation}
