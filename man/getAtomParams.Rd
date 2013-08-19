\name{getBase}
\Rdversion{1.1}
\alias{getBase}
\alias{getScatterLength}
\alias{getScatterFactor}

\title{Fractional coordinates and scattering parameters for a given chemical element}
\description{Get fractional coordinates, neutron scattering length and X-ray scattering factor.} 
\usage{		   
getBase(name)
getScatterLength(name)
getScatterFactor(name)
}
\arguments{
  \item{name}{name of chemical element; see details.}. 
}

\value{
  \code{getBase}: numeric matrix in which each row represents fractional coordinates of the corresponding atom within the primitive cell.
  
  \code{getScatterLength}: numeric describing neutron scattering length for the given element.
  
  \code{getScatterFactor}: list containing X-ray scattering factor parameters for the given element.
}
\seealso{
  \code{\link{createAtom}}
}
\details{
  \code{getScatterLength} contains data for neutron scattering lengths of the elements and their isotopes taken from \url{http://www.ncnr.nist.gov/resources/n-lengths/list.html}.   Parameter \code{name} should be given as in the corresponding table. 

  \code{getScatterFactor} contains data for X-ray scattering factors of the elements and their isotopes taken from Waasmaier et al., 1995.  Parameter \code{name} should be given as in table 1.  The atomic scattering factor is calculated using the method developed by  Waasmaier et al. that implies approximation by a function
    \deqn{f(s)=a_1 * exp(-b_1 s) + a_2 * exp(-b_2 s) + a_3 * exp(-b_3 s) + a_4 * exp(-b_4 s) + c}, with \eqn{s = \frac{Q}{4\pi}}.

  \code{getBase()} function contains information about fractional coordinates for certain simple structure.  Parameter \code{name} can be:
  
  "Cu" to specify fractional coordinates of atoms in monoatomic fcc lattice,

  "Fe" to specify fractional coordinates of atoms in monoatomic bcc lattice,

  "Na" to specify fractional coordinates of Na atoms in NaCl rock salt structure,

  "Cl" to specify fractional coordinates of Cl atoms in NaCl rock salt structure,

  "Ca" to specify fractional coordinates of Ca atoms in CaTiO3 ideal perovskite structure,

  "Ti" to specify fractional coordinates of Ti atoms in CaTiO3 ideal perovskite structure,

  "O3" to specify fractional coordinates of O atoms in CaTiO3 ideal perovskite structure,

  "Mg" to specify fractional coordinates of atoms in monoatomic hcp lattice,

  "Zn" to specify fractional coordinates of Zn atoms in ZnS wurtzite structure,

  "S" to specify fractional coordinates of S atoms in ZnS wurtzite structure.

}
\references{
Waasmaier D. and Kirfel A. (1995): New analytical scattering-factor functions for free atoms and ions. \emph{Acta Cryst.} A\bold{51}, 416--431.
}
\examples{
## get fractional coordinates for Zn
Cd_base <- getBase("Zn")  
## get scattering parameters for Cd
Cd_scL <- getScatterLength("Cd") 
Cd_scF <- getScatterFactor("Cd")

Se_base <- getBase("S")
Se_scL <- getScatterLength("Se")
Se_scF <- getScatterFactor("Se")
}
\keyword{simulation}

