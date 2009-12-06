\name{rss}
\Rdversion{1.1}
\title{Calculate the RSS associated with the fit PDF or total scattering
structure function data}
\alias{rss}
\description{Given a nanoparticle model and PDF or total scattering
  structure function data, this function calculates the PDF or total
  scattering structure function associated with the model.  It then
  determines the residual sum of squares (RSS) between the PDF or total
  scattering structure function and the data. } 
\usage{
rss(par, data, avRes, dataType="PDF",  simPP=NA, simPar=NA,
                PDF.fixed=NA, TotalScatt.fixed=NA,
                verbose,  parscale=NA, skel=NA, con=TRUE,
                punish=FALSE, avResrcenter=FALSE,
                analBroad=FALSE) 
}
\arguments{
  \item{par}{A numeric vector of parameter values. These should be given along
    with the appropriate arguement \code{skel}.  See the examples for
    details. 
  }
  \item{data}{A numeric vector of PDF or total scattering structure
    function data. 
  }
  \item{avRes}{A numeric value; if the model for the data is stochastic,
    this is the number of particles to simulate; the model PDF or total
    scattering structure function is calculated for each particle and
    then averaged.
  }
  \item{dataType}{Either \code{"PDF"} if \code{data} represent PDF data or
    \code{"TotalScatt"} if \code{data} represent total scattering
    structure function data. 
  }
  \item{simPP}{Either \code{NA} or
    the atomic positions in a nanoparticle given as a numeric
    matrix in which each row gives the (3D) position coordinates.  A
    matrix of this form is the return value of the function
    \code{\link{simPart}}.  Default value of \code{NA} means that the
    particle is to be re-simulated each call to \code{rss}
    (necessary for stochastic model functions). 
  }
  \item{simPar}{Arguments are passed to the function
    \code{\link{simPart}} to simulate the nanoparticles if these values
    are not to be optimized; any of the
    following arguments may be specified: \code{latticep}, 
    \code{latticepshell}, \code{sym}, \code{rcenter}.
  }
  \item{PDF.fixed}{Arguments are passed to the function
    \code{\link{calcPDF}}; any of the
    following arguments may be specified: \code{calpha}, 
    \code{dr}, \code{minR}, \code{maxR}, \code{p}.} 
  \item{TotalScatt.fixed}{Arguments are passed to the function
    \code{\link{calcTotalScatt}}; any of the
    following arguments may be specified: \code{dQ}, 
    \code{minQ}, \code{maxQ}, \code{a1}, \code{b1},
    \code{a2}, \code{b2},\code{a3}, \code{b3},
    \code{a4}, \code{b4}, \code{c}. }
  \item{verbose}{Logical indicating whether to print information.}
  \item{parscale}{Either \code{NA} or a numeric
    vector of the same length as \code{par} indicating values by which
    the values in \code{par} should be divided for the purpose of
    parameter scaling.
  }
  \item{skel}{an object of class \code{relistable}.  First the
    parameters to be optimized should be written as a named list of
    form \code{parameters <-list(a=1,b=2,c=3)}.
    Then the argument \code{par} can
    be given as \code{par=unlist(parameters)} and \code{skel} can be
    given as \code{skel=as.relistable(parameters)}.}
  \item{con}{Logical indicating whether to reset
    \code{sigmashell} to be equal to \code{sigmacore} if
    \code{sigmashell > sigmacore} and to reset   
    \code{r} to be equal to \code{rcore} if
    \code{rcore > r}.}
  \item{punish}{If the inequalities described above for \code{con} are
  violated, return a large RSS value (\code{10e15}). }
 
\item{avResrcenter}{Logical; if \code{TRUE} (which is possibly desirable
  only if \code{simPar$rcenter=TRUE}) 
  then a new particle with a random center is
  simulated each of \code{avRes} times.}
\item{analBroad}{If \code{TRUE} then broadening to account for thermal
motion is accomplished analytically; if \code{FALSE} then it is
accomplished using the function \code{\link{displacePart}}. This option
only applies if \code{dataType="PDF"}; if \code{dataType="TotalScatt"}
then \code{link{displacePart}} is always applied. }
}
\value{A numeric value indicating the RSS. }
\examples{
a <- simPart()
x <- calcPDF(a)
xx1 <- broadPDF(x, sigma=.0125)

## plot PDF data  
plot(xx1$r,xx1$gr, type="l")

## in order to estimate the parameters that were used to
## simulate the particles, the DEoptim package may be
## used.  Install it, remove the comment symbols '#' below,
## and use a call like: 

##library(DEoptim)

##ax <- as.relistable(list(latticep=0,r=10,sigma=0)) 

##res <- DEoptim(rss,
##               lower=c(latticep=3.5, r=5.5,sigma=.001),
##               upper=c(latticep=5.6, r=12,sigma=.1),
##               control=DEoptim.control(NP=50, VTR=.5), 
##               data=xx1$gr, dataType="PDF", punish=TRUE,
##               simPar = simPartPar(), 
##               PDF.fixed=PDFPar(), 
##               verbose=TRUE, skel=ax,
##               analBroad=TRUE)

## now res$optim contains estimates for the particle radius,
## displacement variance sigma, and lattice parameter 
}

\keyword{manip}

