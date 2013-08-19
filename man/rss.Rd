\name{rss}
\Rdversion{1.1}
\title{Function to calculate the RSS between the simulated PDF or the total scattering structure function and target data}
\alias{rss}
\alias{simPartPar}
\alias{PDFPar}
\alias{TotalScattPar}
\description{Given a nanoparticle model and PDF or total scattering structure function data, this function calculates the PDF or total scattering structure function associated with the model. It then determines the residual sum of squares (RSS) between the simulated signals and the data. } 
\usage{
rss(par, dataG=NA, dataS=NA, dataSAS=NA, part=NA, type="neutron",
    simPar=NA, PDF.fixed=list(), TotalScatt.fixed=list(), 
    parscale=NA, skel=NA, con=TRUE, oneDW=FALSE, 
    punish=FALSE, gammaR=NA, rvector=NA, fixed=NA,
    wG=0.03, wSAS=0.03, avRes=1, pareto=FALSE)					 
}
\arguments{
  \item{par}{numeric vector of parameter values. These should be given along with the appropriate argument \code{skel}. See details and examples.}
  \item{dataG}{numeric vector, which, if not \code{NA}, determines the PDF data.}
  \item{dataS}{numeric vector, which, if not \code{NA}, determines the total scattering structure function data.}
  \item{dataSAS}{numeric vector, which, if not \code{NA}, determines the total scattering structure function data in SAS region.}
  \item{part}{either \code{NA} or the atomic positions in a nanoparticle given as a numeric matrix in which each row gives the (3D) position coordinates.  A matrix of this form is the return value of the function \code{\link{simPart}}.  Default value of \code{NA} means that the particle is to be re-simulated in each call to \code{rss} (necessary for stochastic model functions).}  
  \item{type}{character indicating type of scattering. Either "X-ray" or "neutron".}
  \item{simPar}{arguments are passed to the function \code{\link{simPart}} to simulate the nanoparticles if these values are not to be optimized; any of the following arguments may be specified: \code{atoms}, \code{atomsShell}, \code{sym},\code{symShell}, \code{latticep}, \code{latticepShell}, \code{r}, \code{rcore}, \code{box},  \code{ellipse}, \code{shell}, \code{rcenter}, \code{move}, \code{rotShell}, \code{rcenterShell} and \code{distr}. The latter should be set up to "normal" or "lognormal" to simulate polydisperse particles with normal or log-normal distribution, respectively.  See \code{\link{simPart}} for the other parameters notation details.}		   
  \item{PDF.fixed}{arguments are passed to the functions \code{\link{calcPDF}} and \code{\link{broadPDF}}. Any of the following arguments may be specified: 
    \code{dr}, \code{minR}, \code{maxR}, \code{scatterLength}, \code{scatterFactor}, \code{delta}, \code{n}, \code{Qmax}, \code{termRip}, \code{maxRTermRip}, \code{N1}, \code{N2}, and \code{Qdamp}; see \code{\link{calcPDF}} and \code{\link{GrSAS}}. If \code{termRip} is \code{TRUE} termination ripples are simulated with parameters \code{dr}, \code{Qmax} and \code{maxRTermRip} (see \code{\link{termRip}}). If \code{Qdamp} is not \code{NA} the PDF is multiplied by the function \eqn{exp(-Q_{damp}^2 r^2/2)}.} 		   
  \item{TotalScatt.fixed}{arguments are passed to the functions \code{\link{calcTotalScatt}}; any of the following arguments may be specified: 
    \code{dQ}, \code{minQ}, \code{maxQ}, \code{minQ_SAS}, \code{maxQ_SAS}, \code{dQ_SAS}, \code{scatterLength}, \code{scatterFactor}, \code{delta}, \code{n}, \code{dr}, \code{del}, \code{eps}, \code{kind}, \code{N1}, \code{N2}, \code{f2}, \code{f2}, \code{convolution}, \code{Qdamp}, \code{SASscale}, and \code{paramSASQ}; see \code{\link{calcTotalScatt}}. If \code{convolution} is \code{TRUE} function \code{\link{gaussConvol}} is called with \code{Qdamp} parameter. If \code{paramSASQ} is \code{TRUE} the total scattering function in SAS region is calculated using functions \code{\link{IqSAS}} or \code{\link{IqSASP}} (with \code{f1} and \code{f2} parameters indicating average scattering lengths in the particle core and shell). Parameter \code{SASscale} sets weighting scheme for total scattering function in SAS region residual. Should be "normal", "log" (for logarithmic scale) or "weighted" (to normalize the residual at each grid point with the function value). }
  \item{parscale}{either \code{NA} or a numeric vector of the same length as \code{par} indicating values by which the values in \code{par} should be divided for the purpose of parameter scaling. }
  \item{skel}{an object of class \code{relistable}. First the parameters to be optimized should be written as a named list of form \code{parameters<-list(a=1,b=2,c=3)}. Then the argument \code{par} can be given as \code{par=unlist(parameters)} and \code{skel} can be given as \code{skel=as.relistable(parameters)}.}
  \item{con}{logical indicating whether to reset \code{r} to be equal to \code{rcore} if \code{rcore > r}.}
  \item{oneDW}{logical indicating whether to use similar \code{sigma} parameter for all atom types within the particle.}
  \item{punish}{if the inequalities described above for \code{con} are violated, return a large RSS value (\code{10e15}). }
  \item{gammaR}{either "param" or a numeric vector of gamma baseline term data. If "param" the baseline is simulated with functions \code{\link{GrSAS}} and \code{\link{GrSASCS}}.} 
  \item{rvector}{numeric vector of grid points at which the PDF was evaluated.}
  \item{fixed}{list; if not \code{NA} indicates parameter names.}
  \item{wG, wSAS}{weights for the PDF and total scattering function in SAS region residuals.}
  \item{avRes}{numeric value; if the model for the data is stochastic, he number of particles to simulate; the model PDF or total scattering structure function is calculated for each particle and then averaged.}
  \item{pareto}{logical; if \code{TRUE} the result value is a vector of residuals for the given data sets that can be used for the subsequent Pareto analysis.}
}
\details{
The following parameters may be fitted in \code{rss} (i.e., specified in \code{par} or \code{skel}): 
\code{r}, \code{rsigma} (sets standard deviation in the normal or log-normal distribution with mean value \code{r}), \code{box}, \code{ellipse},  \code{bkg} (background line that should be extracted from the SAS data), \code{shell}, \code{rcore}, \code{scaleSq} (scale factor for the total scattering function or the PDF), \code{scaleSASq} (scale factor for the total scattering function in SAS region or the gamma baseline term), \code{pStack}, \code{pDimer}, \code{delta}, \code{sigma}, \code{latticep} and \code{latticepShell}. For wurtzite structure to fit \code{u} parameter (second atom z-coordinate) third element in the vectors \code{latticep} and \code{latticepShell} may be specified.
}
\note{
Functions
					
\code{simPartPar(sym = "fcc", symShell=NA, latticep = 4.08, r=10,  }

\code{    atoms=list(), atomsShell=list(), distr = "lognormal",}

\code{    rcenter=FALSE, latticepShell=NA, rcore=NA, shell=NA, move=TRUE,}

\code{    box=NA, ellipse=NA, rotShell=FALSE, rcenterShell=FALSE),}
			   
\code{PDFPar(dr=.01, minR=1, maxR=20, termRip=FALSE, Qmax=30, maxRTermRip=20,}

\code{    scatterLength=NA, n=2, delta=NA, Qdamp=NA, scatterFactor=NA, N1=4, N2=4),}

and

\code{TotalScattPar(dQ=.01, minQ=0.771, maxQ=35, minQ_SAS=0.001, maxQ_SAS=0.771, }

\code{    dQ_SAS=0.005, scatterLength=NA, scatterFactor=NA, dr = 0.001,} 

\code{    del = 0.01, eps=1e-3,  kind="fastHist",  N1=4, N2=4, f1=5, f2=5,}

\code{    SASscale="normal", convolution = FALSE, Qdamp=0.0457, delta=NA, n=2, paramSASQ=FALSE) }

return values suitable for \code{simPar}, \code{PDF.fixed} and \code{TotalScatt.fixed} parameters, respectively.
}

\value{If \code{pareto=FALSE} the R-factor: the square root of the sum of squared distances between the target data and simulated signal divided by a the sum of squared target data values. If \code{pareto=TRUE} a vector with R-factors as components.}

\examples{

part <- res_gQ <- res_gQ_WAS <- res_gQ_SAS <- list() 

## prepare polydisperse sample consisting of 30 particles:
size <- sort(rlnorm(30, meanlog = log(10), sdlog = log(1.1)))
base_Cu <- getBase("Cu")
Au <- createAtom("Au", base=base_Cu)

for(i in 1:30) {
  cat("r = ", size[i], "\n")
  part[[i]] <- simPart(atoms=list(Au), r=size[i], latticep=4.08, 
      rcenter=TRUE) 
  res_gQ_WAS[[i]] <- calcTotalScatt(part[[i]],minQ=.771, maxQ=20, 
      dr=0.01, dQ=0.02, sigma=c(0.02))							
  res_gQ_SAS[[i]] <- calcTotalScatt(part[[i]],minQ=.001, 
      dr=0.01, maxQ=.771, dQ=0.01, sigma=c(0.02))
  cat(i,"\n")
}

## calculate average values:
gQ_av_WAS <- 0
for(i in 1:length(res_gQ_WAS)) {
  gQ_av_WAS <- res_gQ_WAS[[i]]$gQ + gQ_av_WAS
}
gQ_av_WAS <- gQ_av_WAS/length(res_gQ_WAS)

gQ_av_SAS <- 0
for(i in 1:length(res_gQ_SAS)) {
  gQ_av_SAS <- res_gQ_SAS[[i]]$gQ + gQ_av_SAS
}
gQ_av_SAS <- gQ_av_SAS/length(res_gQ_SAS)

## calculate PDF and gamma baseline term
resSAS  <- calcQDepPDF(minR=0, maxR=30, dr=0.01, minQ=.001,
    maxQ=.771, verbose=100, 
    preTotalScat=list(Q=res_gQ_SAS[[1]]$Q,gQ=gQ_av_SAS))
  
mr <- which(res_gQ_WAS[[1]]$Q > 17)[1]
mxr <- which(res_gQ_WAS[[1]]$Q > 19)[1]
cuto <- res_gQ_WAS[[1]]$Q[mr:mxr][which(abs(gQ_av_WAS[mr:mxr])
    ==min(abs(gQ_av_WAS[mr:mxr]) ))[1] ]
resWAS <- calcQDepPDF(minR=0, maxR=30, dr=0.01, minQ=.771, 
    maxQ=cuto, verbose=100, 
    preTotalScat=list(Q=res_gQ_WAS[[1]]$Q,gQ=gQ_av_WAS))

## set boundaries for fitting procedure
iter_0 <- as.relistable(list(latticep=0, r=0, sigma=c(0), rsigma=0))
boundsL <- c(latticep=3.5,  r=10.0, sigma=c(.01), rsigma=1.1)
boundsU <- c(latticep=4.5,  r=14.0,  sigma=c(.03), rsigma=1.3)

## in order to estimate the parameters that were used to
## simulate the particles, the DEoptim package may be
## used. Install it, remove the comment symbols '#' below,
## and use a call like:
#library(DEoptim)

#resDE <- DEoptim(rss, lower = boundsL, upper = boundsU, 
#    oneDW=FALSE, type="neutron",                   
#    control=DEoptim.control(CR=0.85, F=0.7, NP=40, itermax=30, 
#        parallelType=1, packages = list("PerformanceAnalytics"),
#        parVar=list("rss", "simPart", "Unique", "calcPDF", 
#        "IqSASP", "GrSAS", "GrSASCS", "broadPDF", "termRip",
#        "getSymEl", "IqSAS", "displacePart", "calcTotalScatt")),       
#    dataS=NA, dataSAS=gQ_av_SAS, dataG=resWAS$gr, 
#    simPar = simPartPar(atoms=list(Au), rcenter=TRUE, 
#        move=TRUE, rot=FALSE),
#    gammaR = resSAS$gr, rvector = resSAS$r, skel=iter_0, 
#    PDF.fixed=PDFPar(minR=0, maxR=30, dr=.01, 
#        scatterLength=Au$scatterLength),
#    TotalScatt.fixed=TotalScattPar(minQ=0.771, maxQ=20, 
#        dQ=0.02, dQ_SAS=0.01, minQ_SAS=.001, maxQ_SAS=.771, 
#        scatterLength=Au$scatterLength, kind="fastHist", 
#        SASscale="normal", convolution = FALSE),
#    wG=1.0, wSAS=0.05, avRes=10, pareto=FALSE) 
#
## now resDE$optim contains estimates for the lattice parameter,
## particle radius, displacement variance sigma, and standard 
## deviation rsigma.
## show results:
#resDE$optim$bestmem

## package mco may be used to construct a Pareto front of 
## solutions. Note that calculations may take considerable time!

#library(mco)
#resMCO <- nsga2(rss, 4, 2, 
#    dataS=NA, dataSAS=gQ_av_SAS, dataG=resWAS$gr,  oneDW=FALSE,
#    simPar = simPartPar(rcenter=TRUE, move=TRUE, rot=FALSE),
#    gammaR = resSAS$gr, rvector = resSAS$r, skel=iter_0,
#    PDF.fixed=PDFPar(minR=0, maxR=30, dr=.01, 
#        scatterLength=scatterLength),
#    TotalScatt.fixed=TotalScattPar(minQ=0.771, maxQ=20, dQ_SAS=0.01, 
#        dQ=.02, minQ_SAS=.001, maxQ_SAS=.771, scatterLength=scatterLength,
#        kind="fastHist", SASscale="normal", convolution = FALSE),
#    wG=1.0, wSAS=0.2, avRes=10, pareto=TRUE,							 
#    constraints=NULL, cdim=0, popsize=40, generations=c(40), 
#    cprob=0.85, lower.bounds=boundsL, upper.bounds=boundsU)

## show results
#plot(resMCO, xlab="SAS(Q) residual", ylab="G(r) residual")			
#paretoSet(resMCO)



}
\keyword{fit}

