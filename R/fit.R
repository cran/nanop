##
## Objective function
## this function can be minimized to find the parameter values that
## minimize the RSS for a given model and given data
## In our experience this is best done with the package DEoptim. 
rss <- function(par, data, avRes, dataType="PDF",  simPP=NA, simPar=NA,
                PDF.fixed=NA, TotalScatt.fixed=NA,
                verbose,  parscale=NA, skel=NA, con=TRUE,
                punish=FALSE, avResrcenter=FALSE,
                analBroad=FALSE) {
  if(!is.na(parscale)) 
    par <- relist( par / parscale, skel)
  else 
    par <- relist(par, skel)

  if(con) {
    if(length(c(par$rcore,par$r)) == 2) ## test for non-param
      if(par$rcore>par$r)
        par$rcore <- par$r
    if(length(c(par$sigmacore,par$sigmashell)) == 2)
      if(par$sigmacore>par$sigmashell)
        par$sigmashell <- par$sigmacore
     
  }
  else if(punish) { ## only makes sense now for core shell
    if(length(c(par$rcore,par$r,par$sigmacore,par$sigmashell)) == 4)
      if(par$rcore>par$r || par$sigmacore>par$sigmashell)
        return(10e15)
  }
  if(any(unlist(par)<0))
    return(10e15)

  r <- if(is.null(par$r)) NA else par$r
  rcore <- if(is.null(par$rcore)) NA else par$rcore
  sigc <- if(is.null(par$sigmacore)) NA else par$sigmacore
  sig <- if(is.null(par$sigma)) NA else par$sigma
  sigs <- if(is.null(par$sigmashell)) NA else par$sigmashell
  delta <- if(is.null(par$delta)) 0 else par$delta
  
  if(length(par$latticep)>0)
    latticep <- par$latticep
  else latticep <- simPar$latticep
   
  if(length(par$latticepshell)>0)
    latticepshell <- par$latticepshell
  else latticepshell <- simPar$latticepshell

  if(is.na(simPP)) 
    part <- simPart(sym=simPar$sym, latticep=latticep, rcore=rcore,
                    latticepshell=latticepshell, rcenter=simPar$rcenter, 
                    r=r)
  if(dataType=="PDF") {
    if(!analBroad) {
      mod <- 0
      for(i in 1:avRes) {
        if(avResrcenter)
          part <- simPart(sym=simPar$sym, latticep=latticep, rcore=rcore,
                          latticepshell=latticepshell, rcenter=simPar$rcenter, 
                          r=r)
        
        dPart <- displacePart(part, sigma=sig, sigmacore=sigc,
                              sigmashell=sigs, rcore=rcore, rcenter=simPar$rcenter)
        mod <- mod + calcPDF(dPart,
                             calpha=PDF.fixed$calpha,
                             dr=PDF.fixed$dr,
                             minR=PDF.fixed$minR,
                             maxR=PDF.fixed$maxR,
                             p=PDF.fixed$p)$gr    
      }
      mod <- mod/avRes
    }
    else{
      mod <- calcPDF(part, calpha=PDF.fixed$calpha,
                       dr=PDF.fixed$dr,
                       minR=PDF.fixed$minR,
                       maxR=PDF.fixed$maxR,
                     p=PDF.fixed$p, foranalcs = !is.na(rcore))
      mod <- broadPDF(mod, sigma=sig, sigmacore=sigc, sigmashell=sigs,
                      rcore=rcore, delta=delta)$gr
    }
  }
  else if(dataType=="TotalScatt") {
    mod <- 0
    for(i in 1:avRes) {
       if(avResrcenter)
         part <- simPart(sym=simPar$sym, latticep=latticep, rcore=rcore,
                         latticepshell=latticepshell, rcenter=simPar$rcenter, 
                         r=r)
       
       dPart <- displacePart(part, sigma=sig, sigmacore=sigc,
                             sigmashell=sigs, rcore=rcore,
                             rcenter=simPar$rcenter)
       mod <- mod + calcTotalScatt(dPart, dQ=TotalScatt.fixed$dQ,
                                   minQ=TotalScatt.fixed$minQ,
                                   maxQ=TotalScatt.fixed$maxQ,
                                   a1=TotalScatt.fixed$a1,
                                   b1=TotalScatt.fixed$b1,
                                   a2=TotalScatt.fixed$a2,
                                   b2=TotalScatt.fixed$b2,
                                   a3=TotalScatt.fixed$a3,
                                   b3=TotalScatt.fixed$b3,
                                   a4=TotalScatt.fixed$a4,
                                   b4=TotalScatt.fixed$b4,
                                   c=TotalScatt.fixed$c)$gQ
       
     }
    mod <- mod/avRes
  }
  rss <- sum((data-mod)^2)
  if(dataType =="PDF")
    rss <- rss * 100000
  else rss <- rss / 10e7
  if(verbose) {
    pp <- signif(c(unlist(par), rss),4)
    names(pp) <- c(names(par), "**RSS**")
    print.table(pp)
  }
  rss
}
