broadPDF <- function(pdfob, sigma=NA, sigmacore=NA, sigmashell=NA,
                     rcore=NA, delta=0)          {
  dsig <- function(peakind, r, sigi, sigj, gr, delta) 
    dnorm(x=r,mean=r[peakind],sd=sqrt(sigi+sigj-(delta/r[peakind])*
                                sqrt(sigi)*sqrt(sigj)))*gr[peakind]
   
  r <- pdfob$r
  if(!is.na(sigma)) { # uniform 
    gr <- pdfob$gr 
    nz <- which(gr != 0)
   
    np <- rowSums(apply(as.matrix(nz), 1, dsig, r=r,
                        sigi=sigma, sigj=sigma, gr=gr, delta=delta)) 
    xx <- sum(gr) / sum(np)
    pdfob$gr <- np * xx 
    
  }
  else { # core-shell
    gr <- pdfob$gr
    
    grC <- pdfob$grC 
    nzC <- which(grC != 0)
   
    grS <- pdfob$grS 
    nzS <- which(grS != 0)
      
    grCS <- pdfob$grCS 
    nzCS <- which(grCS != 0)

    if(length(nzS)>0) 
      npC <- rowSums(apply(as.matrix(nzC), 1, dsig, r=r, sigi=sigmacore,
                           sigj = sigmacore, gr=grC, delta=delta))
    else
      npC <- rep(0, length(gr))
    scC <- sum(pdfob$grC) / sum(npC) 
    if(length(nzS)>0) 
       npS <- rowSums(apply(as.matrix(nzS), 1, dsig, r=r,
                            sigi=sigmacore, sigj=sigmashell, gr=grS,
                            delta=delta)) 
    else npS <- rep(0, length(gr))
    scS <- sum(pdfob$grS) / sum(npS) 
    if(length(nzCS)>0) 
      npCS <- rowSums(apply(as.matrix(nzCS), 1, dsig, r=r,
                            sigi=sigmacore, sigj=sigmashell,
                            gr=grCS, delta=delta))
    else npCS <- rep(0, length(gr))
   
    scCS <- sum(pdfob$grCS) / sum(npCS)
    pdfob$gr  <- scCS * npCS + scC * npC + scS * npS
    
    pdfob$grC <- scC * npC 
    pdfob$grCS <- scCS * npCS
    pdfob$grS <- scS  * npS
  }
  pdfob
}
